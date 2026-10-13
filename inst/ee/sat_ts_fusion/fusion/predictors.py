"""
Shared predictor construction for sdm_extraction.py (training),
sdm_wall_to_wall.py (prediction) and the optional sdm_modis_export.py
(MODIS metric tiles).

Training and prediction stacks are both built by assemble_stack(), so the
predictors cannot drift apart between scripts (they did before 2026-10-06:
see CLAUDE.md, "Train/predict consistency fixes").

Nothing here calls the Earth Engine API at import time, so the pure-Python
helpers can be unit-tested without ee.Initialize().
"""
import math
import re
import ee

CRS = 'EPSG:4326'
SENTINEL = -9999

# =============================================================================
# SOURCES
# =============================================================================

# native_scale: default --scale (m). Users may choose a coarser scale, never
#   a finer one.
# coverage: 'jan_dec' = composites need imagery in January and December of a
#   year; 'annual' = one image per year (time_start 1 Jan).
# exportable: True if metric tiles can be pre-exported (sdm_modis_export.py)
#   for wall-to-wall --use_exported_metrics. Wall-to-wall computes metrics in
#   memory by default for every source: for MODIS that costs 221 EECU-s per
#   2 deg test tile vs 254 export + 74-87 prediction, identical values
#   (2026-10-09), so export pays off only when 2+ models map the same
#   species and year.
# Coverage years are read from the collections, not listed here.
SOURCES = {
    'modis': {'native_scale': 250, 'coverage': 'jan_dec', 'exportable': True,
              'collections': ['MODIS/061/MOD13Q1', 'MODIS/061/MYD13Q1']},
    # Google Satellite Embedding V1 (AlphaEarth Foundations): 64 unit-length
    # float bands A00-A63, 10 m UTM tiles (~164 km), 2017-2025 as of
    # 2026-10-07, masked over open ocean.
    # HLS v2.0: 30 m, computed in memory for training and prediction (no
    # export stage). Coarse kernel grids are computed directly at the grid
    # scale (benchmark 2026-10-07: 38x cheaper than aggregating 30 m, r 0.96).
    'hls_l30': {'native_scale': 30, 'coverage': 'jan_dec', 'exportable': False,
                'coarse_from': 'direct', 'nir_band': 'B5',
                'collections': ['NASA/HLS/HLSL30/v002']},
    'hls_s30': {'native_scale': 30, 'coverage': 'jan_dec', 'exportable': False,
                'coarse_from': 'direct', 'nir_band': 'B8A',
                'collections': ['NASA/HLS/HLSS30/v002']},
    'embed': {'native_scale': 10, 'coverage': 'annual', 'exportable': False,
              'coarse_from': 'direct',
              'collections': ['GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL']},
}
EMBED_BANDS = [f'A{i:02d}' for i in range(64)]


def resolve_scale(source, scale):
    native = SOURCES[source]['native_scale']
    if scale is None:
        return native
    if scale < native:
        raise ValueError(f'--scale {scale} is finer than {source} native '
                         f'resolution ({native} m)')
    return int(scale)


def source_tag(source, scale):
    """Used in asset folder names, e.g. 'modis_250m'."""
    return f'{source}_{scale}m'


# =============================================================================
# COVERAGE
# 'jan_dec' sources: a year is usable only if there is imagery in both
# January and December (full-year composites). 'annual' sources: a year is
# usable if its annual image exists. Occurrence years outside coverage stop
# the run (decided 2026-10-07: fail loudly, no silent clamping).
# =============================================================================

def covered_years(flags):
    """flags: {year: (start_ok, end_ok)} -> sorted usable years."""
    return sorted(y for y, (jan, dec) in flags.items() if jan and dec)


def coverage_years(source, region=None, first=1999, last=None):
    """Usable years for `source`, from one server call."""
    import datetime
    last = last or datetime.date.today().year
    col = source_collection(source, region)
    annual = SOURCES[source]['coverage'] == 'annual'

    def probe(y):
        y = ee.Number(y)
        if annual:
            any_ = col.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y.add(1), 1, 1)).limit(1).size().gt(0)
            return ee.List([y, any_, any_])
        # limit(1): HLS has millions of images; existence is all we need
        jan = col.filterDate(ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 2, 1)).limit(1).size().gt(0)
        dec = col.filterDate(ee.Date.fromYMD(y, 12, 1), ee.Date.fromYMD(y.add(1), 1, 1)).limit(1).size().gt(0)
        return ee.List([y, jan, dec])

    rows = ee.List.sequence(first, last).map(probe).getInfo()
    return covered_years({int(y): (bool(j), bool(d)) for y, j, d in rows})


def whole_years(values):
    """Occurrence `year` values as ints. Assets uploaded from CSV or
    shapefiles can store years as doubles (2019.0), which would put a '.'
    in batch asset IDs. Raises ValueError for missing or non-integer
    years instead of rounding them."""
    bad = sorted({repr(v) for v in values
                  if not (isinstance(v, int) or
                          (isinstance(v, float) and v.is_integer()))})
    if bad:
        raise ValueError(f'Occurrence year values must be whole numbers; '
                         f'found {", ".join(bad[:10])}')
    return [int(v) for v in values]


def out_of_coverage(years, usable):
    """{year: count} for years not in `usable`. `years` is an iterable of
    occurrence years (one entry per point)."""
    usable = set(usable)
    bad = {}
    for y in years:
        if y not in usable:
            bad[y] = bad.get(y, 0) + 1
    return dict(sorted(bad.items()))


def coverage_message(source, usable, bad):
    span = f'{usable[0]}-{usable[-1]}' if usable else 'none'
    lines = [f'{sum(bad.values())} occurrence(s) fall outside {source} '
             f'coverage (full years {span}):']
    lines += [f'  {y}: {n}' for y, n in bad.items()]
    lines.append('Choose a source that covers these years, or remove these '
                 'points from the occurrence asset.')
    return '\n'.join(lines)


# =============================================================================
# NEIGHBOURHOOD SCALES
# Base scale: the analysis --scale, smoothed with a Gaussian of radius
# scale/2 (the original 250 m scale used 125 m). Larger scales: Gaussian
# radii below, sigma = radius/3, each computed on a grid of radius /
# KERNEL_GRID_FRACTION (never finer than the base) and resampled bilinear.
# Benchmark 2026-10-07, 10 km kernel on a 300 m base: 56x cheaper than at
# the base scale, p95 |diff| 0.0013 NDVI (~1% of across-tile SD).
# =============================================================================

KERNEL_RADII_M       = [1000, 5000, 10000]
KERNEL_GRID_FRACTION = 10   # provisional; 10 tested only for the 10 km kernel

# Where coarse grids are averaged from (`coarse_from`). reduceResolution
# makes every grid pixel read all fine pixels beneath it, so averaging from
# a 10 m base is ruinous: a 1 km grid pixel reads 10,000 pixels per band and
# one extraction batch of 42 embedding points cost 36,174 EECU-s
# (2026-10-07). Options:
#   None    - average from the base scale (original behaviour)
#   metres  - average from max(base, metres): the layer's own native scale
#   'direct'- read the image at the grid scale (GEE's pyramid). Only for a
#             source whose asset pyramids by mean and that training and
#             prediction both read in memory (not slope, which is computed
#             from the DEM at the requested scale, nor LCLUC fractions,
#             which would come from a mode-pyramided class map).
# At a 250 m base every option except 'direct' reduces to None, so MODIS
# models are unchanged.
STATIC_NATIVE_M = 30    # SRTM, LCLUC; WorldClim/CHILI/mTPI/GHSL are coarser


def scale_suffix(metres):
    m = int(round(metres))
    return f'_{m // 1000}km' if m >= 1000 and m % 1000 == 0 else f'_{m}m'


def kernel_scales(base_scale):
    """[(suffix, radius_m, grid_m)], base scale first."""
    out = [(scale_suffix(base_scale), base_scale / 2, base_scale)]
    for r in KERNEL_RADII_M:
        if r > base_scale:
            out.append((scale_suffix(r), r, max(base_scale, r / KERNEL_GRID_FRACTION)))
    return out


SCALE_SUFFIX_RE = re.compile(r'_(\d+m|\d+km)$')


def base_name(band):
    """'NDVI_p50_10km' -> 'NDVI_p50'."""
    return SCALE_SUFFIX_RE.sub('', band)


# Masking: the base scale keeps the pixel's own validity (masked where the
# pixel is masked). Larger scales are defined wherever the neighbourhood has
# data (skipMasked=False), e.g. mean tree height within 10 km at a non-tree
# pixel. Before 2026-10-07 every scale was masked by the centre pixel, so
# lcluc_*_ht at 1/5/10 km were -9999 (filled with 0) at any non-tree point.
def add_neighborhood_scales(image, base_scale, crs=CRS, coarse_from=None):
    base  = image.reproject(crs=crs, scale=base_scale)
    fine  = base_scale if coarse_from in (None, 'direct') else max(base_scale, coarse_from)
    fine_img = base if fine == base_scale else image.reproject(crs=crs, scale=fine)
    names = image.bandNames()
    parts = []
    for sfx, radius, grid in kernel_scales(base_scale):
        kernel = ee.Kernel.gaussian(radius, radius / 3, 'meters')
        if radius < base_scale:      # the base scale itself
            smooth = base.reduceNeighborhood(ee.Reducer.mean(), kernel)
        elif grid == base_scale:
            smooth = base.reduceNeighborhood(ee.Reducer.mean(), kernel, skipMasked=False)
        else:
            if coarse_from == 'direct':
                coarse = image.reproject(crs=crs, scale=grid)
            else:
                coarse = (fine_img.reduceResolution(ee.Reducer.mean(),
                                                    maxPixels=math.ceil((grid / fine) ** 2) + 1)
                                  .reproject(crs=crs, scale=grid))
            smooth = (coarse.reduceNeighborhood(ee.Reducer.mean(), kernel, skipMasked=False)
                            .reproject(crs=crs, scale=grid).resample('bilinear'))
        parts.append(smooth.rename(names.map(lambda b: ee.String(b).cat(sfx))))
    return ee.Image.cat(parts)


# =============================================================================
# STATIC PREDICTORS (all sources)
# Continuous layers are resampled bilinear so coarse inputs (WorldClim
# ~1 km) do not show block edges at fine scales (decided 2026-10-07).
# LCLUC is categorical-derived and keeps nearest-neighbour.
# =============================================================================

CONTINUOUS_RESAMPLE = 'bilinear'   # None = nearest (pre-2026-10-07 behaviour)
EPOCH_YEARS = [2000, 2005, 2010, 2015, 2020]


def nearest_epoch(year):
    """LCLUC/GHSL epoch used for an occurrence/target year."""
    if year < 2003: return 2000
    if year < 2008: return 2005
    if year < 2013: return 2010
    if year < 2018: return 2015
    return 2020


def lcluc_bands(lcluc):
    is_terra_short = lcluc.gte(0).And(lcluc.lte(24))
    is_terra_trees = lcluc.gte(25).And(lcluc.lte(48))
    is_wet_short   = lcluc.gte(100).And(lcluc.lte(124))
    is_wet_trees   = lcluc.gte(125).And(lcluc.lte(148))
    is_water       = lcluc.gte(200).And(lcluc.lte(207))
    is_cropland    = lcluc.eq(244)
    is_builtup     = lcluc.eq(250)
    return ee.Image.cat([
        is_terra_short.rename('lcluc_short_veg').toFloat(),
        is_terra_trees.rename('lcluc_trees').toFloat(),
        is_wet_short.rename('lcluc_wet_short').toFloat(),
        is_wet_trees.rename('lcluc_wet_trees').toFloat(),
        is_water.rename('lcluc_water').toFloat(),
        is_cropland.rename('lcluc_cropland').toFloat(),
        is_builtup.rename('lcluc_builtup').toFloat(),
        lcluc.multiply(4).add(3).updateMask(is_terra_short).toFloat().rename('lcluc_terra_short_pct'),
        lcluc.subtract(100).multiply(4).add(3).updateMask(is_wet_short).toFloat().rename('lcluc_wet_short_pct'),
        lcluc.subtract(22).updateMask(is_terra_trees).toFloat().rename('lcluc_terra_ht'),
        lcluc.subtract(122).updateMask(is_wet_trees).toFloat().rename('lcluc_wet_ht'),
    ])


def static_layers(year):
    """(continuous, categorical) static predictor images for `year`."""
    ep   = nearest_epoch(year)
    srtm = ee.Image('USGS/SRTMGL1_003')
    cont = ee.Image.cat([
        ee.Image('WORLDCLIM/V1/BIO'),
        srtm.select('elevation').rename('elevation'),
        ee.Terrain.slope(srtm).rename('slope'),
        ee.Image('CSP/ERGo/1_0/Global/SRTM_CHILI').rename('chili'),
        ee.Image('CSP/ERGo/1_0/Global/SRTM_mTPI').rename('mtpi'),
        ee.Image(f'JRC/GHSL/P2023A/GHS_POP/{ep}').select([0])
          .rename('pop_density').max(ee.Image(0)),
    ]).toFloat()
    if CONTINUOUS_RESAMPLE:
        cont = cont.resample(CONTINUOUS_RESAMPLE)
    cat = lcluc_bands(ee.Image(f'projects/glad/GLCLU2020/v2/LCLUC_{ep}'))
    return cont, cat


def assemble_stack(source, source_metrics, year, scale, crs=CRS):
    """Full predictor stack. `source_metrics` is the imagery-source image
    (MODIS metrics computed in memory for training, or read from exported
    metric tiles for prediction). Masked pixels are left masked."""
    cont, cat = static_layers(year)
    src_from = SOURCES[source].get('coarse_from')
    return ee.Image.cat([
        add_neighborhood_scales(cont, scale, crs, STATIC_NATIVE_M),
        add_neighborhood_scales(source_metrics, scale, crs, src_from),
        add_neighborhood_scales(cat, scale, crs, STATIC_NATIVE_M),
    ]).toFloat()


# =============================================================================
# VEGETATION INDEX SOURCES: MODIS 23 x 16-day and HLS monthly NDVI/EVI
# composites, shared gap-fill and metrics
# =============================================================================

IDX_BANDS = ['NDVI', 'EVI']
INTERVALS = [
    (1,16,'d001'), (17,32,'d017'), (33,48,'d033'), (49,64,'d049'),
    (65,80,'d065'), (81,96,'d081'), (97,112,'d097'), (113,128,'d113'),
    (129,144,'d129'), (145,160,'d145'), (161,176,'d161'), (177,192,'d177'),
    (193,208,'d193'), (209,224,'d209'), (225,240,'d225'), (241,256,'d241'),
    (257,272,'d257'), (273,288,'d273'), (289,304,'d289'), (305,320,'d305'),
    (321,336,'d321'), (337,352,'d337'), (353,366,'d353'),
]
N_INTERVALS  = len(INTERVALS)
LABELS       = [iv[2] for iv in INTERVALS]
METRIC_BANDS = [f'{idx}_{m}' for idx in IDX_BANDS for m in ['p25', 'p50', 'p75', 'amp']]

# HLS: calendar-month composites (decided 2026-10-07)
MONTH_LABELS = [f'm{m:02d}' for m in range(1, 13)]


def band_names(labels):
    """Interleaved NDVI/EVI band names for a list of interval labels."""
    return [f'{idx}_{lab}' for lab in labels for idx in IDX_BANDS]


BAND_NAMES = band_names(LABELS)


def _blank_composite(names):
    return ee.Image.constant([0, 0]).rename(names).toInt16().updateMask(0)


def source_collection(source, region=None):
    """Collection for `source`: MODIS merged + QA-masked (bands NDVI, EVI);
    HLS cloud-masked NDVI/EVI; embeddings as published (bands A00-A63)."""
    if source == 'embed':
        col = ee.ImageCollection(SOURCES['embed']['collections'][0])
        return col.filterBounds(region) if region is not None else col
    if source in HLS_SOURCES:
        return hls_collection(source, region)
    if source != 'modis':
        raise ValueError(f'unknown source {source}')

    def mask(img):
        good = img.select('SummaryQA').eq(0)
        return (img.select(IDX_BANDS).updateMask(good)
                   .copyProperties(img, ['system:time_start']))

    cols = [ee.ImageCollection(c) for c in SOURCES['modis']['collections']]
    if region is not None:
        cols = [c.filterBounds(region) for c in cols]
    return cols[0].merge(cols[1]).map(mask)


def modis_interval_composites(year, region=None):
    """List of 23 images, each with bands NDVI_dXXX, EVI_dXXX (int16,
    masked where the interval has no good observation)."""
    col = source_collection('modis', region).filterDate(
        ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 12, 31))
    out = []
    for doy_start, doy_end, label in INTERVALS:
        ic    = col.filter(ee.Filter.calendarRange(doy_start, doy_end, 'day_of_year'))
        names = [f'{idx}_{label}' for idx in IDX_BANDS]
        out.append(ee.Image(ee.Algorithms.If(
            ic.size().gt(0),
            ic.median().rename(names).toInt16(),
            _blank_composite(names))))
    return out


def modis_annual_image(year, region=None):
    return ee.Image.cat(modis_interval_composites(year, region)).toInt16()


# =============================================================================
# HLS: Landsat (L30) or Sentinel-2 (S30), the user's choice, never merged.
# Per image: Fmask bits 1-4 masked (cloud, adjacent cloud/shadow, shadow,
# snow/ice), NDVI/EVI scaled x10000 to int16 so the MODIS gap-fill code and
# units apply unchanged. Monthly medians; overlapping MGRS tiles duplicate
# dates, which does not change a median.
# =============================================================================

HLS_SOURCES    = ('hls_l30', 'hls_s30')
HLS_CLEAR_BITS = [1, 2, 3, 4]
# S30 footprints in UTM zones 01/60 are infinite (L30 checked: not at a
# Costa Rica point, but filtered the same way), so filterBounds anywhere
# returns them (1,743 S30 images at a Costa Rica point for 2024, 103 with
# data). Drop those zones unless the region reaches them.
HLS_ANTIMERIDIAN_ZONES = ['01', '60']
HLS_ANTIMERIDIAN_LON   = 174


def hls_masked_vi(img, nir_band):
    fm  = img.select('Fmask')
    bad = ee.Image(0)
    for b in HLS_CLEAR_BITS:
        bad = bad.Or(fm.rightShift(b).bitwiseAnd(1))
    nir, red, blue = img.select(nir_band), img.select('B4'), img.select('B2')
    # Negative surface reflectance (mostly water) gave NDVI medians down to
    # -1.13 (2026-10-08); mask it, and clamp both indices to MODIS's valid
    # range (EVI's denominator can approach 0).
    ok   = nir.gte(0).And(red.gte(0)).And(blue.gte(0))
    ndvi = nir.subtract(red).divide(nir.add(red)).clamp(-1, 1)
    evi  = img.expression('2.5 * (n - r) / (n + 6 * r - 7.5 * b + 1)',
                          {'n': nir, 'r': red, 'b': blue}).clamp(-1, 1)
    return (ee.Image.cat([ndvi, evi]).multiply(10000).round().toInt16()
              .rename(IDX_BANDS).updateMask(bad.Not().And(ok))
              .copyProperties(img, ['system:time_start']))


def hls_collection(source, region=None):
    cfg = SOURCES[source]
    col = ee.ImageCollection(cfg['collections'][0])
    if region is not None:
        col = col.filterBounds(region)
        lons = ee.List(ee.Geometry(region).bounds(1000).coordinates().get(0)).map(
            lambda c: ee.List(c).get(0))
        near = (ee.Number(lons.reduce(ee.Reducer.min())).lt(-HLS_ANTIMERIDIAN_LON)
                .Or(ee.Number(lons.reduce(ee.Reducer.max())).gt(HLS_ANTIMERIDIAN_LON)))
        # system:index is 'T{MGRS tile}_{datetime}' in both collections; L30
        # has no MGRS_TILE_ID property (filtering on it dropped every image)
        far_zone = ee.Filter.Or(*[ee.Filter.stringStartsWith('system:index', f'T{z}')
                                  for z in HLS_ANTIMERIDIAN_ZONES])
        col = ee.ImageCollection(ee.Algorithms.If(near, col, col.filter(far_zone.Not())))
    return col.map(lambda img: hls_masked_vi(img, cfg['nir_band']))


def hls_monthly_composites(source, year, region=None):
    col = source_collection(source, region).filterDate(
        ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
    out = []
    for m, label in enumerate(MONTH_LABELS, start=1):
        names = [f'{idx}_{label}' for idx in IDX_BANDS]
        # merging a fully masked image keeps the bands when the month is empty
        out.append(col.filter(ee.Filter.calendarRange(m, m, 'month'))
                      .merge(ee.ImageCollection([_blank_composite(IDX_BANDS)]))
                      .median().rename(names).toInt16())
    return out


def hls_annual_image(source, year, region=None):
    return ee.Image.cat(hls_monthly_composites(source, year, region)).toInt16()


def prior_years(year, gap_years, first_usable):
    """Prior years used for gap-fill, most recent first."""
    return [year - k for k in range(1, gap_years + 1) if year - k >= first_usable]


DEFAULT_GAP_YEARS = 2


def recorded_gap_years(values):
    """gap_years stamped by extraction on training points and carried by
    fitting onto the features asset, so prediction gap-fills the way the
    model was trained. Data from before 2026-10-09 has no stamp; it was
    extracted with the default. Returns (gap_years, stamped)."""
    found = {int(v) for v in values if v is not None and v == v}  # skip NaN
    if not found:
        return DEFAULT_GAP_YEARS, False
    if len(found) > 1:
        raise ValueError(f'mixed gap_years {sorted(found)} in training data')
    return found.pop(), True


# Gap-fill rule (identical for in-memory training and staged export, and
# shared by the MODIS 16-day and HLS monthly series):
#   1. Prior-year fill only for "long gaps": interval missing AND an adjacent
#      interval also missing (first/last interval count as their own
#      neighbour). Most recent prior year first; NDVI and EVI replaced
#      together, keyed on NDVI validity.
#   2. Remaining gaps: per-pixel linear interpolation between the nearest
#      valid intervals within max_dist on each side.
# Masked pixels must be filled with unmask(): ee.Image.where() does NOT write
# into masked pixels.
INTERP_MAX_DIST = 3  # MODIS intervals (48 days); original design intent, provisional
# HLS months: 2 (~60 days), as in the 2026-10-07 benchmark (83% of pixels
# with all 12 months after fill). Provisional.
HLS_INTERP_MAX_DIST = 2


def long_gap_bands(n):
    """Bitmask bands holding n long-gap flags, 12 per int16 band."""
    return [f'long_gap_{b}' for b in range(math.ceil(n / 12))]


LONG_GAP_BANDS = long_gap_bands(N_INTERVALS)


def _valid(img):
    return img.mask().gt(0)


def long_gap_masks(obs):
    n = len(obs)
    return [obs[i].Not().And(obs[max(0, i-1)].Not().Or(obs[min(n-1, i+1)].Not()))
            for i in range(n)]


def interp_series(series, names, max_dist=INTERP_MAX_DIST):
    """Fill gaps in a list of single-band images by per-pixel linear
    interpolation between the nearest valid neighbours <= max_dist away."""
    n, out = len(series), []
    for i in range(n):
        pv = pd = nv = nd = None
        for d in range(1, max_dist + 1):
            if i - d >= 0:
                c    = series[i-d]
                dist = ee.Image.constant(d).updateMask(_valid(c))
                pv   = c if pv is None else pv.unmask(c)
                pd   = dist if pd is None else pd.unmask(dist)
            if i + d < n:
                c    = series[i+d]
                dist = ee.Image.constant(d).updateMask(_valid(c))
                nv   = c if nv is None else nv.unmask(c)
                nd   = dist if nd is None else nd.unmask(dist)
        if pv is not None and nv is not None:
            both = (pv.toFloat().multiply(nd).add(nv.toFloat().multiply(pd))
                      .divide(pd.add(nd)))
            val  = both.unmask(pv.toFloat()).unmask(nv.toFloat())
        else:
            val  = (pv if pv is not None else nv).toFloat()
        out.append(series[i].unmask(val.round().toInt16()).toInt16().rename(names[i]))
    return out


def _long_gap_bit(img, i, lg_bands):
    return img.select(lg_bands[i // 12]).rightShift(i % 12).bitwiseAnd(1)


def first_merge_image(primary, year, labels=LABELS):
    """Target-year stack + source_year + long-gap bitmask bands (int16,
    12 intervals per band)."""
    n    = len(labels)
    lgb  = long_gap_bands(n)
    obs  = [_valid(primary.select(f'NDVI_{lab}')) for lab in labels]
    lg   = long_gap_masks(obs)
    bits = []
    for b in range(len(lgb)):
        acc = ee.Image.constant(0)
        for i in range(12*b, min(12*(b+1), n)):
            # Cast first: leftShift on a 0/1 boolean image loses the bits
            acc = acc.add(lg[i].toInt16().leftShift(i % 12))
        bits.append(acc.toInt16().rename(lgb[b]))
    qual = obs[0].multiply(year).toInt16().rename('source_year')
    return primary.select(band_names(labels)).addBands(qual).addBands(bits)


def merge_prior_image(current, prior, year, labels=LABELS):
    """Fill long gaps in `current` from one prior year's annual stack."""
    lgb = long_gap_bands(len(labels))
    parts, use0 = [], None
    for i, lab in enumerate(labels):
        names = [f'{idx}_{lab}' for idx in IDX_BANDS]
        cur   = current.select(names)
        pri   = prior.select(names)
        use   = (_long_gap_bit(current, i, lgb).And(_valid(cur.select(names[0])).Not())
                                               .And(_valid(pri.select(names[0]))))
        parts.append(cur.unmask(pri.updateMask(use)))
        if i == 0:
            use0 = use
    qual = current.select('source_year').where(use0, ee.Image.constant(year).toInt16())
    return (ee.Image.cat(parts).addBands(qual)
              .addBands(current.select(lgb)))


def interp_image(merged, labels=LABELS, max_dist=INTERP_MAX_DIST):
    """Interpolate remaining gaps; drops the long-gap bands. source_year is
    set to -99 where interval 0 was filled by interpolation."""
    names_all = band_names(labels)
    per_idx = {}
    for idx in IDX_BANDS:
        names        = [f'{idx}_{lab}' for lab in labels]
        per_idx[idx] = dict(zip(names, interp_series([merged.select(n) for n in names],
                                                     names, max_dist)))
    stack      = ee.Image.cat([per_idx[b.split('_')[0]][b] for b in names_all])
    was_gap    = _valid(merged.select(names_all[0])).Not()
    now_filled = _valid(stack.select(names_all[0]))
    qual = merged.select('source_year').where(
        was_gap.And(now_filled), ee.Image.constant(-99).toInt16())
    return stack.addBands(qual)


def modis_gapfilled_image(year, priors, region=None):
    """In-memory equivalent of export stages 1-2 (same functions, no assets)."""
    cur = first_merge_image(modis_annual_image(year, region), year)
    for y in priors:
        cur = merge_prior_image(cur, modis_annual_image(y, region), y)
    return interp_image(cur)


def hls_gapfilled_image(source, year, priors, region=None):
    """Same rule as MODIS on 12 monthly composites."""
    cur = first_merge_image(hls_annual_image(source, year, region), year, MONTH_LABELS)
    for y in priors:
        cur = merge_prior_image(cur, hls_annual_image(source, y, region), y, MONTH_LABELS)
    return interp_image(cur, MONTH_LABELS, HLS_INTERP_MAX_DIST)


def metrics_image(gapfilled, labels=LABELS):
    """p25/p50/p75/amp per index over the intervals (METRIC_BANDS)."""
    bands = []
    for idx in IDX_BANDS:
        col = ee.ImageCollection([gapfilled.select([f'{idx}_{lab}']).rename('v') for lab in labels])
        p25 = col.reduce(ee.Reducer.percentile([25])).rename(f'{idx}_p25')
        p50 = col.median().rename(f'{idx}_p50')
        p75 = col.reduce(ee.Reducer.percentile([75])).rename(f'{idx}_p75')
        bands += [p25, p50, p75, p75.subtract(p25).rename(f'{idx}_amp')]
    return ee.Image.cat(bands)


def source_predictors(source, year, priors, region=None):
    """Imagery-source predictors computed in memory: for training (all
    sources) and for prediction unless --use_exported_metrics."""
    if source == 'modis':
        return metrics_image(modis_gapfilled_image(year, priors, region))
    if source in HLS_SOURCES:
        return metrics_image(hls_gapfilled_image(source, year, priors, region), MONTH_LABELS)
    if source == 'embed':
        return (source_collection('embed', region)
                  .filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year + 1, 1, 1))
                  .mosaic().select(EMBED_BANDS))
    raise ValueError(f'unknown source {source}')
