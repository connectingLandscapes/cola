"""
SDM MODIS Export Pipeline
==========================
Exports gap-filled 16-day NDVI/EVI composites at 250m as EE assets,
ready for use in wall-to-wall SDM prediction.

All outputs are derived from SPECIES + TARGET_YEAR (model-independent),
so the same metrics can be reused across multiple model runs.

Pipeline stages (run in order):
  Stage 1 — export_annual:  build + export 23 x 16-day composites per year
  Stage 2 — gap_fill:       sequential merge + linear interpolation
  Stage 3 — reduce_to_metrics: compute p25/p50/p75/amp, export final asset

Each stage skips tiles whose output assets already exist, so the script
is fully resumable — rerun with the same arguments after any failure.

Usage:
  python sdm_modis_export.py --species puma --target_year 2025 \
      --stage export_annual --run_mode full [--ee_project geersprocessing] ...
"""

import ee, time, math, argparse
from datetime import datetime

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Export Pipeline')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--species',      type=str, required=True)
    p.add_argument('--target_year',  type=int, required=True)

    # ── Pipeline control ──────────────────────────────────────────────────────
    p.add_argument('--stage',        type=str, required=True,
                   choices=['export_annual', 'gap_fill', 'reduce_to_metrics'])
    p.add_argument('--run_mode',     type=str, default='full',
                   choices=['test', 'full'])
    p.add_argument('--max_concurrent', type=int, default=3,
                   help='Max concurrent GEE batch tasks. Leave headroom '
                        'for other work (default 3).')

    # ── Temporal / spatial (shared across scripts) ────────────────────────────
    p.add_argument('--gap_years',    type=int, default=2)
    p.add_argument('--min_year',     type=int, default=2000)
    p.add_argument('--max_year',     type=int, default=2025)
    p.add_argument('--crs',          type=str, default='EPSG:4326')
    p.add_argument('--scale',        type=int, default=250)
    p.add_argument('--tile_degrees', type=int, default=2)

    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--gee_assets',   type=str,
                   default='projects/geersprocessing/assets')
    p.add_argument('--range_asset',  type=str, default=None)

    return p.parse_args()

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

SPECIES      = args.species
TARGET_YEAR  = args.target_year
STAGE        = args.stage
RUN_MODE     = args.run_mode
MAX_CONCURRENT = args.max_concurrent
GAP_YEARS    = args.gap_years
MIN_YEAR     = args.min_year
MAX_YEAR     = args.max_year
CRS          = args.crs
SCALE        = args.scale
TILE_DEGREES = args.tile_degrees
GEE_ASSETS   = args.gee_assets

RANGE_ASSET    = args.range_asset or f'{GEE_ASSETS}/{SPECIES}_iucn_range'
BASE_FOLDER    = f'{GEE_ASSETS}/{SPECIES}_SDM_modis_{TARGET_YEAR}'
TEMP_FOLDER    = f'{BASE_FOLDER}/temp'
METRICS_FOLDER = f'{BASE_FOLDER}/metrics'

# =============================================================================
# 2. FIXED PARAMETERS
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
N_INTERVALS = len(INTERVALS)
BAND_NAMES  = [f'{idx}_{iv[2]}' for iv in INTERVALS for idx in IDX_BANDS]

# =============================================================================
# 3. DESCRIPTION SANITIZER
# GEE task descriptions allow only: a-z A-Z 0-9 . , : ; _ -
# Tile names use '+' for positive coordinates — must be stripped.
# Asset IDs are unaffected (paths allow '+').
# =============================================================================

def safe_name(s):
    """Replace '+' with 'p' in strings used as GEE task descriptions or
    asset IDs. GEE disallows '+' in both contexts.
    e.g. tile_-082_+00 -> tile_-082_p00"""
    return s.replace('+', 'p')

# =============================================================================
# 4. TASK MANAGEMENT
# =============================================================================

def submit_and_wait(task_list, max_concurrent=3, check_interval=30):
    total_eecu = 0
    completed, failed = [], []
    queue, running    = list(task_list), []
    while queue or running:
        while len(running) < max_concurrent and queue:
            item = queue.pop(0)
            item['task'].start()
            running.append(item)
            print(f'  Started: {item["desc"]} ({len(running)}/{max_concurrent})')
        still = []
        for item in running:
            s = item['task'].status()
            if s['state'] == 'COMPLETED':
                eecu = s.get('batch_eecu_usage_seconds', 0)
                total_eecu += eecu
                print(f'  OK {item["desc"]} -- {eecu:.1f} EECU-s')
                completed.append(item)
            elif s['state'] == 'FAILED':
                err = s.get('error_message', '')
                print(f'  FAIL {item["desc"]} -- {err}')
                if not item.get('retried'):
                    item['retried'] = True
                    try:
                        item['task'].start()
                        still.append(item)
                        print(f'  Retrying...')
                    except Exception as e:
                        print(f'  Retry error: {e}')
                        failed.append(item)
                else:
                    failed.append(item)
            else:
                still.append(item)
        running = still
        if running or queue:
            time.sleep(check_interval)
    return completed, failed, total_eecu

# =============================================================================
# 5. MODIS COLLECTION BUILDER
# =============================================================================

def get_modis_collection(year, aoi):
    start = ee.Date.fromYMD(year, 1, 1)
    end   = ee.Date.fromYMD(year, 12, 31)

    def mask_and_scale(img):
        qa   = img.select('SummaryQA')
        good = qa.eq(0)
        ndvi = img.select('NDVI').multiply(1).updateMask(good).rename('NDVI')
        evi  = img.select('EVI').multiply(1).updateMask(good).rename('EVI')
        return (ndvi.addBands(evi)
                    .copyProperties(img, ['system:time_start', 'system:index']))

    terra = (ee.ImageCollection('MODIS/061/MOD13Q1')
               .filterDate(start, end).filterBounds(aoi).map(mask_and_scale))
    aqua  = (ee.ImageCollection('MODIS/061/MYD13A1')
               .filterDate(start, end).filterBounds(aoi).map(mask_and_scale))
    return terra.merge(aqua)

# =============================================================================
# 6. BUILD ANNUAL 16-DAY STACK
# =============================================================================

def build_annual_modis(year, aoi):
    col   = get_modis_collection(year, aoi)
    bands = []
    for doy_start, doy_end, label in INTERVALS:
        interval_col = col.filter(
            ee.Filter.calendarRange(doy_start, doy_end, 'day_of_year'))
        for idx in IDX_BANDS:
            final_name = f'{idx}_{label}'
            composite  = ee.Image(ee.Algorithms.If(
                interval_col.size().gt(0),
                interval_col.select(idx).median().rename(final_name).toInt16(),
                ee.Image.constant(0).rename(final_name).toInt16()
                            .updateMask(ee.Image(0))
            ))
            bands.append(composite)
    img = bands[0]
    for b in bands[1:]:
        img = img.addBands(b)
    return img.toInt16()

# =============================================================================
# 7. STAGE 1: EXPORT ANNUAL STACKS
# =============================================================================

def run_export_annual(tiles, max_concurrent=3):
    years     = list(range(TARGET_YEAR, TARGET_YEAR - GAP_YEARS - 1, -1))
    task_list = []

    for tile in tiles:
        for year in years:
            yr       = max(MIN_YEAR, min(MAX_YEAR, year))
            asset_id = _annual_asset_id(tile['name'], yr)

            if asset_exists(asset_id):
                print(f'  Skipping (exists): {asset_id.split("/")[-1]}')
                continue

            img  = build_annual_modis(yr, tile['geom'])
            desc = safe_name(f'{SPECIES}_modis_{yr}_{tile["name"]}')
            task = ee.batch.Export.image.toAsset(
                image=img, description=desc, assetId=asset_id,
                region=tile['geom'], scale=SCALE, crs=CRS, maxPixels=1e11)
            task_list.append({'task': task, 'desc': desc,
                              'tile': tile['name'], 'year': yr})

    if not task_list:
        print('  All annual assets already exist.')
        return []

    print(f'Submitting {len(task_list)} annual export tasks '
          f'({len(tiles)} tiles x {GAP_YEARS+1} years)...')
    completed, failed, total_eecu = submit_and_wait(task_list, max_concurrent)
    print(f'  Completed: {len(completed)}, Failed: {len(failed)}, '
          f'EECU-s: {total_eecu:,.1f}')
    return completed

# =============================================================================
# 8. STAGE 2: GAP-FILL (barrier-synchronized parallel)
#
# All tiles advance through each step together:
#   Barrier 1: merge primary year       (all tiles)
#   Barrier 2..N: merge each prior year (all tiles)
#   Barrier N+1: export prefinal        (all tiles)
#   Barrier N+2: export gapfilled       (all tiles, includes linear interp)
#
# Resumable at each barrier — if an intermediate asset already exists for a
# tile it is skipped at that step. Tiles that fail at any barrier are dropped
# from subsequent steps and reported at the end.
# =============================================================================

def _build_interp_image(prefinal_id):
    """Build linear-interpolated gap-filled image from a prefinal asset."""
    prefinal     = ee.Image(prefinal_id)
    n_bands      = len(BAND_NAMES)
    bands_interp = []

    for i, bname in enumerate(BAND_NAMES):
        band      = prefinal.select([bname])
        still_gap = band.mask().Not()
        pv, pd = None, 0
        for p in range(i-1, max(-1, i-4), -1):
            if pv is None: pv, pd = prefinal.select([BAND_NAMES[p]]), i - p
        nv, nd = None, 0
        for n in range(i+1, min(n_bands, i+4)):
            if nv is None: nv, nd = prefinal.select([BAND_NAMES[n]]), n - i
        if pv is not None and nv is not None:
            t = pd + nd
            interp = pv.multiply(nd/t).add(nv.multiply(pd/t)).toInt16()
        elif pv is not None: interp = pv.toInt16()
        elif nv is not None: interp = nv.toInt16()
        else:                interp = band
        bands_interp.append(
            band.unmask(interp.updateMask(still_gap)).rename(bname))

    interp_stack = bands_interp[0]
    for b in bands_interp[1:]:
        interp_stack = interp_stack.addBands(b)
    was_gap     = prefinal.select(BAND_NAMES[0]).mask().Not()
    now_filled  = interp_stack.select(BAND_NAMES[0]).mask()
    interp_qual = prefinal.select('source_year').where(
        was_gap.And(now_filled), ee.Image.constant(-99).toInt16())
    return interp_stack.addBands(interp_qual)


def run_gap_fill(tiles, max_concurrent=3):
    years = list(range(TARGET_YEAR, TARGET_YEAR - GAP_YEARS - 1, -1))

    # ── Filter to tiles that need processing ─────────────────────────────────
    pending = []
    skipped_done, skipped_missing = 0, 0
    for tile in tiles:
        if asset_exists(_gapfilled_asset_id(tile['name'])):
            skipped_done += 1
            continue
        missing = [y for y in years
                   if not asset_exists(_annual_asset_id(tile['name'], y))]
        if missing:
            print(f'  Skipping {tile["name"]}: missing annual assets {missing}')
            skipped_missing += 1
            continue
        pending.append(tile)

    print(f'  Tiles to gap-fill: {len(pending)} | '
          f'Already done: {skipped_done} | '
          f'Missing inputs: {skipped_missing}')

    if not pending:
        print('  Nothing to do.')
        return

    # Track current merge asset ID per tile — updated after each barrier
    current_ids = {}
    all_eecu    = 0
    failed_names = set()

    # ── Barrier 1: export primary year + source_year band ────────────────────
    print(f'\n  Barrier 1/{GAP_YEARS+2}: merge year {years[0]} '
          f'({len(pending)} tiles)...')
    task_list = []
    for tile in pending:
        start_id = _temp_merge_id(tile['name'], years[0])
        if asset_exists(start_id):
            print(f'    Resuming (exists): {tile["name"]} y{years[0]}_merged')
            current_ids[tile['name']] = start_id
            continue
        primary = ee.Image(_annual_asset_id(tile['name'], years[0]))
        qual    = (primary.select(BAND_NAMES[0]).mask()
                          .multiply(years[0]).toInt16().rename('source_year'))
        desc    = safe_name(f'{SPECIES}_modismerge_{years[0]}_{tile["name"]}')
        task    = ee.batch.Export.image.toAsset(
            image=primary.addBands(qual), description=desc,
            assetId=start_id, region=tile['geom'],
            scale=SCALE, crs=CRS, maxPixels=1e11)
        task_list.append({'task': task, 'desc': desc,
                          'tile': tile, 'out_id': start_id,
                          'del_id': _annual_asset_id(tile['name'], years[0])})

    if task_list:
        completed, failed, eecu = submit_and_wait(task_list, max_concurrent)
        all_eecu += eecu
        for item in completed:
            current_ids[item['tile']['name']] = item['out_id']
            try: ee.data.deleteAsset(item['del_id'])
            except Exception: pass
        failed_names |= {item['tile']['name'] for item in failed}

    # ── Barriers 2..N: merge each prior year ─────────────────────────────────
    for step, year in enumerate(years[1:], start=2):
        yr           = max(MIN_YEAR, min(MAX_YEAR, year))
        active       = [t for t in pending if t['name'] not in failed_names
                        and t['name'] in current_ids]
        print(f'\n  Barrier {step}/{GAP_YEARS+2}: merge year {yr} '
              f'({len(active)} tiles)...')
        task_list = []
        for tile in active:
            merged_id = _temp_merge_id(tile['name'], yr)
            if asset_exists(merged_id):
                print(f'    Resuming (exists): {tile["name"]} y{yr}_merged')
                current_ids[tile['name']] = merged_id
                continue
            current_id = current_ids[tile['name']]
            prior_id   = _annual_asset_id(tile['name'], yr)
            current    = ee.Image(current_id)
            prior      = ee.Image(prior_id)
            still_gap  = current.select(BAND_NAMES[0]).mask().Not()
            filled     = current.select(BAND_NAMES).unmask(prior.select(BAND_NAMES))
            new_qual   = current.select('source_year').where(
                still_gap.And(prior.select(BAND_NAMES[0]).mask()),
                ee.Image.constant(yr).toInt16())
            merged     = filled.addBands(new_qual)
            desc       = safe_name(f'{SPECIES}_modismerge_{yr}_{tile["name"]}')
            task       = ee.batch.Export.image.toAsset(
                image=merged, description=desc, assetId=merged_id,
                region=tile['geom'], scale=SCALE, crs=CRS, maxPixels=1e11)
            task_list.append({'task': task, 'desc': desc,
                              'tile': tile, 'out_id': merged_id,
                              'old_current': current_id, 'old_prior': prior_id})

        if task_list:
            completed, failed, eecu = submit_and_wait(task_list, max_concurrent)
            all_eecu += eecu
            for item in completed:
                current_ids[item['tile']['name']] = item['out_id']
                for old_id in [item['old_current'], item['old_prior']]:
                    try: ee.data.deleteAsset(old_id)
                    except Exception: pass
            failed_names |= {item['tile']['name'] for item in failed}

    # ── Barrier N+1: export prefinal ─────────────────────────────────────────
    active = [t for t in pending if t['name'] not in failed_names
              and t['name'] in current_ids]
    print(f'\n  Barrier {GAP_YEARS+2}/{GAP_YEARS+2}: export prefinal '
          f'({len(active)} tiles)...')
    task_list = []
    for tile in active:
        prefinal_id = _prefinal_asset_id(tile['name'])
        if asset_exists(prefinal_id):
            print(f'    Resuming (exists): {tile["name"]} prefinal')
            continue
        current_id = current_ids[tile['name']]
        desc       = safe_name(f'{SPECIES}_modisprefinal_{tile["name"]}')
        task       = ee.batch.Export.image.toAsset(
            image=ee.Image(current_id), description=desc,
            assetId=prefinal_id, region=tile['geom'],
            scale=SCALE, crs=CRS, maxPixels=1e11)
        task_list.append({'task': task, 'desc': desc,
                          'tile': tile, 'out_id': prefinal_id,
                          'del_id': current_id})

    if task_list:
        completed, failed, eecu = submit_and_wait(task_list, max_concurrent)
        all_eecu += eecu
        for item in completed:
            try: ee.data.deleteAsset(item['del_id'])
            except Exception: pass
        failed_names |= {item['tile']['name'] for item in failed}

    # ── Barrier N+2: linear interpolation + export gapfilled ─────────────────
    active = [t for t in pending if t['name'] not in failed_names
              and asset_exists(_prefinal_asset_id(t['name']))]
    print(f'\n  Barrier {GAP_YEARS+3}/{GAP_YEARS+2}+1: '
          f'interp + gapfilled export ({len(active)} tiles)...')
    task_list = []
    for tile in active:
        final_id    = _gapfilled_asset_id(tile['name'])
        prefinal_id = _prefinal_asset_id(tile['name'])
        if asset_exists(final_id):
            print(f'    Resuming (exists): {tile["name"]} gapfilled')
            try: ee.data.deleteAsset(prefinal_id)
            except Exception: pass
            continue
        final_img = _build_interp_image(prefinal_id)
        desc      = safe_name(f'{SPECIES}_modisgapfill_{tile["name"]}')
        task      = ee.batch.Export.image.toAsset(
            image=final_img, description=desc, assetId=final_id,
            region=tile['geom'], scale=SCALE, crs=CRS, maxPixels=1e11)
        task_list.append({'task': task, 'desc': desc,
                          'tile': tile, 'prefinal_id': prefinal_id})

    if task_list:
        completed, failed, eecu = submit_and_wait(task_list, max_concurrent)
        all_eecu += eecu
        for item in completed:
            try: ee.data.deleteAsset(item['prefinal_id'])
            except Exception: pass
        failed_names |= {item['tile']['name'] for item in failed}

    # ── Summary ───────────────────────────────────────────────────────────────
    n_done   = len(pending) - len(failed_names)
    print(f'\n  Gap-fill complete: {n_done}/{len(pending)} tiles succeeded, '
          f'{all_eecu:,.1f} EECU-s')
    if failed_names:
        print(f'  ⚠️  {len(failed_names)} tiles failed — rerun to retry:')
        for name in sorted(failed_names):
            print(f'    {name}')

# =============================================================================
# 9. STAGE 3: COMPUTE METRICS
# =============================================================================

def run_reduce_to_metrics(tiles, max_concurrent=3):
    task_list = []

    for tile in tiles:
        metrics_id = f'{METRICS_FOLDER}/{tile["name"]}_metrics'

        if asset_exists(metrics_id):
            print(f'  Skipping (exists): {tile["name"]}')
            continue

        gapfill_id = _gapfilled_asset_id(tile['name'])
        if not asset_exists(gapfill_id):
            print(f'  Skipping {tile["name"]}: run gap_fill stage first')
            continue

        img = ee.Image(gapfill_id)
        all_metric_bands = []
        for idx in IDX_BANDS:
            idx_col = ee.ImageCollection(
                [img.select([f'{idx}_{iv[2]}']).rename('IDX') for iv in INTERVALS])
            p25 = idx_col.reduce(ee.Reducer.percentile([25])).rename(f'{idx}_p25')
            p50 = idx_col.median()                            .rename(f'{idx}_p50')
            p75 = idx_col.reduce(ee.Reducer.percentile([75])).rename(f'{idx}_p75')
            amp = p75.subtract(p25)                           .rename(f'{idx}_amp')
            all_metric_bands.extend([p25, p50, p75, amp])

        qual = img.select('source_year').rename('source_year')
        metrics = all_metric_bands[0]
        for b in all_metric_bands[1:]:
            metrics = metrics.addBands(b)
        metrics = metrics.addBands(qual)

        desc = safe_name(f'{SPECIES}_modismetrics_{tile["name"]}')
        task = ee.batch.Export.image.toAsset(
            image=metrics, description=desc, assetId=metrics_id,
            region=tile['geom'], scale=SCALE, crs=CRS, maxPixels=1e10)
        task_list.append({'task': task, 'desc': desc, 'tile': tile['name']})

    if not task_list:
        print('  All metric assets already exist.')
        return []

    print(f'Submitting {len(task_list)} metric tasks...')
    completed, failed, total_eecu = submit_and_wait(task_list, max_concurrent)
    print(f'  Completed: {len(completed)}, Failed: {len(failed)}, '
          f'EECU-s: {total_eecu:,.1f}')

    for item in completed:
        try:
            ee.data.deleteAsset(_gapfilled_asset_id(item['tile']))
            print(f'  Deleted gapfilled: {item["tile"]}')
        except Exception: pass

    return completed

# =============================================================================
# 10. ASSET PATH HELPERS
# =============================================================================

def _annual_asset_id(tile_name, year):
    return f'{TEMP_FOLDER}/{tile_name}_y{year}_annual'

def _temp_merge_id(tile_name, year):
    return f'{TEMP_FOLDER}/{tile_name}_y{year}_merged'

def _prefinal_asset_id(tile_name):
    return f'{TEMP_FOLDER}/{tile_name}_prefinal'

def _gapfilled_asset_id(tile_name):
    return f'{TEMP_FOLDER}/{tile_name}_gapfilled'

def asset_exists(asset_id):
    try:
        ee.data.getAsset(asset_id)
        return True
    except Exception:
        return False

def ensure_folders():
    for folder in [BASE_FOLDER, TEMP_FOLDER, METRICS_FOLDER]:
        try:
            ee.data.createAsset({'type': 'Folder'}, folder)
            print(f'  Created: {folder}')
        except Exception as e:
            e_str = str(e).lower()
            if 'already exists' in e_str or 'cannot overwrite' in e_str:
                print(f'  Exists:  {folder}')
            else:
                print(f'  Error:   {folder}: {e}')

def get_test_tiles():
    d = TILE_DEGREES
    return [{'geom': ee.Geometry.Rectangle([-84, 9, -84+d, 9+d]),
             'name': f'test_costa_rica_{d}deg'}]

def build_tiles(range_fc):
    bounds  = range_fc.geometry().bounds()
    coords  = bounds.coordinates().get(0).getInfo()
    lons    = [c[0] for c in coords]
    lats    = [c[1] for c in coords]
    min_lon = math.floor(min(lons) / TILE_DEGREES) * TILE_DEGREES
    min_lat = math.floor(min(lats) / TILE_DEGREES) * TILE_DEGREES
    max_lon = math.ceil(max(lons)  / TILE_DEGREES) * TILE_DEGREES
    max_lat = math.ceil(max(lats)  / TILE_DEGREES) * TILE_DEGREES
    tiles   = []
    lat = min_lat
    while lat < max_lat:
        lon = min_lon
        while lon < max_lon:
            geom = ee.Geometry.Rectangle(
                [lon, lat, lon+TILE_DEGREES, lat+TILE_DEGREES])
            if range_fc.geometry().intersects(
                    geom, ee.ErrorMargin(1000)).getInfo():
                tiles.append({'geom': geom,
                              'name': safe_name(f'tile_{int(lon):+04d}_{int(lat):+03d}')})
            lon += TILE_DEGREES
        lat += TILE_DEGREES
    print(f'Built {len(tiles)} tiles ({TILE_DEGREES} deg) over range')
    return tiles

def next_stage():
    stages = ['export_annual', 'gap_fill', 'reduce_to_metrics']
    idx    = stages.index(STAGE) if STAGE in stages else 0
    return stages[idx+1] if idx < len(stages)-1 else 'done'

# =============================================================================
# 11. MAIN
# =============================================================================

def main():
    print('=' * 60)
    print(f'SDM MODIS EXPORT')
    print(f'Species: {SPECIES} | Year: {TARGET_YEAR}')
    print(f'Stage: {STAGE} | Mode: {RUN_MODE} | Max concurrent: {MAX_CONCURRENT}')
    print(f'Base folder:    {BASE_FOLDER}')
    print(f'Metrics folder: {METRICS_FOLDER}')
    print('=' * 60)

    ensure_folders()

    if RUN_MODE == 'test':
        tiles = get_test_tiles()
        print(f'\nTest tile: {tiles[0]["name"]}')
    else:
        range_fc = ee.FeatureCollection(RANGE_ASSET)
        tiles    = build_tiles(range_fc)

    print()
    if STAGE == 'export_annual':
        run_export_annual(tiles, max_concurrent=MAX_CONCURRENT)
    elif STAGE == 'gap_fill':
        run_gap_fill(tiles, max_concurrent=MAX_CONCURRENT)
    elif STAGE == 'reduce_to_metrics':
        run_reduce_to_metrics(tiles, max_concurrent=MAX_CONCURRENT)

    print(f'\nStage complete. Next: --stage {next_stage()}')

if __name__ == '__main__':
    main()
