"""
SDM MODIS Predictor Extraction
================================
Extracts predictor variables at species occurrence locations using
MODIS 250m NDVI/EVI. All outputs are named after SPECIES and MODEL_ID
so multiple runs never overwrite each other.

Usage:
  python sdm_modis_extraction.py --species puma --model_id puma_modis_m2 \
      --run_mode full [--ee_project geersprocessing] [--batch_size 25] ...
"""

import ee, time, os, sys, argparse
from datetime import datetime

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Predictor Extraction')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',      type=str, default='geersprocessing',
                   help='GEE project ID')
    p.add_argument('--species',         type=str, required=True,
                   help='Species name, e.g. puma')
    p.add_argument('--model_id',        type=str, required=True,
                   help='Unique model run ID, e.g. puma_modis_m2')

    # ── Run control ───────────────────────────────────────────────────────────
    p.add_argument('--run_mode',        type=str, default='single',
                   choices=['single', 'full', 'resume'],
                   help='Execution mode')
    p.add_argument('--batch_year',      type=int, default=2015,
                   help='Year for single-batch test run')
    p.add_argument('--batch_num',       type=int, default=1,
                   help='Batch number for single-batch test run')
    p.add_argument('--batch_size',      type=int, default=25,
                   help='Points per batch')
    p.add_argument('--max_concurrent',  type=int, default=5,
                   help='Max concurrent GEE tasks')

    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--occurrence_asset', type=str, default=None,
                   help='GEE asset path for occurrence points '
                        '(default: {gee_assets}/{species}_occurrences_batched)')
    p.add_argument('--gee_assets',      type=str,
                   default='projects/geersprocessing/assets',
                   help='Root GEE assets folder (no trailing slash)')

    # ── Temporal range (shared across scripts) ────────────────────────────────
    p.add_argument('--min_year',        type=int, default=2000)
    p.add_argument('--max_year',        type=int, default=2025)
    p.add_argument('--gap_years',       type=int, default=2,
                   help='Prior years used for gap-fill')

    # ── Spatial / scale ───────────────────────────────────────────────────────
    p.add_argument('--target_scale',    type=int, default=250,
                   help='Export scale in metres')
    p.add_argument('--point_buffer',    type=int, default=15000,
                   help='Buffer around each point for predictor stack (m)')

    return p.parse_args()

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# All output paths derived from species + model_id — change those two and
# everything else follows automatically.
# =============================================================================

SPECIES          = args.species
MODEL_ID         = args.model_id
RUN_MODE         = args.run_mode
BATCH_YEAR       = args.batch_year
BATCH_NUM        = args.batch_num
BATCH_SIZE       = args.batch_size
MAX_CONCURRENT   = args.max_concurrent
MIN_YEAR         = args.min_year
MAX_YEAR         = args.max_year
GAP_YEARS        = args.gap_years
TARGET_SCALE     = args.target_scale
POINT_BUFFER     = args.point_buffer

GEE_ASSETS       = args.gee_assets

# Input asset
OCCURRENCE_ASSET = (args.occurrence_asset or
                    f'{GEE_ASSETS}/{SPECIES}_occurrences_batched')

# Output folder — derived from SPECIES + MODEL_ID
EXPORT_FOLDER    = f'{GEE_ASSETS}/{SPECIES}_SDM_modis_exports_{MODEL_ID}'

# Log file — derived from SPECIES + MODEL_ID
LOG_FILE         = f'completed_batches_{SPECIES}_{MODEL_ID}.txt'

# =============================================================================
# 2. FIXED PARAMETERS (not user-facing)
# =============================================================================

SR_SCALE    = 10000
IDX_BANDS   = ['NDVI', 'EVI']

INTERVALS = [
    (1,16,'d001'), (17,32,'d017'), (33,48,'d033'), (49,64,'d049'),
    (65,80,'d065'), (81,96,'d081'), (97,112,'d097'), (113,128,'d113'),
    (129,144,'d129'), (145,160,'d145'), (161,176,'d161'), (177,192,'d177'),
    (193,208,'d193'), (209,224,'d209'), (225,240,'d225'), (241,256,'d241'),
    (257,272,'d257'), (273,288,'d273'), (289,304,'d289'), (305,320,'d305'),
    (321,336,'d321'), (337,352,'d337'), (353,366,'d353'),
]
N_INTERVALS = len(INTERVALS)
ALL_BANDS   = [f'{idx}_{iv[2]}' for iv in INTERVALS for idx in IDX_BANDS]

# Gaussian kernels — smooth neighborhood transitions.
# sigma = radius / 3: weight at kernel edge ~1% of centre weight.
KERNEL_250M = ee.Kernel.gaussian(125,    125/3,   'meters')
KERNEL_1KM  = ee.Kernel.gaussian(1000,   1000/3,  'meters')
KERNEL_5KM  = ee.Kernel.gaussian(5000,   5000/3,  'meters')
KERNEL_10KM = ee.Kernel.gaussian(10000,  10000/3, 'meters')

# Epoch-matched LCLUC and GHSL
LCLUC = {y: ee.Image(f'projects/glad/GLCLU2020/v2/LCLUC_{y}')
         for y in [2000,2005,2010,2015,2020]}
GHSL  = {y: ee.Image(f'JRC/GHSL/P2023A/GHS_POP/{y}')
         for y in [2000,2005,2010,2015,2020]}

def _nearest(d, year):
    year = max(MIN_YEAR, min(MAX_YEAR, year))
    if year < 2003: return d[2000]
    if year < 2008: return d[2005]
    if year < 2013: return d[2010]
    if year < 2018: return d[2015]
    return d[2020]

def get_nearest_lcluc(year): return _nearest(LCLUC, year)
def get_nearest_ghsl(year):  return _nearest(GHSL,  year)

# =============================================================================
# 3. TASK MANAGEMENT
# =============================================================================

def wait_for_tasks(tasks, check_interval=60):
    pending = list(tasks)
    while pending:
        still = []
        for task in pending:
            s = task.status()
            if s['state'] == 'COMPLETED':
                eecu = s.get('batch_eecu_usage_seconds', 0)
                print(f"OK {s['description']} ({eecu:.1f} EECU-s)")
            elif s['state'] == 'FAILED':
                print(f"FAIL {s['description']}: {s.get('error_message','')}")
            else:
                still.append(task)
        pending = still
        if pending:
            print(f"{datetime.now().strftime('%H:%M')} -- {len(pending)} running...")
            time.sleep(check_interval)
    print('All tasks complete.')


def submit_batch_queue(task_list, max_concurrent=5, check_interval=60,
                       log_file=None):
    completed, queue, running = [], list(task_list), []
    while queue or running:
        while len(running) < max_concurrent and queue:
            item = queue.pop(0)
            item['task'].start()
            running.append(item)
            print(f"Started: {item['description']} ({len(running)}/{max_concurrent})")
        still = []
        for item in running:
            s = item['task'].status()
            if s['state'] == 'COMPLETED':
                eecu = s.get('batch_eecu_usage_seconds', 0)
                print(f"OK {item['description']} -- {eecu:.1f} EECU-s")
                completed.append(item)
                if log_file:
                    log_completed(log_file, item['description'])
            elif s['state'] == 'FAILED':
                err = s.get('error_message', '')
                print(f"FAIL {item['description']} -- {err}")
                if not item.get('retried'):
                    item['retried'] = True
                    try:
                        item['task'].start()
                        still.append(item)
                    except Exception as e:
                        e_str = str(e).lower()
                        if 'request_id' in e_str or 'already started' in e_str:
                            still.append(item)
                        else:
                            print(f"Retry error: {e}")
                else:
                    print(f"Skipping after retry failure.")
            else:
                still.append(item)
        running = still
        if running or queue:
            time.sleep(check_interval)
    return completed


def load_completed(log_file):
    try:
        with open(log_file, 'r') as f:
            return set(line.strip() for line in f if line.strip())
    except FileNotFoundError:
        return set()


def log_completed(log_file, description):
    with open(log_file, 'a') as f:
        f.write(description + '\n')

# =============================================================================
# 4. MODIS COLLECTION BUILDER
# =============================================================================

def get_modis_collection(year, aoi):
    """Merge Terra + Aqua 16-day 250m VI, good pixels only."""
    start = ee.Date.fromYMD(year, 1, 1)
    end   = ee.Date.fromYMD(year, 12, 31)

    def mask_and_scale(img):
        qa   = img.select('SummaryQA')
        good = qa.eq(0)
        ndvi = img.select('NDVI').updateMask(good).rename('NDVI')
        evi  = img.select('EVI').updateMask(good).rename('EVI')
        return (ndvi.addBands(evi)
                    .copyProperties(img, ['system:time_start', 'system:index']))

    terra = (ee.ImageCollection('MODIS/061/MOD13Q1')
               .filterDate(start, end).filterBounds(aoi)
               .map(mask_and_scale))
    aqua  = (ee.ImageCollection('MODIS/061/MYD13A1')
               .filterDate(start, end).filterBounds(aoi)
               .map(mask_and_scale))
    return terra.merge(aqua)

# =============================================================================
# 5. INTERVAL COMPOSITE
# =============================================================================

def build_interval_composite(year, iv_idx, year_col):
    doy_start, doy_end, label = INTERVALS[iv_idx]
    interval_col = year_col.filter(
        ee.Filter.calendarRange(doy_start, doy_end, 'day_of_year'))
    n_images = interval_col.size()
    bands    = []
    for idx in IDX_BANDS:
        final_name = f'{idx}_{label}'
        composite  = ee.Image(ee.Algorithms.If(
            n_images.gt(0),
            interval_col.select(idx).median().rename(final_name).toInt16(),
            ee.Image.constant(0).rename(final_name).toInt16()
                        .updateMask(ee.Image(0))
        ))
        bands.append(composite)
    img = bands[0]
    for b in bands[1:]:
        img = img.addBands(b)
    return (img.set('iv_idx', iv_idx, 'label', label,
                    'doy_start', doy_start, 'doy_end', doy_end, 'year', year,
                    'system:time_start',
                    ee.Date.fromYMD(year, 1, 1)
                           .advance(doy_start - 1, 'day').millis()))

# =============================================================================
# 6. GAP-FILL
# =============================================================================

def build_year_composites(year, aoi):
    col = get_modis_collection(year, aoi)
    return [build_interval_composite(year, i, col) for i in range(N_INTERVALS)]


def gap_fill_composites(target_composites, target_year, aoi):
    eff_gap = min(GAP_YEARS, target_year - MIN_YEAR)
    prior   = []
    for y in range(1, eff_gap + 1):
        c = get_modis_collection(target_year - y, aoi)
        prior.append([build_interval_composite(target_year - y, i, c)
                      for i in range(N_INTERVALS)])
    filled = []
    for i in range(N_INTERVALS):
        tgt      = ee.Image(target_composites[i])
        tobs     = tgt.select(f'{IDX_BANDS[0]}_{INTERVALS[i][2]}').mask()
        prev_obs = ee.Image(target_composites[max(0, i-1)]).select(
            f'{IDX_BANDS[0]}_{INTERVALS[max(0,i-1)][2]}').mask()
        next_obs = ee.Image(target_composites[min(N_INTERVALS-1, i+1)]).select(
            f'{IDX_BANDS[0]}_{INTERVALS[min(N_INTERVALS-1,i+1)][2]}').mask()
        long_gap = tobs.Not().And(prev_obs.Not().Or(next_obs.Not()))
        fb, fo = tgt, tobs
        for y in range(len(prior)):
            pc  = ee.Image(prior[y][i])
            po  = pc.select(f'{IDX_BANDS[0]}_{INTERVALS[i][2]}').mask()
            use = long_gap.And(fo.Not()).And(po)
            fb  = fb.where(use, pc)
            fo  = fo.Or(use)
        filled.append(ee.Image(fb).copyProperties(
            tgt, ['iv_idx','label','doy_start','doy_end','year','system:time_start']))

    interp = []
    for i in range(N_INTERVALS):
        comp      = ee.Image(filled[i])
        all_bands = [f'{idx}_{INTERVALS[i][2]}' for idx in IDX_BANDS]
        gap       = comp.select(all_bands[0]).mask().Not()
        pv, pd = None, 0
        for p in range(i-1, max(-1, i-4), -1):
            if pv is None: pv, pd = ee.Image(filled[p]), i - p
        nv, nd = None, 0
        for n in range(i+1, min(N_INTERVALS, i+4)):
            if nv is None: nv, nd = ee.Image(filled[n]), n - i
        prior_bands = ([f'{idx}_{INTERVALS[i-pd][2]}' for idx in IDX_BANDS]
                       if pv is not None else None)
        next_bands  = ([f'{idx}_{INTERVALS[i+nd][2]}' for idx in IDX_BANDS]
                       if nv is not None else None)
        if pv is not None and nv is not None:
            t  = pd + nd
            ib = (pv.select(prior_bands).rename(all_bands).multiply(nd/t)
                    .add(nv.select(next_bands).rename(all_bands).multiply(pd/t))
                    .toInt16())
        elif pv is not None:
            ib = pv.select(prior_bands).rename(all_bands).toInt16()
        elif nv is not None:
            ib = nv.select(next_bands).rename(all_bands).toInt16()
        else:
            ib = comp.select(all_bands)
        fb = comp.select(all_bands).where(gap, ib)
        interp.append(ee.Image(fb).copyProperties(
            comp, ['iv_idx','label','doy_start','doy_end','year','system:time_start']))
    return interp

# =============================================================================
# 7. METRICS
# =============================================================================

def extract_metrics(filled):
    all_bands = []
    for idx in IDX_BANDS:
        col = ee.ImageCollection([
            ee.Image(filled[i]).select([f'{idx}_{INTERVALS[i][2]}']).rename(idx)
            for i in range(N_INTERVALS)
        ])
        p25 = col.reduce(ee.Reducer.percentile([25])).rename(f'{idx}_p25')
        p75 = col.reduce(ee.Reducer.percentile([75])).rename(f'{idx}_p75')
        p50 = col.median()                            .rename(f'{idx}_p50')
        amp = p75.subtract(p25)                       .rename(f'{idx}_amp')
        all_bands.extend([p25, p50, p75, amp])
    result = all_bands[0]
    for b in all_bands[1:]:
        result = result.addBands(b)
    return result


def build_modis_stack(year, aoi):
    target = build_year_composites(year, aoi)
    filled = gap_fill_composites(target, year, aoi)
    return extract_metrics(filled)

# =============================================================================
# 8. NEIGHBORHOOD SCALES
# =============================================================================

def add_neighborhood_scales(image, aoi=None):
    img   = image.clip(aoi.buffer(20000)) if aoi else image
    bands = img.bandNames()
    sfx   = lambda s: bands.map(lambda b: ee.String(b).cat(s))
    return (img.reduceNeighborhood(ee.Reducer.mean(), KERNEL_250M) .rename(sfx('_250m'))
               .addBands(img.reduceNeighborhood(ee.Reducer.mean(), KERNEL_1KM)  .rename(sfx('_1km')))
               .addBands(img.reduceNeighborhood(ee.Reducer.mean(), KERNEL_5KM)  .rename(sfx('_5km')))
               .addBands(img.reduceNeighborhood(ee.Reducer.mean(), KERNEL_10KM) .rename(sfx('_10km'))))

# =============================================================================
# 9. LCLUC STACK
# =============================================================================

def get_lcluc_stack(lcluc, aoi):
    is_terra_short = lcluc.gte(0).And(lcluc.lte(24))
    is_terra_trees = lcluc.gte(25).And(lcluc.lte(48))
    is_wet_short   = lcluc.gte(100).And(lcluc.lte(124))
    is_wet_trees   = lcluc.gte(125).And(lcluc.lte(148))
    is_water       = lcluc.gte(200).And(lcluc.lte(207))
    is_cropland    = lcluc.eq(244)
    is_builtup     = lcluc.eq(250)
    terra_short_pct = (lcluc.multiply(4).add(3).updateMask(is_terra_short)
                            .toFloat().rename('lcluc_terra_short_pct'))
    wet_short_pct   = (lcluc.subtract(100).multiply(4).add(3).updateMask(is_wet_short)
                            .toFloat().rename('lcluc_wet_short_pct'))
    terra_tree_ht   = (lcluc.subtract(22).updateMask(is_terra_trees)
                            .toFloat().rename('lcluc_terra_ht'))
    wet_tree_ht     = (lcluc.subtract(122).updateMask(is_wet_trees)
                            .toFloat().rename('lcluc_wet_ht'))
    stack = ee.Image.cat([
        is_terra_short.rename('lcluc_short_veg').toFloat(),
        is_terra_trees.rename('lcluc_trees').toFloat(),
        is_wet_short.rename('lcluc_wet_short').toFloat(),
        is_wet_trees.rename('lcluc_wet_trees').toFloat(),
        is_water.rename('lcluc_water').toFloat(),
        is_cropland.rename('lcluc_cropland').toFloat(),
        is_builtup.rename('lcluc_builtup').toFloat(),
        terra_short_pct, wet_short_pct, terra_tree_ht, wet_tree_ht,
    ])
    return add_neighborhood_scales(stack, aoi)

# =============================================================================
# 10. FULL PREDICTOR STACK
# =============================================================================

def build_predictor_stack(year, aoi):
    yr   = max(MIN_YEAR, min(MAX_YEAR, year))
    srtm = ee.Image('USGS/SRTMGL1_003')
    topo = (srtm.select('elevation').rename('elevation')
                .addBands(ee.Terrain.slope(srtm).rename('slope'))
                .addBands(ee.Image('CSP/ERGo/1_0/Global/SRTM_CHILI').rename('chili'))
                .addBands(ee.Image('CSP/ERGo/1_0/Global/SRTM_mTPI').rename('mtpi')))
    ghsl = (get_nearest_ghsl(yr).select([0]).rename('pop_density')
                                .max(ee.Image(0)).toFloat())
    return (add_neighborhood_scales(ee.Image('WORLDCLIM/V1/BIO'), aoi)
            .addBands(add_neighborhood_scales(topo, aoi))
            .addBands(add_neighborhood_scales(build_modis_stack(yr, aoi), aoi))
            .addBands(get_lcluc_stack(get_nearest_lcluc(yr), aoi))
            .addBands(add_neighborhood_scales(ghsl, aoi))
            .toFloat())

# =============================================================================
# 11. BATCH BUILDER
# =============================================================================

def build_batches(occurrences, batch_size=25):
    print('Fetching unique years from asset...')
    years = occurrences.aggregate_array('year').distinct().sort().getInfo()
    years = [y for y in years if y >= MIN_YEAR]
    print(f'Unique years (MODIS era): {years}')
    batches = []
    for yi, year in enumerate(years):
        print(f'  Fetching IDs for {year} ({yi+1}/{len(years)})...')
        ids = sorted(occurrences.filter(ee.Filter.eq('year', year))
                                .aggregate_array('system:index').getInfo())
        for i in range(-(-len(ids) // batch_size)):
            batches.append({'year': year, 'batch': i+1,
                            'ids': ids[i*batch_size:(i+1)*batch_size]})
    sizes = [len(b['ids']) for b in batches]
    print(f'Total batches: {len(batches)}, pts: min={min(sizes)} max={max(sizes)}')
    return batches

# =============================================================================
# 12. EXTRACTION TASK BUILDER
# =============================================================================

def build_extraction_task(occurrences, batch_info):
    year, batch_num, ids = batch_info['year'], batch_info['batch'], batch_info['ids']
    pts = (occurrences
           .filter(ee.Filter.eq('year', year))
           .filter(ee.Filter.inList('system:index', ids)))

    def ep(point):
        aoi = point.geometry().buffer(POINT_BUFFER).bounds()
        return point.set(
            build_predictor_stack(year, aoi)
            .unmask(-9999)
            .reduceRegion(reducer=ee.Reducer.first(),
                          geometry=point.geometry(),
                          scale=TARGET_SCALE,
                          tileScale=8))

    desc = f'SDM_MODIS_{SPECIES}_{year}_b{batch_num:02d}'
    task = ee.batch.Export.table.toAsset(
        collection=pts.map(ep),
        description=desc,
        assetId=f'{EXPORT_FOLDER}/{desc}')
    return {'task': task, 'description': desc}

# =============================================================================
# 13. OUTPUT FOLDER SETUP
# =============================================================================

def ensure_folder():
    try:
        ee.data.createAsset({'type': 'Folder'}, EXPORT_FOLDER)
        print(f'  Created: {EXPORT_FOLDER}')
    except Exception as e:
        e_str = str(e).lower()
        if 'already exists' in e_str or 'cannot overwrite' in e_str:
            print(f'  Exists:  {EXPORT_FOLDER}')
        else:
            print(f'  Error:   {EXPORT_FOLDER}: {e}')

# =============================================================================
# 14. MAIN
# =============================================================================

if __name__ == '__main__':
    occ = ee.FeatureCollection(OCCURRENCE_ASSET)
    print(f'Species:          {SPECIES}')
    print(f'Model ID:         {MODEL_ID}')
    print(f'Export folder:    {EXPORT_FOLDER}')
    print(f'Log file:         {LOG_FILE}')
    print(f'Total points:     {occ.size().getInfo()}')
    ensure_folder()

    if RUN_MODE == 'single':
        print(f'\nSingle-batch: year={BATCH_YEAR}, batch={BATCH_NUM}')
        pts = (occ.filter(ee.Filter.eq('year', BATCH_YEAR))
                  .filter(ee.Filter.eq('batch', BATCH_NUM)))
        print(f'Points: {pts.size().getInfo()}')

        def ep(point):
            aoi = point.geometry().buffer(POINT_BUFFER).bounds()
            return point.set(
                build_predictor_stack(BATCH_YEAR, aoi)
                .unmask(-9999)
                .reduceRegion(reducer=ee.Reducer.first(),
                              geometry=point.geometry(),
                              scale=TARGET_SCALE,
                              tileScale=8))

        desc = f'SDM_MODIS_{SPECIES}_{BATCH_YEAR}_b{BATCH_NUM:02d}'
        task = ee.batch.Export.table.toAsset(
            collection=pts.map(ep),
            description=desc,
            assetId=f'{EXPORT_FOLDER}/{desc}')
        task.start()
        print(f'Task started: {desc}')
        wait_for_tasks([task])

    else:
        is_resume  = (RUN_MODE == 'resume')
        mode_label = 'Resume' if is_resume else 'Full'
        print(f'Mode: {mode_label} | batch_size={BATCH_SIZE}, '
              f'max_concurrent={MAX_CONCURRENT}')

        batches = build_batches(occ, BATCH_SIZE)

        if is_resume:
            already_done = load_completed(LOG_FILE)
            def batch_desc(b):
                return f"SDM_MODIS_{SPECIES}_{b['year']}_b{b['batch']:02d}"
            batches = [b for b in batches if batch_desc(b) not in already_done]
            print(f'Resume: {len(batches)} batches remaining.')
        else:
            try:
                os.remove(LOG_FILE)
                print(f'Cleared log: {LOG_FILE}')
            except FileNotFoundError:
                pass

        if not batches:
            print('All batches complete.')
            sys.exit(0)

        print('Building task list...')
        task_list = [build_extraction_task(occ, b) for b in batches]
        print(f'{len(task_list)} tasks. Starting submission...')
        completed = submit_batch_queue(task_list, MAX_CONCURRENT,
                                       log_file=LOG_FILE)
        total_eecu = sum(t['task'].status().get('batch_eecu_usage_seconds', 0)
                         for t in completed)
        print(f'\nDone: {len(completed)}/{len(task_list)} batches, '
              f'{total_eecu:,.0f} EECU-s')
        if len(completed) < len(task_list):
            print(f'Set --run_mode resume and rerun.')
