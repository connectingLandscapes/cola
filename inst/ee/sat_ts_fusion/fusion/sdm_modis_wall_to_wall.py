"""
SDM MODIS Wall-to-Wall Prediction Script
==========================================
Applies a pre-trained RF classifier to pre-exported MODIS metrics tiles
to produce a wall-to-wall habitat suitability map.

All inputs derived from MODEL_ID (classifier, features) and
SPECIES + TARGET_YEAR (metrics tiles). Output folder also
derived from MODEL_ID + TARGET_YEAR.

IMPORTANT: Gaussian kernel parameters must be identical to those used in
  sdm_modis_extraction.py — predictor values must match between training
  and prediction or model performance will degrade.

Usage:
  python sdm_modis_wall_to_wall.py --species puma --model_id puma_modis_m2 \
      --target_year 2025 [--run_mode full] [--ee_project geersprocessing] ...
"""

import ee, time, math, argparse

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Wall-to-Wall Prediction')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--species',      type=str, required=True)
    p.add_argument('--model_id',     type=str, required=True,
                   help='Model run ID — must match sdm_model_fitting.py')
    p.add_argument('--target_year',  type=int, required=True,
                   help='Year to map — must match sdm_modis_export.py')

    # ── Run control ───────────────────────────────────────────────────────────
    p.add_argument('--run_mode',      type=str, default='test',
                   choices=['test', 'full'])
    p.add_argument('--max_concurrent', type=int, default=3,
                   help='Max concurrent GEE batch tasks. Leave headroom '
                        'for other work (default 3).')

    # ── Spatial / scale (shared across scripts) ───────────────────────────────
    p.add_argument('--crs',          type=str, default='EPSG:4326')
    p.add_argument('--scale',        type=int, default=250)
    p.add_argument('--tile_degrees', type=int, default=2,
                   help='Must match TILE_DEGREES used in sdm_modis_export.py')

    # ── Temporal (shared across scripts) ─────────────────────────────────────
    p.add_argument('--min_year',     type=int, default=2000)
    p.add_argument('--max_year',     type=int, default=2025)

    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--gee_assets',   type=str,
                   default='projects/geersprocessing/assets')
    p.add_argument('--range_asset',  type=str, default=None,
                   help='Default: {gee_assets}/{species}_iucn_range')

    return p.parse_args()
## Output to be saved here:  desc = safe_name(f'{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}')

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

SPECIES      = args.species
MODEL_ID     = args.model_id
TARGET_YEAR  = args.target_year
RUN_MODE        = args.run_mode
MAX_CONCURRENT  = args.max_concurrent
CRS             = args.crs
SCALE        = args.scale
TILE_DEGREES = args.tile_degrees
MIN_YEAR     = args.min_year
MAX_YEAR     = args.max_year
GEE_ASSETS   = args.gee_assets

RANGE_ASSET  = (args.range_asset or f'{GEE_ASSETS}/{SPECIES}_iucn_range')

# Classifier / features — derived from MODEL_ID
CLASSIFIER_ASSET = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier'
STRINGS_ASSET    = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier_strings'
FEATURES_ASSET   = f'{GEE_ASSETS}/{MODEL_ID}_selected_features'

# Metrics — derived from SPECIES + TARGET_YEAR (model-independent)
METRICS_FOLDER   = f'{GEE_ASSETS}/{SPECIES}_SDM_modis_{TARGET_YEAR}/metrics'

# Output — derived from MODEL_ID + TARGET_YEAR
OUTPUT_FOLDER    = f'{GEE_ASSETS}/{MODEL_ID}_prediction_{TARGET_YEAR}'

# =============================================================================
# 2. FIXED PARAMETERS
# =============================================================================

# Gaussian kernels — must be identical to sdm_modis_extraction.py.
# sigma = radius / 3: weight at kernel edge ~1% of centre weight.
KERNEL_250M = ee.Kernel.gaussian(125,    125/3,   'meters')
KERNEL_1KM  = ee.Kernel.gaussian(1000,   1000/3,  'meters')
KERNEL_5KM  = ee.Kernel.gaussian(5000,   5000/3,  'meters')
KERNEL_10KM = ee.Kernel.gaussian(10000,  10000/3, 'meters')

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
# 3. NEIGHBORHOOD SCALES
# =============================================================================

def add_neighborhood_scales(image):
    bands = image.bandNames()
    sfx   = lambda s: bands.map(lambda b: ee.String(b).cat(s))
    return (image.reduceNeighborhood(ee.Reducer.mean(), KERNEL_250M) .rename(sfx('_250m'))
                 .addBands(image.reduceNeighborhood(ee.Reducer.mean(), KERNEL_1KM)  .rename(sfx('_1km')))
                 .addBands(image.reduceNeighborhood(ee.Reducer.mean(), KERNEL_5KM)  .rename(sfx('_5km')))
                 .addBands(image.reduceNeighborhood(ee.Reducer.mean(), KERNEL_10KM) .rename(sfx('_10km'))))

# =============================================================================
# 4. LCLUC STACK
# =============================================================================

def get_lcluc_stack(lcluc):
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
    return add_neighborhood_scales(stack)

# =============================================================================
# 5. PREDICTOR STACK
# =============================================================================

def build_predictor_stack(year, tile_name):
    yr = max(MIN_YEAR, min(MAX_YEAR, year))

    # MODIS metrics from pre-exported asset
    modis_metrics = ee.Image(f'{METRICS_FOLDER}/{tile_name}_metrics')
    modis_bands   = modis_metrics.bandNames()
    modis_sfx     = lambda s: modis_bands.map(lambda b: ee.String(b).cat(s))
    modis_stack   = (modis_metrics.rename(modis_sfx('_250m'))
                                  .addBands(modis_metrics.reduceNeighborhood(
                                      ee.Reducer.mean(), KERNEL_1KM).rename(modis_sfx('_1km')))
                                  .addBands(modis_metrics.reduceNeighborhood(
                                      ee.Reducer.mean(), KERNEL_5KM).rename(modis_sfx('_5km')))
                                  .addBands(modis_metrics.reduceNeighborhood(
                                      ee.Reducer.mean(), KERNEL_10KM).rename(modis_sfx('_10km'))))

    # Drop source_year scale variants — not a predictor
    all_bands  = modis_stack.bandNames()
    keep_bands = all_bands.map(lambda b: ee.Algorithms.If(
        ee.String(b).index('source_year').gte(0), None, b)).removeAll([None])
    modis_stack = modis_stack.select(keep_bands)

    # WorldClim
    wc_stack = add_neighborhood_scales(ee.Image('WORLDCLIM/V1/BIO'))

    # Topography
    srtm       = ee.Image('USGS/SRTMGL1_003')
    topo_stack = (add_neighborhood_scales(
                      srtm.select('elevation').rename('elevation')
                          .addBands(ee.Terrain.slope(srtm).rename('slope')))
                  .addBands(add_neighborhood_scales(
                      ee.Image('CSP/ERGo/1_0/Global/SRTM_CHILI').rename('chili')))
                  .addBands(add_neighborhood_scales(
                      ee.Image('CSP/ERGo/1_0/Global/SRTM_mTPI').rename('mtpi'))))

    # LCLUC + population
    lcluc_stack = get_lcluc_stack(get_nearest_lcluc(yr))
    ghsl_stack  = add_neighborhood_scales(
        get_nearest_ghsl(yr).select([0]).rename('pop_density')
        .max(ee.Image(0)).toFloat())

    return (wc_stack
            .addBands(topo_stack)
            .addBands(modis_stack)
            .addBands(lcluc_stack)
            .addBands(ghsl_stack)
            .toFloat()
            .unmask(-9999))

# =============================================================================
# 6. TASK MANAGEMENT
# =============================================================================

def submit_and_wait(task_list, max_concurrent=3, check_interval=60):
    """Submit tasks up to max_concurrent at a time, wait for all to complete."""
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
                print(f'  OK {item["desc"]} ({eecu:.1f} EECU-s)')
                completed.append(item)
            elif s['state'] == 'FAILED':
                print(f'  FAIL {item["desc"]}: {s.get("error_message","")}')
                failed.append(item)
            else:
                still.append(item)
        running = still

        if running or queue:
            from datetime import datetime
            print(f'  {datetime.now().strftime("%H:%M")} -- '
                  f'{len(running)} running, {len(queue)} queued...')
            time.sleep(check_interval)

    print(f'  Completed: {len(completed)} | Failed: {len(failed)} | '
          f'Total EECU-s: {total_eecu:,.1f}')
    return completed, failed, total_eecu

# =============================================================================
# 7. TILING
# =============================================================================

def safe_name(s):
    """Replace '+' with 'p' — GEE disallows '+' in asset IDs and descriptions.
    Must match tile naming used in sdm_modis_export.py."""
    return s.replace('+', 'p')

def get_test_tiles():
    d = TILE_DEGREES
    return [
        {'geom': ee.Geometry.Rectangle([-84, 9, -84+d, 9+d]),
         'name': f'test_costa_rica_{d}deg'},
        {'geom': ee.Geometry.Rectangle([-84, 8, -84+d, 8+d]),
         'name': f'test_costa_rica_south_{d}deg'},
    ]

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

# =============================================================================
# 8. ASSET UTILITIES
# =============================================================================

def asset_exists(asset_id):
    try:
        ee.data.getAsset(asset_id)
        return True
    except Exception:
        return False

def ensure_folder(folder):
    try:
        ee.data.createAsset({'type': 'Folder'}, folder)
        print(f'  Created: {folder}')
    except Exception as e:
        e_str = str(e).lower()
        if 'already exists' in e_str or 'cannot overwrite' in e_str:
            print(f'  Exists:  {folder}')
        else:
            print(f'  Error:   {folder}: {e}')

# =============================================================================
# 9. MAIN
# =============================================================================

def main():
    print('=' * 60)
    print(f'SDM MODIS WALL-TO-WALL PREDICTION')
    print(f'Species: {SPECIES} | Model ID: {MODEL_ID} | Year: {TARGET_YEAR}')
    print(f'Mode: {RUN_MODE} | Tile: {TILE_DEGREES}°')
    print(f'Metrics:  {METRICS_FOLDER}')
    print(f'Output:   {OUTPUT_FOLDER}')
    print('=' * 60)

    # Load classifier
    print(f'\nLoading classifier from {STRINGS_ASSET}...')
    try:
        from geemap import ml
        strings_fc  = ee.FeatureCollection(STRINGS_ASSET)
        output_mode = strings_fc.first().get('output_mode').getInfo()
        ee_strings  = (strings_fc.sort('tree_index')
                                  .aggregate_array('tree_string').getInfo())
        classifier  = ml.strings_to_classifier(ee_strings).setOutputMode(output_mode)
        print(f'Loaded {len(ee_strings)} trees, mode={output_mode}')
    except Exception as e:
        print(f'Warning: strings load failed ({e}), falling back to Classifier.load()')
        classifier = ee.Classifier.load(CLASSIFIER_ASSET)

    # Load selected features
    print(f'Loading features from {FEATURES_ASSET}...')
    selected_features = (ee.FeatureCollection(FEATURES_ASSET)
                           .first().get('selected_features')
                           .getInfo().split(','))
   
    ## Remove column
    print(f'Features ({len(selected_features)}): {selected_features}')

    # Load range + build tiles
    if RUN_MODE == 'test':
        tiles = get_test_tiles()
        print(f'\nTest tiles: {[t["name"] for t in tiles]}')
    else:
        print(f'Loading range from {RANGE_ASSET}...')
        tiles = build_tiles(ee.FeatureCollection(RANGE_ASSET))

    ensure_folder(OUTPUT_FOLDER)

    print(f'\nBuilding prediction tasks (max_concurrent={MAX_CONCURRENT})...')
    task_list, skipped = [], 0

    for tile in tiles:
        asset_id = (f'{OUTPUT_FOLDER}/'
                    f'{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}')
        
        if asset_exists(asset_id):
            print(f'  Skipping (exists): {tile["name"]}')
            skipped += 1
            continue
        
        metrics_id = f'{METRICS_FOLDER}/{tile["name"]}_metrics'
        if not asset_exists(metrics_id):
            print(f'  Skipping (no metrics): {tile["name"]} | metrics_id: {metrics_id}')
            skipped += 1
            continue
        
        stack = build_predictor_stack(TARGET_YEAR, tile['name'])
        
        # Check all selected features are present
        missing = ee.List(selected_features).removeAll(stack.bandNames())
        n_miss  = missing.size().getInfo()
        if n_miss > 0:
            print(f'  WARNING {tile["name"]}: {n_miss} features missing: '
                  f'{missing.getInfo()}')
        

        # Sentinel fraction check
        sel_stack     = stack.select(selected_features)
        # print('all good 2')

        sentinel_frac = (sel_stack.eq(-9999).reduce(ee.Reducer.mean())
                                  .reduceRegion(reducer=ee.Reducer.mean(),
                                                geometry=tile['geom'],
                                                scale=SCALE * 4,
                                                maxPixels=1e6).getInfo())
        mean_sentinel = (sum(sentinel_frac.values()) / len(sentinel_frac)
                         if sentinel_frac else 0)
        if mean_sentinel > 0.01:
            print(f'  WARNING {tile["name"]}: {mean_sentinel:.1%} sentinels '
                  f'in selected features')
        else:
            print(f'  Queued ({mean_sentinel:.2%} sentinels) — {tile["name"]}')

        # Classify
        suitability = sel_stack.classify(classifier).rename('suitability').toFloat()

        # GLAD land/water mask
        suitability = (suitability
                       .updateMask(ee.Image('projects/glad/OceanMask').lte(1))
                       .updateMask(ee.Image('projects/glad/GLCLU2020/v2/LCLUC')
                                    .lt(208).Or(
                                   ee.Image('projects/glad/GLCLU2020/v2/LCLUC')
                                    .gt(211))))

        desc = safe_name(f'{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}')
        task = ee.batch.Export.image.toAsset(
            image=suitability, description=desc, assetId=asset_id,
            region=tile['geom'], scale=SCALE, crs=CRS, maxPixels=1e10)
        task_list.append({'task': task, 'desc': desc})

    print(f'\nQueued: {len(task_list)} | Skipped: {skipped}')
    if task_list:
        completed, failed, total_eecu = submit_and_wait(
            task_list, max_concurrent=MAX_CONCURRENT)
        print(f'\nComplete. Output: {OUTPUT_FOLDER}')
        if failed:
            print(f'⚠️  {len(failed)} tiles failed — rerun to retry.')

if __name__ == '__main__':
    main()
