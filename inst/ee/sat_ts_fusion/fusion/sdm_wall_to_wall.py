"""
SDM Wall-to-Wall Prediction Script
==================================
Applies a pre-trained RF classifier tile by tile to produce a wall-to-wall
habitat suitability map.

The predictor stack is built by predictors.assemble_stack(), the same
function extraction uses for training data. Imagery-source metrics are
computed in memory, as in extraction, unless --use_exported_metrics reads
MODIS metric tiles from sdm_modis_export.py (cheaper only when 2+ models
map the same species and year; identical values, 2026-10-09). The predictor
source and scale are read from the model's features asset, written by
sdm_model_fitting.py from the training data, so they cannot be mismatched.

Usage:
  python sdm_wall_to_wall.py --species puma --model_id puma_modis_m3 \
      --target_year 2025 [--run_mode full] [--ee_project geersprocessing] ...
"""

import ee, time, argparse, json
from datetime import datetime
import predictors as P
import grid as G
import cost as C

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM Wall-to-Wall Prediction')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--species',      type=str, required=True)
    p.add_argument('--model_id',     type=str, required=True,
                   help='Model run ID - must match sdm_model_fitting.py')
    p.add_argument('--target_year',  type=int, required=True,
                   help='Year to map (with --use_exported_metrics, must match '
                        'sdm_modis_export.py)')

    # ── Run control ───────────────────────────────────────────────────────────
    p.add_argument('--run_mode',      type=str, default='test',
                   choices=['test', 'full'])
    p.add_argument('--max_concurrent', type=int, default=3,
                   help='Max concurrent GEE batch tasks. Leave headroom '
                        'for other work (default 3).')
    p.add_argument('--estimate_only', action='store_true',
                   help='Print the cost estimate and exit')
    p.add_argument('--tile_degrees', type=float, default=None,
                   help='0 or omitted: chosen from the model scale; with '
                        '--use_exported_metrics, must match export')

    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--gee_assets',   type=str,
                   default='projects/geersprocessing/assets')
    p.add_argument('--range_asset',  type=str, default=None,
                   help='Default: {gee_assets}/{species}_iucn_range')
    p.add_argument('--use_exported_metrics', action='store_true',
                   help='Read MODIS metrics from sdm_modis_export.py tiles '
                        'instead of computing them in memory')
    p.add_argument('--metrics_folder', type=str, default=None,
                   help='With --use_exported_metrics: override the metric tiles '
                        'folder (default derived from species, source/scale '
                        'and target year)')

    return p.parse_args()

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

SPECIES        = args.species
MODEL_ID       = args.model_id
TARGET_YEAR    = args.target_year
RUN_MODE       = args.run_mode
MAX_CONCURRENT = args.max_concurrent
CRS            = P.CRS
GEE_ASSETS     = args.gee_assets

RANGE_ASSET      = args.range_asset or f'{GEE_ASSETS}/{SPECIES}_iucn_range'
CLASSIFIER_ASSET = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier'
STRINGS_ASSET    = f'{GEE_ASSETS}/{MODEL_ID}_RF_classifier_strings'
FEATURES_ASSET   = f'{GEE_ASSETS}/{MODEL_ID}_selected_features'
OUTPUT_FOLDER    = f'{GEE_ASSETS}/{MODEL_ID}_prediction_{TARGET_YEAR}'

# =============================================================================
# 2. MODEL
# =============================================================================

def load_classifier():
    print(f'\nLoading classifier from {STRINGS_ASSET}...')
    try:
        from geemap import ml
        strings_fc  = ee.FeatureCollection(STRINGS_ASSET)
        output_mode = strings_fc.first().get('output_mode').getInfo()
        ee_strings  = (strings_fc.sort('tree_index')
                                  .aggregate_array('tree_string').getInfo())
        classifier  = ml.strings_to_classifier(ee_strings).setOutputMode(output_mode)
        print(f'Loaded {len(ee_strings)} trees, mode={output_mode}')
        return classifier
    except Exception as e:
        print(f'Warning: strings load failed ({e}), falling back to Classifier.load()')
        return ee.Classifier.load(CLASSIFIER_ASSET)


def load_features():
    """selected features, fill image, source, scale from the features asset."""
    print(f'Loading features from {FEATURES_ASSET}...')
    info     = ee.FeatureCollection(FEATURES_ASSET).first().toDictionary().getInfo()
    features = info['selected_features'].split(',')
    print(f'Features ({len(features)}): {features}')

    if 'predictor_source' in info and 'scale' in info:
        source, scale = info['predictor_source'], int(info['scale'])
    else:
        source, scale = 'modis', P.SOURCES['modis']['native_scale']
        print(f'WARNING: no predictor_source/scale on {FEATURES_ASSET} (model '
              f'fit before 2026-10-07) - assuming {source} {scale} m.')

    # Masked predictors get the same fill values fitting used for -9999
    # sentinels (training median, or 0 for LCLUC sub-class bands).
    if 'fill_values' in info:
        fills  = json.loads(info['fill_values'])
        values = [float(fills[f]) for f in features]
        print('Fill values loaded from features asset')
    else:
        values = [P.SENTINEL] * len(features)
        print(f'WARNING: no fill_values on {FEATURES_ASSET} (model fit before '
              f'2026-10-06) - masked predictors passed as {P.SENTINEL}, which '
              f'does not match how training data was cleaned. Refit the model.')
    fill_image = ee.Image.constant(values).rename(features).toFloat()

    gap_years, stamped = P.recorded_gap_years([info.get('gap_years')])
    if not stamped:
        print(f'WARNING: no gap_years on {FEATURES_ASSET} (model fit before '
              f'2026-10-09) - assuming extraction used the default {gap_years}.')
    return features, fill_image, source, scale, gap_years

# =============================================================================
# 3. TASK MANAGEMENT
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
            print(f'  {datetime.now().strftime("%H:%M")} -- '
                  f'{len(running)} running, {len(queue)} queued...')
            time.sleep(check_interval)

    print(f'  Completed: {len(completed)} | Failed: {len(failed)} | '
          f'Total EECU-s: {total_eecu:,.1f}')
    return completed, failed, total_eecu

# =============================================================================
# 4. METRIC TILES
# =============================================================================

def list_metric_tiles(metrics_folder):
    """Names of tiles with metric assets (one paged listing, not per-tile checks)."""
    names, token = set(), None
    while True:
        params = {'parent': metrics_folder, 'pageSize': 1000}
        if token: params['pageToken'] = token
        try:
            resp = ee.data.listAssets(params)
        except Exception:
            return names
        names |= {a['name'].split('/')[-1][:-len('_metrics')]
                  for a in resp.get('assets', [])}
        token = resp.get('nextPageToken')
        if not token:
            return names


def neighbour_tile_names(tile, existing, tile_degrees):
    """Tile plus its grid neighbours that have metric assets. Kernels reach
    10 km, so neighbours avoid truncated neighbourhoods at tile edges."""
    names = [tile['name']]
    for x, y in G.neighbour_cells(tile['lon'], tile['lat'], tile_degrees):
        n = G.tile_name(x, y, tile_degrees)
        if n in existing and n not in names:
            names.append(n)
    return names


def metrics_mosaic(metrics_folder, tile_names):
    return (ee.ImageCollection([ee.Image(f'{metrics_folder}/{n}_metrics')
                                for n in tile_names])
              .mosaic().select(P.METRIC_BANDS))

# =============================================================================
# 5. ASSET UTILITIES
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

def check_features_buildable(features, source, scale, prior, region):
    """Stop before any task if the model needs predictors the stack cannot
    build. Seen 2026-10-20: a model listed 'yyyy', a column of the occurrence
    data that the extraction of the time copied into the training tables;
    every tile would have failed at select(). Exported metric tiles have the
    same band names as the in-memory metrics."""
    metrics = P.source_predictors(source, TARGET_YEAR, prior, region)
    names   = P.assemble_stack(source, metrics, TARGET_YEAR, scale, CRS).bandNames()
    missing = ee.List(features).removeAll(names).getInfo()
    if missing:
        raise SystemExit(
            f'Model {MODEL_ID} uses {len(missing)} predictor(s) that the '
            f'{source} stack does not have: {missing}. These usually come '
            f'from extra columns of the occurrence data (training data '
            f'extracted before 2026-10-12, when extraction began keeping only '
            f'year and the response column). Re-extract and refit the model.')

# =============================================================================
# 6. MAIN
# =============================================================================

def main():
    classifier = load_classifier()
    features, fill_image, source, scale, gap_years = load_features()
    tile_degrees   = args.tile_degrees or G.auto_tile_degrees(scale)
    exported       = args.use_exported_metrics
    if exported and not P.SOURCES[source]['exportable']:
        raise SystemExit(f'--use_exported_metrics: {source} has no export stage')
    if args.metrics_folder and not exported:
        raise SystemExit('--metrics_folder needs --use_exported_metrics')
    metrics_folder = (args.metrics_folder or
                      f'{GEE_ASSETS}/{SPECIES}_SDM_{P.source_tag(source, scale)}_{TARGET_YEAR}/metrics')
    usable = P.coverage_years(source)
    if TARGET_YEAR not in usable:
        raise SystemExit(f'--target_year {TARGET_YEAR} is outside {source} '
                         f'coverage ({usable[0]}-{usable[-1]})')

    print('=' * 60)
    print(f'SDM WALL-TO-WALL PREDICTION')
    print(f'Species: {SPECIES} | Model ID: {MODEL_ID} | Year: {TARGET_YEAR}')
    print(f'Source: {source} | Scale: {scale} m | Tiles: {tile_degrees} deg | Mode: {RUN_MODE}')
    print(f'Metrics:  {metrics_folder if exported else "computed in memory from " + source}')
    if exported:
        # Gap-fill was fixed when the metric tiles were exported.
        print(f'Gap-fill: as exported (model trained with --gap_years {gap_years}; '
              f'export should match)')
    elif P.SOURCES[source]['coverage'] == 'annual':
        print('Gap-fill: none (annual source)')
    else:
        print(f'Gap-fill: prior years {P.prior_years(TARGET_YEAR, gap_years, usable[0])} '
              f'(--gap_years {gap_years} from training)')
    print(f'Output:   {OUTPUT_FOLDER}')
    print('=' * 60)

    if RUN_MODE == 'test':
        tiles = G.test_tiles(tile_degrees)
        print(f'\nTest tiles: {[t["name"] for t in tiles]}')
    else:
        print(f'Loading range from {RANGE_ASSET}...')
        tiles = G.build_tiles(ee.FeatureCollection(RANGE_ASSET), tile_degrees)

    # Band names do not depend on the tile, so one check covers them all.
    if tiles:
        check_features_buildable(features, source, scale,
                                 P.prior_years(TARGET_YEAR, gap_years, usable[0]),
                                 tiles[0]['geom'])

    print(C.report([C.prediction(source, len(tiles), tile_degrees, scale,
                                 max_concurrent=MAX_CONCURRENT, exported=exported)],
                   'Cost estimate (all tiles; tiles already done are skipped)'))
    if args.estimate_only:
        return

    ensure_folder(OUTPUT_FOLDER)
    if exported:
        existing = list_metric_tiles(metrics_folder)
        print(f'Metric tiles available: {len(existing)}')

    print(f'\nBuilding prediction tasks (max_concurrent={MAX_CONCURRENT})...')
    task_list, skipped = [], 0

    for tile in tiles:
        asset_id = f'{OUTPUT_FOLDER}/{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}'
        if asset_exists(asset_id):
            print(f'  Skipping (exists): {tile["name"]}')
            skipped += 1
            continue
        if exported:
            if tile['name'] not in existing:
                print(f'  Skipping (no metrics): {tile["name"]}')
                skipped += 1
                continue
            metrics = metrics_mosaic(metrics_folder,
                                     neighbour_tile_names(tile, existing, tile_degrees))
        else:
            # Same call (and prior years) extraction uses; region only limits
            # which source images are loaded, buffered beyond the largest
            # kernel. Embeddings ignore the prior years.
            metrics = P.source_predictors(
                source, TARGET_YEAR, P.prior_years(TARGET_YEAR, gap_years, usable[0]),
                tile['geom'].buffer(1.5 * max(P.KERNEL_RADII_M)))
        stack   = P.assemble_stack(source, metrics, TARGET_YEAR, scale, CRS)

        missing = ee.List(features).removeAll(stack.bandNames()).getInfo()
        if missing:
            print(f'  WARNING {tile["name"]}: {len(missing)} features missing: {missing}')

        # Masked-predictor fraction check (these pixels get fill values)
        sel_stack   = stack.select(features)
        masked_frac = (sel_stack.mask().eq(0).reduce(ee.Reducer.mean())
                                .reduceRegion(reducer=ee.Reducer.mean(),
                                              geometry=tile['geom'],
                                              scale=scale * 4,
                                              maxPixels=1e6).getInfo())
        mean_masked = (sum(masked_frac.values()) / len(masked_frac)
                       if masked_frac else 0)
        if mean_masked > 0.01:
            print(f'  WARNING {tile["name"]}: {mean_masked:.1%} masked (filled) '
                  f'in selected features')
        else:
            print(f'  Queued ({mean_masked:.2%} masked) - {tile["name"]}')

        suitability = (sel_stack.unmask(fill_image).classify(classifier)
                                .rename('suitability').toFloat())

        # GLAD land/water mask
        lcluc = ee.Image('projects/glad/GLCLU2020/v2/LCLUC')
        suitability = (suitability
                       .updateMask(ee.Image('projects/glad/OceanMask').lte(1))
                       .updateMask(lcluc.lt(208).Or(lcluc.gt(211))))

        desc = G.safe_name(f'{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}')
        task = ee.batch.Export.image.toAsset(
            image=suitability, description=desc, assetId=asset_id,
            region=tile['geom'], scale=scale, crs=CRS, maxPixels=1e10)
        task_list.append({'task': task, 'desc': desc})

    print(f'\nQueued: {len(task_list)} | Skipped: {skipped}')
    if task_list:
        completed, failed, total_eecu = submit_and_wait(
            task_list, max_concurrent=MAX_CONCURRENT)
        print(f'\nComplete. Output: {OUTPUT_FOLDER}')
        if failed:
            print(f'WARNING: {len(failed)} tiles failed - rerun to retry.')

if __name__ == '__main__':
    main()
