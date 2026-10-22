"""
SDM MODIS Export Pipeline
==========================
OPTIONAL. Exports gap-filled MODIS NDVI/EVI metric tiles for
sdm_wall_to_wall.py --use_exported_metrics. Wall-to-wall computes the
same metrics in memory by default; exporting first is cheaper only when 2+
models map the same species and year (2 deg test tile, 2026-10-09: export
254 + prediction 74-87 EECU-s per model, vs 221 per model in memory).
Compositing, gap-fill and metrics come from predictors.py, the same
functions extraction uses in memory, so exported metrics equal training
metrics (verified 2026-10-06: 0 differences over 3,000 px).

All outputs are derived from SPECIES + source/scale + TARGET_YEAR
(model-independent), so the same metrics serve multiple model runs.

Pipeline stages (run in order):
  Stage 1 — export_annual:     23 x 16-day composites per year per tile
  Stage 2 — gap_fill:          long-gap prior-year merges + interpolation
  Stage 3 — reduce_to_metrics: p25/p50/p75/amp + source_year

Each stage skips tiles whose outputs already exist and resumes tiles part-way
through gap_fill, so rerun with the same arguments after any failure.

Usage:
  python sdm_modis_export.py --species puma --target_year 2025 \
      --stage export_annual --run_mode full [--scale 250] ...
"""

import ee, time, argparse
import predictors as P
import grid as G
import cost as C

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Export Pipeline')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--species',      type=str, required=True)
    p.add_argument('--target_year',  type=int, required=True)
    p.add_argument('--scale',        type=int, default=None,
                   help='Analysis scale in metres (default: MODIS native)')

    # ── Pipeline control ──────────────────────────────────────────────────────
    p.add_argument('--stage',        type=str, required=True,
                   choices=['export_annual', 'gap_fill', 'reduce_to_metrics'])
    p.add_argument('--run_mode',     type=str, default='full',
                   choices=['test', 'full'])
    p.add_argument('--max_concurrent', type=int, default=3,
                   help='Max concurrent GEE batch tasks. Leave headroom '
                        'for other work (default 3).')

    # ── Temporal / spatial ────────────────────────────────────────────────────
    p.add_argument('--gap_years',    type=int, default=2)
    p.add_argument('--estimate_only', action='store_true',
                   help='Print the cost estimate for this stage and exit')
    p.add_argument('--tile_degrees', type=float, default=None,
                   help='0 or omitted: chosen from --scale '
                        '(grid.auto_tile_degrees); must match wall-to-wall')
    p.add_argument('--keep_intermediate', action='store_true',
                   help='Keep the annual, merged and gap-filled assets in '
                        '{base}/temp instead of deleting each once the next '
                        'step has used it (for inspection; uses GEE storage)')

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

SOURCE         = 'modis'
SPECIES        = args.species
TARGET_YEAR    = args.target_year
SCALE          = P.resolve_scale(SOURCE, args.scale)
TAG            = P.source_tag(SOURCE, SCALE)
STAGE          = args.stage
RUN_MODE       = args.run_mode
MAX_CONCURRENT = args.max_concurrent
GAP_YEARS      = args.gap_years
CRS            = P.CRS
TILE_DEGREES   = args.tile_degrees or G.auto_tile_degrees(SCALE)
GEE_ASSETS     = args.gee_assets
KEEP_INTERMEDIATE = args.keep_intermediate

RANGE_ASSET    = args.range_asset or f'{GEE_ASSETS}/{SPECIES}_iucn_range'
BASE_FOLDER    = f'{GEE_ASSETS}/{SPECIES}_SDM_{TAG}_{TARGET_YEAR}'
TEMP_FOLDER    = f'{BASE_FOLDER}/temp'
METRICS_FOLDER = f'{BASE_FOLDER}/metrics'

# =============================================================================
# 2. TASK MANAGEMENT
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


def export_image(img, desc, asset_id, tile, max_pixels=1e11):
    return ee.batch.Export.image.toAsset(
        image=img, description=desc, assetId=asset_id, region=tile['geom'],
        scale=SCALE, crs=CRS, maxPixels=max_pixels)

# =============================================================================
# 3. YEARS
# =============================================================================

def gap_fill_years():
    """Target year first, then prior years within source coverage. Stops if
    the target year itself is not covered."""
    usable = P.coverage_years(SOURCE)
    if TARGET_YEAR not in usable:
        raise SystemExit(f'--target_year {TARGET_YEAR} is outside {SOURCE} '
                         f'full-year coverage ({usable[0]}-{usable[-1]})')
    return [TARGET_YEAR] + P.prior_years(TARGET_YEAR, GAP_YEARS, usable[0])

# =============================================================================
# 4. STAGE 1: EXPORT ANNUAL STACKS
# =============================================================================

def run_export_annual(tiles, years, max_concurrent=3):
    task_list = []

    for tile in tiles:
        s = gap_fill_state(tile['name'], years)
        if s in ('done', 'metrics'):
            print(f'  Skipping (gap-fill done): {tile["name"]}')
            continue
        # Annuals up to the latest existing merge were consumed and deleted
        for yr in years[s + 1:]:
            asset_id = _annual_asset_id(tile['name'], yr)
            if asset_exists(asset_id):
                print(f'  Skipping (exists): {asset_id.split("/")[-1]}')
                continue
            desc = G.safe_name(f'{SPECIES}_modis_{yr}_{tile["name"]}')
            task = export_image(P.modis_annual_image(yr, tile['geom']), desc, asset_id, tile)
            task_list.append({'task': task, 'desc': desc,
                              'tile': tile['name'], 'year': yr})

    if not task_list:
        print('  All annual assets already exist.')
        return []

    print(f'Submitting {len(task_list)} annual export tasks '
          f'({len(tiles)} tiles x up to {len(years)} years)...')
    completed, failed, total_eecu = submit_and_wait(task_list, max_concurrent)
    print(f'  Completed: {len(completed)}, Failed: {len(failed)}, '
          f'EECU-s: {total_eecu:,.1f}')
    return completed

# =============================================================================
# 5. STAGE 2: GAP-FILL (barrier-synchronized parallel)
#
# All tiles advance through each step together:
#   Step 0:      target year + long-gap flags      -> merged_{y0}
#   Step 1..N:   fill long gaps from prior year yk -> merged_{yk}
#   Final:       linear interpolation              -> gapfilled
# Each step deletes the assets it consumed once its export succeeds,
# unless --keep_intermediate. Kept assets do not affect resume: the
# state below looks for the furthest step first.
#
# Resumable: each tile restarts after its latest existing merged asset.
# Tiles that fail at any step are dropped and reported at the end.
# =============================================================================

def gap_fill_state(tile_name, years):
    """'metrics' / 'done' if past gap-fill, else the index into `years` of the
    latest existing merged asset (-1 if none)."""
    if asset_exists(f'{METRICS_FOLDER}/{tile_name}_metrics'):
        return 'metrics'
    if asset_exists(_gapfilled_asset_id(tile_name)):
        return 'done'
    for k in range(len(years) - 1, -1, -1):
        if asset_exists(_temp_merge_id(tile_name, years[k])):
            return k
    return -1


def run_gap_fill(tiles, years, max_concurrent=3):
    last = len(years) - 1

    # ── Resume state per tile ────────────────────────────────────────────────
    step, pending = {}, []
    skipped_done, skipped_missing = 0, 0
    for tile in tiles:
        s = gap_fill_state(tile['name'], years)
        if s in ('done', 'metrics'):
            skipped_done += 1
            continue
        missing = [y for y in years[s + 1:]
                   if not asset_exists(_annual_asset_id(tile['name'], y))]
        if missing:
            print(f'  Skipping {tile["name"]}: missing annual assets {missing} '
                  f'(run export_annual)')
            skipped_missing += 1
            continue
        step[tile['name']] = s
        pending.append(tile)

    print(f'  Tiles to gap-fill: {len(pending)} | '
          f'Already done: {skipped_done} | '
          f'Missing inputs: {skipped_missing}')
    if not pending:
        print('  Nothing to do.')
        return

    all_eecu, failed_names = 0, set()

    def run_step(task_list):
        nonlocal all_eecu
        if not task_list:
            return
        completed, failed, eecu = submit_and_wait(task_list, max_concurrent)
        all_eecu += eecu
        for item in completed:
            step[item['tile']] = item['new_step']
            if KEEP_INTERMEDIATE:
                continue
            for old_id in item['delete']:
                try: ee.data.deleteAsset(old_id)
                except Exception: pass
        failed_names.update(item['tile'] for item in failed)

    # ── Merge steps ──────────────────────────────────────────────────────────
    for k, yr in enumerate(years):
        active = [t for t in pending if t['name'] not in failed_names
                  and step[t['name']] == k - 1]
        print(f'\n  Step {k}/{last}: merge year {yr} ({len(active)} tiles)...')
        task_list = []
        for tile in active:
            name   = tile['name']
            annual = _annual_asset_id(name, yr)
            if k == 0:
                img, delete = P.first_merge_image(ee.Image(annual), yr), [annual]
            else:
                prev   = _temp_merge_id(name, years[k-1])
                img    = P.merge_prior_image(ee.Image(prev), ee.Image(annual), yr)
                delete = [prev, annual]
            desc = G.safe_name(f'{SPECIES}_modismerge_{yr}_{name}')
            task = export_image(img, desc, _temp_merge_id(name, yr), tile)
            task_list.append({'task': task, 'desc': desc, 'tile': name,
                              'new_step': k, 'delete': delete})
        run_step(task_list)

    # ── Interpolation + gapfilled export ─────────────────────────────────────
    active = [t for t in pending if t['name'] not in failed_names
              and step[t['name']] == last]
    print(f'\n  Final: interp + gapfilled export ({len(active)} tiles)...')
    task_list = []
    for tile in active:
        name   = tile['name']
        merged = _temp_merge_id(name, years[last])
        desc   = G.safe_name(f'{SPECIES}_modisgapfill_{name}')
        task   = export_image(P.interp_image(ee.Image(merged)), desc,
                              _gapfilled_asset_id(name), tile)
        task_list.append({'task': task, 'desc': desc, 'tile': name,
                          'new_step': 'done', 'delete': [merged]})
    run_step(task_list)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_done = len(pending) - len(failed_names)
    print(f'\n  Gap-fill complete: {n_done}/{len(pending)} tiles succeeded, '
          f'{all_eecu:,.1f} EECU-s')
    if failed_names:
        print(f'  WARNING: {len(failed_names)} tiles failed - rerun to retry:')
        for name in sorted(failed_names):
            print(f'    {name}')

# =============================================================================
# 6. STAGE 3: COMPUTE METRICS
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

        img     = ee.Image(gapfill_id)
        metrics = P.metrics_image(img).addBands(img.select('source_year'))
        desc    = G.safe_name(f'{SPECIES}_modismetrics_{tile["name"]}')
        task    = export_image(metrics, desc, metrics_id, tile, max_pixels=1e10)
        task_list.append({'task': task, 'desc': desc, 'tile': tile['name']})

    if not task_list:
        print('  All metric assets already exist.')
        return []

    print(f'Submitting {len(task_list)} metric tasks...')
    completed, failed, total_eecu = submit_and_wait(task_list, max_concurrent)
    print(f'  Completed: {len(completed)}, Failed: {len(failed)}, '
          f'EECU-s: {total_eecu:,.1f}')

    for item in (completed if not KEEP_INTERMEDIATE else []):
        try:
            ee.data.deleteAsset(_gapfilled_asset_id(item['tile']))
            print(f'  Deleted gapfilled: {item["tile"]}')
        except Exception: pass

    return completed

# =============================================================================
# 7. ASSET PATH HELPERS
# =============================================================================

def _annual_asset_id(tile_name, year):
    return f'{TEMP_FOLDER}/{tile_name}_y{year}_annual'

def _temp_merge_id(tile_name, year):
    return f'{TEMP_FOLDER}/{tile_name}_y{year}_merged'

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

def next_stage():
    stages = ['export_annual', 'gap_fill', 'reduce_to_metrics']
    idx    = stages.index(STAGE) if STAGE in stages else 0
    return stages[idx+1] if idx < len(stages)-1 else 'done'

# =============================================================================
# 8. MAIN
# =============================================================================

def main():
    print('=' * 60)
    print(f'SDM MODIS EXPORT')
    print(f'Species: {SPECIES} | Year: {TARGET_YEAR} | Scale: {SCALE} m | '
          f'Tiles: {TILE_DEGREES} deg')
    print(f'Stage: {STAGE} | Mode: {RUN_MODE} | Max concurrent: {MAX_CONCURRENT}')
    print(f'Base folder:    {BASE_FOLDER}')
    print('=' * 60)

    years = gap_fill_years()
    print(f'Years (target, then gap-fill): {years}')
    if not args.estimate_only:
        ensure_folders()

    if RUN_MODE == 'test':
        tiles = G.test_tiles(TILE_DEGREES)
        print(f'\nTest tile: {tiles[0]["name"]}')
    else:
        tiles = G.build_tiles(ee.FeatureCollection(RANGE_ASSET), TILE_DEGREES)

    print()
    print(C.report([C.export(SOURCE, len(tiles), TILE_DEGREES, SCALE,
                             max_concurrent=MAX_CONCURRENT, export_stage=STAGE)],
                   'Cost estimate (all tiles; tiles already done are skipped)'))
    if args.estimate_only:
        return
    print()
    if STAGE == 'export_annual':
        run_export_annual(tiles, years, max_concurrent=MAX_CONCURRENT)
    elif STAGE == 'gap_fill':
        run_gap_fill(tiles, years, max_concurrent=MAX_CONCURRENT)
    elif STAGE == 'reduce_to_metrics':
        run_reduce_to_metrics(tiles, max_concurrent=MAX_CONCURRENT)

    print(f'\nStage complete. Next: --stage {next_stage()}')

if __name__ == '__main__':
    main()
