"""
SDM Predictor Extraction
========================
Extracts predictor variables at species occurrence locations. Predictor
construction lives in predictors.py and is shared with wall-to-wall
prediction. All outputs are named after SPECIES, the predictor source/scale
and MODEL_ID so multiple runs never overwrite each other.

Occurrence years must fall inside the source's coverage; otherwise the run
stops before submitting any task.

Usage:
  python sdm_extraction.py --species puma --model_id puma_modis_m3 \
      --run_mode full [--predictor_source modis] [--scale 250] ...
"""

import ee, time, os, sys, argparse
from collections import Counter
from datetime import datetime
import predictors as P
import cost as C

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM Predictor Extraction')

    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',      type=str, default='geersprocessing',
                   help='GEE project ID')
    p.add_argument('--species',         type=str, required=True,
                   help='Species name, e.g. puma')
    p.add_argument('--model_id',        type=str, required=True,
                   help='Unique model run ID, e.g. puma_modis_m3')
    p.add_argument('--predictor_source', type=str, default='modis',
                   choices=sorted(P.SOURCES))
    p.add_argument('--scale',           type=int, default=None,
                   help='Analysis scale in metres (default: source native; '
                        'may only be coarser)')

    # ── Run control ───────────────────────────────────────────────────────────
    p.add_argument('--estimate_only',   action='store_true',
                   help='Print the cost estimate and exit')
    p.add_argument('--run_mode',        type=str, default='single',
                   choices=['single', 'full', 'resume'],
                   help='Execution mode')
    p.add_argument('--batch_year',      type=int, default=2015,
                   help='Year for single-batch test run')
    p.add_argument('--batch_num',       type=int, default=1,
                   help='Batch number for single run: the n-th group of '
                        '--batch_size points in --batch_year, as full mode '
                        'builds them (not the occurrence "batch" property)')
    p.add_argument('--batch_size',      type=int, default=25,
                   help='Points per batch')
    p.add_argument('--max_concurrent',  type=int, default=5,
                   help='Max concurrent GEE tasks')

    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--occurrence_asset', type=str, default=None,
                   help='GEE asset path for occurrence points '
                        '(default: {gee_assets}/{species}_occurrences_batched)')
    # --column_train is the name the cola Shiny GUI passes.
    p.add_argument('--target_col', '--column_train', dest='target_col',
                   type=str, default='presence',
                   help='Response column on the occurrence asset (same as '
                        'fitting --target_col). Only this and year are '
                        'carried into the training tables, so other '
                        'properties cannot become predictors.')
    p.add_argument('--gee_assets',      type=str,
                   default='projects/geersprocessing/assets',
                   help='Root GEE assets folder (no trailing slash)')
    p.add_argument('--working_dir',     type=str, default='.',
                   help='Local folder for the resume log (created if missing)')

    # ── Temporal ──────────────────────────────────────────────────────────────
    p.add_argument('--gap_years',       type=int, default=P.DEFAULT_GAP_YEARS,
                   help='Prior years used for gap-fill (stamped on the points; '
                        'prediction reuses it via the features asset)')

    return p.parse_args()

args = parse_args()
ee.Initialize(project=args.ee_project)

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

SPECIES          = args.species
MODEL_ID         = args.model_id
SOURCE           = args.predictor_source
SCALE            = P.resolve_scale(SOURCE, args.scale)
TAG              = P.source_tag(SOURCE, SCALE)
RUN_MODE         = args.run_mode
BATCH_YEAR       = args.batch_year
BATCH_NUM        = args.batch_num
BATCH_SIZE       = args.batch_size
MAX_CONCURRENT   = args.max_concurrent
GAP_YEARS        = args.gap_years
GEE_ASSETS       = args.gee_assets
TARGET_COL       = args.target_col
WORKING_DIR      = args.working_dir

OCCURRENCE_ASSET = (args.occurrence_asset or
                    f'{GEE_ASSETS}/{SPECIES}_occurrences_batched')
EXPORT_FOLDER    = f'{GEE_ASSETS}/{SPECIES}_SDM_{TAG}_exports_{MODEL_ID}'
LOG_FILE         = os.path.join(WORKING_DIR,
                                f'completed_batches_{SPECIES}_{TAG}_{MODEL_ID}.txt')

# Imagery is only loaded within this distance of a batch's points; must
# exceed the largest kernel radius so neighbourhoods are complete.
REGION_BUFFER_M  = 1.5 * max(P.KERNEL_RADII_M)

# =============================================================================
# 2. TASK MANAGEMENT
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
                n = ee.FeatureCollection(item['asset_id']).size().getInfo()
                if n != item['n_expected']:
                    # A short table is a silent failure; leave it out of the
                    # log so --run_mode resume redoes it.
                    print(f"SHORT {item['description']} -- {n}/{item['n_expected']} "
                          f"rows; not logged as complete")
                else:
                    print(f"OK {item['description']} -- {n} rows, {eecu:.1f} EECU-s")
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
# 3. COVERAGE CHECK
# =============================================================================

def occurrence_years(occurrences):
    """One int year per point; stops on missing or fractional years."""
    try:
        return P.whole_years(occurrences.aggregate_array('year').getInfo())
    except ValueError as e:
        sys.exit(str(e))


def check_target_col(occurrences):
    """Stop before any task if the response column is missing."""
    props = occurrences.first().propertyNames().getInfo()
    if TARGET_COL not in props:
        sys.exit(f'Occurrence asset has no {TARGET_COL!r} property (set '
                 f'--target_col/--column_train). Properties: '
                 f'{sorted(p for p in props if not p.startswith("system:"))}')


def check_coverage(occurrences):
    """Stop before any task if occurrence years fall outside the source's
    full-year coverage. Returns the usable years."""
    years  = occurrence_years(occurrences)
    usable = P.coverage_years(SOURCE)
    bad    = P.out_of_coverage(years, usable)
    if bad:
        print('\n' + P.coverage_message(SOURCE, usable, bad))
        sys.exit(1)
    print(f'Coverage OK: {len(years)} points, years {min(years)}-{max(years)} '
          f'within {SOURCE} {usable[0]}-{usable[-1]}')
    return usable

# =============================================================================
# 4. EXTRACTION
# =============================================================================

def sample_points(pts, year, first_usable):
    """Predictor values at `pts` (all from `year`), one sampleRegions call so
    nearby points share computed image tiles."""
    region  = pts.geometry().bounds().buffer(REGION_BUFFER_M)
    metrics = P.source_predictors(SOURCE, year,
                                 P.prior_years(year, GAP_YEARS, first_usable), region)
    stack   = P.assemble_stack(SOURCE, metrics, year, SCALE).unmask(P.SENTINEL)
    return (stack.sampleRegions(collection=pts, properties=['year', TARGET_COL],
                                scale=SCALE, projection=ee.Projection(P.CRS),
                                tileScale=8, geometries=True)  # toAsset rejects null geometry
                 .map(lambda f: f.set({'predictor_source': SOURCE, 'scale': SCALE,
                                       'gap_years': GAP_YEARS})))


def build_batches(occurrences, batch_size=25):
    years = sorted(set(occurrence_years(occurrences)))
    print(f'Unique years: {years}')
    batches = []
    for yi, year in enumerate(years):
        print(f'  Fetching IDs for {year} ({yi+1}/{len(years)})...')
        ids = sorted(occurrences.filter(ee.Filter.eq('year', year))
                                .aggregate_array('system:index').getInfo())
        for i in range(-(-len(ids) // batch_size)):
            batches.append({'year': year, 'batch': i+1, 'offset': i*batch_size,
                            'ids': ids[i*batch_size:(i+1)*batch_size]})
    sizes = [len(b['ids']) for b in batches]
    print(f'Total batches: {len(batches)}, pts: min={min(sizes)} max={max(sizes)}')
    return batches


def batch_points(occurrences, year, offset, batch_size):
    """One batch as a slice of the year's points sorted by system:index
    (same order as build_batches' ids). Do NOT use
    Filter.inList('system:index', ids): on such a collection size()/limit()
    see the points but geometry(), aggregate_array() and sampleRegions() see
    none (verified 2026-10-07: sampleRegions returned 0 of 25 rows)."""
    year_pts = occurrences.filter(ee.Filter.eq('year', year)).sort('system:index')
    return ee.FeatureCollection(year_pts.toList(batch_size, offset))


def batch_desc(year, batch_num):
    return f'SDM_{SOURCE.upper()}_{SPECIES}_{year}_b{batch_num:02d}'


def build_extraction_task(pts, year, batch_num, first_usable, n_expected):
    desc     = batch_desc(year, batch_num)
    asset_id = f'{EXPORT_FOLDER}/{desc}'
    task = ee.batch.Export.table.toAsset(
        collection=sample_points(pts, year, first_usable),
        description=desc, assetId=asset_id)
    return {'task': task, 'description': desc, 'asset_id': asset_id,
            'n_expected': n_expected}

# =============================================================================
# 5. OUTPUT FOLDER SETUP
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
# 6. MAIN
# =============================================================================

if __name__ == '__main__':
    occ = ee.FeatureCollection(OCCURRENCE_ASSET)
    print(f'Species:          {SPECIES}')
    print(f'Model ID:         {MODEL_ID}')
    print(f'Source / scale:   {SOURCE} / {SCALE} m')
    print(f'Export folder:    {EXPORT_FOLDER}')
    print(f'Response column:  {TARGET_COL}')
    print(f'Log file:         {LOG_FILE}')

    check_target_col(occ)
    usable = check_coverage(occ)
    if not args.estimate_only:
        ensure_folder()
        os.makedirs(WORKING_DIR, exist_ok=True)

    if RUN_MODE == 'single':
        print(f'\nSingle-batch: year={BATCH_YEAR}, batch={BATCH_NUM}')
        # Same slice full/resume mode builds for this batch number (not the
        # asset's own 'batch' property, which may use another batch size).
        pts = batch_points(occ, BATCH_YEAR, (BATCH_NUM - 1) * BATCH_SIZE, BATCH_SIZE)
        n_pts = pts.size().getInfo()
        if n_pts == 0:
            sys.exit(f'No points in year {BATCH_YEAR} batch {BATCH_NUM} '
                     f'(batch_size {BATCH_SIZE}).')
        print(f'Points: {n_pts}')
        print(C.report([C.extraction(SOURCE, n_pts, SCALE, BATCH_SIZE, 1)]))
        if args.estimate_only:
            sys.exit(0)
        item = build_extraction_task(pts, BATCH_YEAR, BATCH_NUM, usable[0], n_pts)
        item['task'].start()
        print(f'Task started: {item["description"]}')
        wait_for_tasks([item['task']])
        n = ee.FeatureCollection(item['asset_id']).size().getInfo()
        print(f'Rows written: {n}/{n_pts}' + ('' if n == n_pts else '  <-- SHORT'))

    else:
        is_resume  = (RUN_MODE == 'resume')
        print(f'Mode: {"Resume" if is_resume else "Full"} | '
              f'batch_size={BATCH_SIZE}, max_concurrent={MAX_CONCURRENT}')

        batches = build_batches(occ, BATCH_SIZE)

        if is_resume:
            already_done = load_completed(LOG_FILE)
            batches = [b for b in batches
                       if batch_desc(b['year'], b['batch']) not in already_done]
            print(f'Resume: {len(batches)} batches remaining.')
        elif not args.estimate_only:
            try:
                os.remove(LOG_FILE)
                print(f'Cleared log: {LOG_FILE}')
            except FileNotFoundError:
                pass

        if not batches:
            print('All batches complete.')
            sys.exit(0)

        n_left = sum(len(b['ids']) for b in batches)
        print(C.report([C.extraction(SOURCE, n_left, SCALE, BATCH_SIZE, MAX_CONCURRENT)],
                       f'Cost estimate ({len(batches)} batches to run)'))
        if args.estimate_only:
            sys.exit(0)
        print('Building task list...')
        task_list = []
        for b in batches:
            pts = batch_points(occ, b['year'], b['offset'], BATCH_SIZE)
            task_list.append(build_extraction_task(pts, b['year'], b['batch'],
                                                   usable[0], len(b['ids'])))
        print(f'{len(task_list)} tasks. Starting submission...')
        completed = submit_batch_queue(task_list, MAX_CONCURRENT,
                                       log_file=LOG_FILE)
        total_eecu = sum(t['task'].status().get('batch_eecu_usage_seconds', 0)
                         for t in completed)
        print(f'\nDone: {len(completed)}/{len(task_list)} batches, '
              f'{total_eecu:,.0f} EECU-s')
        if len(completed) < len(task_list):
            print(f'Set --run_mode resume and rerun.')
