"""
Cost estimates (GEE batch EECU-seconds, task counts, wall-clock) for the
SDM pipeline, built from measured task costs. No Earth Engine calls.

Every rate below is a measurement; keep its provenance with it and replace
it when better numbers arrive. Rates measured at one scale are extrapolated
to others by output pixel count, which is flagged in the report. HLS
prediction is the exception: its cost per pixel rises at coarser scales,
so it uses rates measured at 30 and 300 m (PREDICT_PER_MPX_BY_SCALE).
"""
import math

M_PER_DEGREE = 111320

# GEE noncommercial monthly batch quotas (EECU-hours), from
# developers.google.com/earth-engine/guides/noncommercial_tiers, 2026-10-08.
# Concurrent-task limits per tier are not published there; the user's
# --max_concurrent stands in for them.
TIERS = {'community': 150, 'contributor': 1_000, 'partner': 100_000}
QUEUE_MAX_TASKS = 3000   # READY tasks per project (guides/usage)

# ── Measured rates ───────────────────────────────────────────────────────────
# Extraction: EECU-s per occurrence point (25-point batches, 2019 puma test
# asset, points spread across the Americas, 2026-10-07/08).
#   embed 10 m : 9,348 EECU-s / 421 pts in 17 tasks (412-660 per batch)
#   modis 250 m: 957 EECU-s / 25 pts (gap-fills 3 years in memory per batch)
#   hls_l30 30 m: 50,483 EECU-s / 421 pts in 17 tasks (2,158-4,617 per batch)
#   hls_s30 30 m: 9,079 and 9,202 EECU-s / 25 pts (same batch twice; the
#                 second loaded 31x fewer images with no saving)
EXTRACT_PER_POINT = {'embed': 22.0, 'modis': 38.0, 'hls_l30': 120.0, 'hls_s30': 365.0}
EXTRACT_TASK_MIN  = 1.2          # median task runtime, minutes

# Wall-to-wall prediction: EECU-s per output megapixel (20 selected
# features, classify + ocean mask).
#   modis 250 m, 2 deg Costa Rica test tile (0.79 Mpx, 27% masked), model
#     m2 (1 MODIS feature of 17), 2026-10-09: metrics in memory (default)
#     221 -> 279 / Mpx; all-land Tocantins tile 318. From exported metric
#     tiles: 74 (2026-10-06), 87 the same day as the in-memory run.
#   embed 10 m, 0.1 deg inland tile (1.24 Mpx): 128 -> 103 / Mpx
#   embed 10 m, 0.1 deg all-ocean tile: 40 (masked tiles are cheaper)
PREDICT_PER_MPX = {'modis': 279.0, 'embed': 103.0, 'hls_l30': 305.0, 'hls_s30': 756.0}
# Prediction from exported metric tiles (sdm_modis_export.py), excluding
# the export itself (EXPORT_PER_MPX).
PREDICT_PER_MPX_EXPORTED = {'modis': 93.0}
# HLS is computed in memory, so its rate depends on how many HLS predictors
# the model selects. Measured 2026-10-08 on one 0.25 deg 30 m tile (0.86
# Mpx): model puma_hlsl30_t2019 (1 HLS feature of 15) 305; all 32 HLS
# predictors 920; all 8 at 30 m 479; whole 172-band stack 2,356.
PREDICT_HLS_ALL_FACTOR = 3.0
# HLS cost per output pixel rises at coarser output scales: each coarse
# pixel still reads every scene beneath it, so pixel-count extrapolation
# from 30 m underestimates. Measured 2026-10-08 at 300 m on 2 deg tiles
# (0.55 Mpx): the same 15 features (30 m base renamed to 300 m) through
# the 30 m classifier, 773 EECU-s on an all-land Tocantins tile (1,403 per
# Mpx); 727 on the half-ocean Costa Rica test tile (cost tracks area, not
# land pixels); all 32 HLS predictors 1,383 (2,510 per Mpx, 1.8x).
# Between measured scales the rate is interpolated in log-log space.
# S30, 2026-10-08: the L30 model's 15 features computed from S30 on the
# same 0.25 deg 30 m tile, classified with the L30 forest: 651 EECU-s
# (756 per Mpx, 2.5x L30), 6.8 min. Its 300 m rate is not measured; it
# is assumed to rise from 30 m by the same factor as L30's (4.6x).
PREDICT_PER_MPX_BY_SCALE = {'hls_l30': {30: 305.0, 300: 1403.0},
                            'hls_s30': {30: 756.0, 300: 756.0 * 1403.0 / 305.0}}
PREDICT_SCALE_ASSUMED = {'hls_s30': 30}   # only this scale measured
# Effective minutes per task in a full run. Wall-clock is the least reliable
# estimate: the June 2026 MODIS run (673 tiles, 3 concurrent, 46.7 h) averaged
# 12.5 min per task, while single test tasks finish in 1-2.4 min.
PREDICT_TASK_MIN = {'modis': 12.5, 'embed': 2.4, 'hls_l30': 5.0, 'hls_s30': 6.8}
PREDICT_TASK_MIN_MEASURED = {'modis': 'full run', 'embed': 'single test task',
                             'hls_l30': 'single test task', 'hls_s30': 'single test task'}

# MODIS staged export (annual composites, gap-fill, metrics), per 2 deg tile
# per target year at 250 m: 104 + 62 + 65 + 23 = 254 EECU-s in 8 tasks of
# ~1 min (2026-10-06 test tile, ~50% ocean, so inland tiles may cost more).
# Optional since 2026-10-09: export + exported prediction (254 + 74-87)
# beats in-memory prediction (221) once 2+ models map the same species and
# year from one export.
EXPORT_PER_MPX   = {'modis': 254.0 / 0.79}
EXPORT_TASKS_PER_TILE = {'modis': 8}
EXPORT_STAGES = {                 # stage: (share of EECU, tasks per tile)
    'export_annual':     (104 / 254, 3),
    'gap_fill':          (127 / 254, 4),
    'reduce_to_metrics': (23 / 254, 1),
}
EXPORT_TASK_MIN  = 1.2

MEASURED_SCALE = {'modis': 250, 'embed': 10, 'hls_l30': 30, 'hls_s30': 30}
UNCERTAINTY = 0.5    # report estimates as -50% / +50%


def tile_mpx(tile_degrees, scale):
    """Output megapixels in one tile (degrees -> metres at the equator)."""
    side = tile_degrees * M_PER_DEGREE / scale
    return side * side / 1e6


def tiles_for_area(area_km2, tile_degrees):
    """Approximate tile count when the exact grid is too large to list."""
    tile_km2 = (tile_degrees * M_PER_DEGREE / 1000) ** 2
    return math.ceil(area_km2 / tile_km2)


def _unmeasured(name, n_tasks, task_min, max_concurrent, source):
    """Stage with tasks/wall-clock but no EECU figure (eecu_s None)."""
    return stage(name, None, n_tasks, task_min, max_concurrent,
                 [f'no measured EECU rate for {source} yet'])


def stage(name, eecu_s, n_tasks, task_min, max_concurrent, notes=()):
    waves = math.ceil(n_tasks / max(1, max_concurrent))
    return {'stage': name, 'eecu_s': eecu_s, 'tasks': n_tasks,
            'hours': waves * task_min / 60, 'notes': list(notes)}


def extraction(source, n_points, scale, batch_size=25, max_concurrent=3):
    if source not in EXTRACT_PER_POINT:
        return _unmeasured('extraction', math.ceil(n_points / batch_size),
                           EXTRACT_TASK_MIN, max_concurrent, source)
    notes = []
    if scale != MEASURED_SCALE[source]:
        notes.append(f'per-point rate measured at {MEASURED_SCALE[source]} m only')
    return stage('extraction', n_points * EXTRACT_PER_POINT[source],
                 math.ceil(n_points / batch_size), EXTRACT_TASK_MIN,
                 max_concurrent, notes)


def export(source, n_tiles, tile_degrees, scale, n_years=1, max_concurrent=3,
           export_stage=None):
    """Staged metric export (sources with an export stage only); all
    stages, or one of EXPORT_STAGES."""
    if source not in EXPORT_PER_MPX:
        raise ValueError(f'{source} has no export stage')
    share, per_tile = ((1.0, EXPORT_TASKS_PER_TILE[source]) if export_stage is None
                       else EXPORT_STAGES[export_stage])
    mpx = n_tiles * tile_mpx(tile_degrees, scale) * n_years
    notes = []
    if scale != MEASURED_SCALE[source]:
        notes.append(f'extrapolated by pixel count from {MEASURED_SCALE[source]} m')
    return stage(export_stage or 'export', mpx * EXPORT_PER_MPX[source] * share,
                 n_tiles * n_years * per_tile, EXPORT_TASK_MIN, max_concurrent, notes)


def _predict_rate(source, scale):
    """EECU-s per output Mpx at `scale`, and any caveat notes."""
    by_scale = PREDICT_PER_MPX_BY_SCALE.get(source)
    if not by_scale:
        if scale != MEASURED_SCALE[source]:
            return PREDICT_PER_MPX[source], [
                f'extrapolated by pixel count from {MEASURED_SCALE[source]} m']
        return PREDICT_PER_MPX[source], []
    pts = sorted(by_scale.items())
    measured = PREDICT_SCALE_ASSUMED.get(source)
    if measured is not None and scale != measured:
        rate, notes = _interp_rate(pts, scale)
        return rate, notes + [f'only {measured} m measured for {source}; other '
                              f'scales assume the hls_l30 scale factor']
    return _interp_rate(pts, scale)


def _interp_rate(pts, scale):
    by_scale = dict(pts)
    if scale in by_scale:
        return by_scale[scale], []
    if scale < pts[0][0] or scale > pts[-1][0]:
        s0, r0 = pts[0] if scale < pts[0][0] else pts[-1]
        return r0, [f'rate measured at {pts[0][0]}-{pts[-1][0]} m only; '
                    f'using the {s0} m rate']
    for (s0, r0), (s1, r1) in zip(pts, pts[1:]):
        if s0 < scale < s1:
            f = math.log(scale / s0) / math.log(s1 / s0)
            return r0 * (r1 / r0) ** f, [f'rate interpolated between {s0} and {s1} m']


def prediction(source, n_tiles, tile_degrees, scale, n_years=1, max_concurrent=3,
               exported=False):
    """exported: read metrics from export tiles instead of computing them."""
    if source not in PREDICT_PER_MPX:
        return _unmeasured('wall-to-wall', n_tiles * n_years, PREDICT_TASK_MIN['modis'],
                           max_concurrent, source)
    mpx = n_tiles * tile_mpx(tile_degrees, scale) * n_years
    if exported:
        if source not in PREDICT_PER_MPX_EXPORTED:
            raise ValueError(f'{source} has no export stage')
        rate, notes = PREDICT_PER_MPX_EXPORTED[source], []
        if scale != MEASURED_SCALE[source]:
            notes.append(f'extrapolated by pixel count from {MEASURED_SCALE[source]} m')
    else:
        rate, notes = _predict_rate(source, scale)
    notes = ['coastal/partial tiles cost less than full land tiles'] + notes
    if source.startswith('hls'):
        notes.append(f'rate is for a model with 1 HLS predictor; up to '
                     f'~{PREDICT_HLS_ALL_FACTOR:.0f}x if it selects many HLS predictors')
    if PREDICT_TASK_MIN_MEASURED[source] != 'full run':
        notes.append('wall-clock from a single test task; full runs can be '
                     'several times slower (MODIS June 2026: ~5x)')
    return stage('wall-to-wall', mpx * rate,
                 n_tiles * n_years, PREDICT_TASK_MIN[source], max_concurrent, notes)


def report(stages, title='Cost estimate'):
    """Plain-text table plus totals against each tier's monthly quota."""
    lines = [title, '-' * len(title),
             f'{"Stage":<14}{"EECU-h":>12}{"range":>22}{"tasks":>10}{"wall-clock":>13}']
    for s in stages:
        if s['eecu_s'] is None:
            lines.append(f'{s["stage"]:<14}{"unmeasured":>12}{"":>22}{s["tasks"]:>10,}'
                         f'{_fmt_hours(s["hours"]):>13}')
            continue
        h = s['eecu_s'] / 3600
        rng = f'{h * (1 - UNCERTAINTY):,.1f}-{h * (1 + UNCERTAINTY):,.1f}'
        lines.append(f'{s["stage"]:<14}{h:>12,.1f}{rng:>22}{s["tasks"]:>10,}'
                     f'{_fmt_hours(s["hours"]):>13}')
    total_h = sum(s['eecu_s'] or 0 for s in stages) / 3600
    total_wall = sum(s['hours'] for s in stages)
    lines.append(f'{"TOTAL":<14}{total_h:>12,.1f}'
                 f'{f"{total_h * (1 - UNCERTAINTY):,.1f}-{total_h * (1 + UNCERTAINTY):,.1f}":>22}'
                 f'{sum(s["tasks"] for s in stages):>10,}{_fmt_hours(total_wall):>13}')
    if any(s['eecu_s'] is None for s in stages):
        lines.append('TOTAL excludes unmeasured stages.')
    lines.append('Share of monthly batch quota: ' + ', '.join(
        f'{t} {100 * total_h / q:,.1f}%' for t, q in TIERS.items()))
    for s in stages:
        for n in s['notes']:
            lines.append(f'  note ({s["stage"]}): {n}')
    return '\n'.join(lines)


def _fmt_hours(h):
    if h < 1:
        return f'{h * 60:.0f} min'
    if h < 48:
        return f'{h:.1f} h'
    return f'{h / 24:.1f} days'
