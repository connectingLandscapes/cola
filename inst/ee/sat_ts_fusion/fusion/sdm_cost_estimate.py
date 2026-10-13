"""
Estimate GEE cost (EECU-hours), task count and wall-clock time for a full
SDM run before launching anything. Runs only cheap metadata queries (point
count/years, range bounds/area); submits no tasks.

  python sdm_cost_estimate.py --species clouded_leopard \\
      --occurrence_asset projects/.../clouded_leopard_occurrences \\
      --predictor_source embed [--scale 30] --target_years 2024 2025

Rates and their provenance live in cost.py.
"""
import argparse

import ee

import cost as C
import grid as G
import predictors as P


def parse_args():
    p = argparse.ArgumentParser(description='SDM pipeline cost estimate')
    p.add_argument('--ee_project',       type=str, default='geersprocessing')
    p.add_argument('--species',          type=str, required=True)
    p.add_argument('--predictor_source', type=str, default='modis',
                   choices=sorted(P.SOURCES))
    p.add_argument('--scale',            type=int, default=None,
                   help='Default: source native resolution')
    p.add_argument('--target_years',     type=int, nargs='+', required=True,
                   help='Years to map wall-to-wall')
    p.add_argument('--tile_degrees',     type=float, default=None,
                   help='Default: chosen from the scale (grid.auto_tile_degrees)')
    p.add_argument('--use_exported_metrics', action='store_true',
                   help='Price the optional MODIS export + prediction from it '
                        '(cheaper only if 2+ models map the same year)')
    p.add_argument('--batch_size',       type=int, default=25)
    p.add_argument('--max_concurrent',   type=int, default=3)
    p.add_argument('--gee_assets',       type=str,
                   default='projects/geersprocessing/assets')
    p.add_argument('--occurrence_asset', type=str, default=None,
                   help='Default: {gee_assets}/{species}_occurrences_batched')
    p.add_argument('--range_asset',      type=str, default=None,
                   help='Default: {gee_assets}/{species}_iucn_range')
    return p.parse_args()


def main():
    args = parse_args()
    ee.Initialize(project=args.ee_project)
    source = args.predictor_source
    scale  = P.resolve_scale(source, args.scale)
    deg    = args.tile_degrees or G.auto_tile_degrees(scale)
    occ    = ee.FeatureCollection(args.occurrence_asset or
                                  f'{args.gee_assets}/{args.species}_occurrences_batched')
    rng    = ee.FeatureCollection(args.range_asset or
                                  f'{args.gee_assets}/{args.species}_iucn_range')

    # Coverage: the same stop-and-tell rule extraction applies
    usable = P.coverage_years(source)
    years  = occ.aggregate_array('year').getInfo()
    bad    = P.out_of_coverage(years, usable)
    bad_targets = [y for y in args.target_years if y not in usable]
    if bad or bad_targets:
        if bad:
            print(P.coverage_message(source, usable, bad))
        if bad_targets:
            print(f'Target year(s) {bad_targets} are outside {source} coverage '
                  f'({usable[0]}-{usable[-1]}).')
        raise SystemExit(1)

    n_tiles, method = G.count_tiles(rng, deg)
    stages = [C.extraction(source, len(years), scale, args.batch_size,
                           args.max_concurrent)]
    exported = args.use_exported_metrics
    if exported and not P.SOURCES[source]['exportable']:
        raise SystemExit(f'--use_exported_metrics: {source} has no export stage')
    if exported:
        stages.append(C.export(source, n_tiles, deg, scale,
                               len(args.target_years), args.max_concurrent))
    stages.append(C.prediction(source, n_tiles, deg, scale,
                               len(args.target_years), args.max_concurrent,
                               exported=exported))

    print(f'Species: {args.species} | Source: {source} | Scale: {scale} m | '
          f'Tiles: {deg} deg')
    print(f'Points: {len(years):,} ({min(years)}-{max(years)}) | '
          f'Range tiles: {n_tiles:,} ({method}) | Target years: {args.target_years}')
    print(f'Max concurrent tasks: {args.max_concurrent}\n')
    print(C.report(stages))
    print('  note: model fitting runs locally (no EECU).')


if __name__ == '__main__':
    main()
