"""
Tile grid shared by sdm_modis_export.py and sdm_wall_to_wall.py.
Tile names must match exactly between the two, so both use these helpers.
No Earth Engine calls at import time.
"""
import math

# Tile size follows the analysis scale so each tile stays near
# TILE_TARGET_PX pixels per side: 2-degree MODIS tiles at 250 m (~890 px)
# are the size actually run end-to-end (June 2026, 2026-10-06 test).
TILE_TARGET_PX      = 1000
TILE_DEGREE_CHOICES = [0.1, 0.25, 0.5, 1, 2]   # 2 is the largest tested
M_PER_DEGREE        = 111320


def auto_tile_degrees(scale):
    want = scale * TILE_TARGET_PX / M_PER_DEGREE
    fits = [d for d in TILE_DEGREE_CHOICES if d <= want]
    return max(fits) if fits else TILE_DEGREE_CHOICES[0]


def safe_name(s):
    """GEE task descriptions and asset IDs disallow '+' (and '.')."""
    return s.replace('+', 'p').replace('.', 'd')


def tile_name(lon, lat, tile_degrees):
    """Integer-degree tiles keep the original format (tile_-084_p08);
    fractional tiles use thousandths of a degree (tile_-084d250_p08d500)."""
    if float(tile_degrees).is_integer():
        return safe_name(f'tile_{int(round(lon)):+04d}_{int(round(lat)):+03d}')
    return safe_name(f'tile_{lon:+08.3f}_{lat:+07.3f}')


def grid_cells(bounds, tile_degrees):
    """Cells [(lon, lat)] of lower-left corners covering bounds
    (min_lon, min_lat, max_lon, max_lat). Integer arithmetic in thousandths
    of a degree so fractional tiles never accumulate float error."""
    step = int(round(tile_degrees * 1000))
    lo_x = math.floor(bounds[0] * 1000 / step) * step
    lo_y = math.floor(bounds[1] * 1000 / step) * step
    hi_x = math.ceil(bounds[2] * 1000 / step) * step
    hi_y = math.ceil(bounds[3] * 1000 / step) * step
    return [(x / 1000, y / 1000)
            for y in range(lo_y, hi_y, step) for x in range(lo_x, hi_x, step)]


def neighbour_cells(lon, lat, tile_degrees):
    """The cell and its 8 neighbours, same thousandths arithmetic."""
    step = int(round(tile_degrees * 1000))
    x0, y0 = int(round(lon * 1000)), int(round(lat * 1000))
    return [((x0 + dx) / 1000, (y0 + dy) / 1000)
            for dy in (-step, 0, step) for dx in (-step, 0, step)]


def build_tiles(range_fc, tile_degrees):
    """Tiles intersecting the range: [{'geom', 'name', 'lon', 'lat'}]."""
    import ee
    coords = range_fc.geometry().bounds().coordinates().get(0).getInfo()
    lons, lats = [c[0] for c in coords], [c[1] for c in coords]
    cells = grid_cells((min(lons), min(lats), max(lons), max(lats)), tile_degrees)
    geoms = [ee.Geometry.Rectangle([x, y, x + tile_degrees, y + tile_degrees])
             for x, y in cells]
    hits = ee.List([range_fc.geometry().intersects(g, ee.ErrorMargin(1000))
                    for g in geoms]).getInfo()          # one call, not one per cell
    tiles = [{'geom': g, 'name': tile_name(x, y, tile_degrees), 'lon': x, 'lat': y}
             for (x, y), g, hit in zip(cells, geoms, hits) if hit]
    print(f'Built {len(tiles)} tiles ({tile_degrees} deg) over range')
    return tiles


def test_tiles(tile_degrees):
    """Fixed Costa Rica test tile(s), on the grid so neighbours work."""
    import ee
    d = tile_degrees
    if float(d).is_integer():          # existing MODIS test assets
        lon, lat, name = -84, 9, f'test_costa_rica_{safe_name(str(d))}deg'
    else:
        # (-84, 9) is open Pacific at 0.1-0.5 deg (0 valid pixels,
        # 2026-10-07); (-84, 10) is inland (Braulio Carrillo NP).
        lon, lat, name = -84, 10, f'test_costa_rica_inland_{safe_name(str(d))}deg'
    return [{'geom': ee.Geometry.Rectangle([lon, lat, lon + d, lat + d]),
             'name': name, 'lon': lon, 'lat': lat}]


EXACT_COUNT_MAX_CELLS = 50000   # grid built and tested server-side


def count_tiles(range_fc, tile_degrees):
    """(n_tiles, method). Exact when the bounding grid is small enough to
    test every cell, otherwise range area / tile area (an undercount of
    edge tiles, which are mostly masked and cheap)."""
    import ee
    coords = range_fc.geometry().bounds().coordinates().get(0).getInfo()
    lons, lats = [c[0] for c in coords], [c[1] for c in coords]
    cells = grid_cells((min(lons), min(lats), max(lons), max(lats)), tile_degrees)
    if len(cells) <= EXACT_COUNT_MAX_CELLS:
        d = tile_degrees
        xs = sorted({x for x, _ in cells})
        ys = sorted({y for _, y in cells})
        def row(y):
            return ee.List(xs).map(lambda x: ee.Feature(ee.Geometry.Rectangle(
                [x, y, ee.Number(x).add(d), ee.Number(y).add(d)], None, False)))
        grid = ee.FeatureCollection(ee.List(ys).map(row).flatten())
        return grid.filterBounds(range_fc.geometry()).size().getInfo(), 'exact'
    area_km2 = range_fc.geometry().area(1000).getInfo() / 1e6
    tile_km2 = (tile_degrees * M_PER_DEGREE / 1000) ** 2
    return math.ceil(area_km2 / tile_km2), f'area-based ({len(cells):,} grid cells)'
