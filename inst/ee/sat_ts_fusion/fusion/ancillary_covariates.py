import ee


def mosaic_gee_ic(ic:ee.ImageCollection = None):
    """Mosaic an image collection and use the projection from the first image in the collection

    :param ic: ee.ImageCollection
    :return: mosaiced ee.Image
    """
    proj = ee.Image(ic.first()).projection()
    orig_band_names = ee.Image(ic.first()).bandNames()
    mosaic_img = ic.mosaic().setDefaultProjection(proj).rename(orig_band_names)
    return mosaic_img




def get_best_dem():
    """Stack several available digital elevation models (ranked from best to worst) and take the first non null value

    :return: DEM image
    """
    # order DEMs (best to worst)
    dem_3dep = ee.Image("USGS/3DEP/10m").select('elevation') # orthometric heights (NAVD88)
    dem_glo30 = ee.Image(mosaic_gee_ic(ee.ImageCollection("COPERNICUS/DEM/GLO30").select(0))).rename('elevation') # orthometric heights (EGM2008)
    dem_alos = ee.Image(mosaic_gee_ic(ee.ImageCollection("JAXA/ALOS/AW3D30/V4_1").select(0))).rename('elevation') # orthometric heights (EGM96)
    dem_nasa = ee.Image("NASA/NASADEM_HGT/001").select('elevation') # orthometric heights (EGM96)
    dem_ic = ee.ImageCollection(ee.List([dem_3dep, dem_glo30, dem_alos, dem_nasa]))
    best_dem = ee.Image(dem_ic.reduce(ee.Reducer.firstNonNull())).reproject(crs='EPSG:4326', crsTransform=None, scale=30).rename('elevation')
    return best_dem



def degrees_to_radians(img:ee.Image = None):
    return img.multiply(3.14159265359).divide(180)



def make_covar_topo_stack(dem:ee.Image = get_best_dem()):
    """Make a multiband image with elevation, slope, and aspect
    :param dem: ee.Image digital elevation model
    :return: an image with these bands: elevation, slope (in radians), aspect (in radians)
    """
    slope_rad = ee.Image(degrees_to_radians(ee.Terrain.slope(dem)))
    aspect_rad = ee.Image(degrees_to_radians(ee.Terrain.aspect(dem)))
    topo_img_stack = (ee.Image.cat([dem, slope_rad, aspect_rad])
                              .rename(['elevation', 'slope_rad', 'aspect_rad']))
    return topo_img_stack.set({'temporal_res': 'static',
                               'nominal_spatial_res': 30})





def get_google_sat_embed(year):
    """Make an annual mosaic of Googles Satellite Embedding dataset

    :param year: year of interest, as integer
    :return: a mosaic image corresponding to a year
    """
    # filter by year of interest and mosaic images in the collection
    sat_embed_y = (ee.ImageCollection("GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL")
                     .filter(ee.Filter.calendarRange(year, year, 'year'))
                     .mosaic())
    # rename bands
    orig_band_names = ee.Image(sat_embed_y).bandNames()
    new_band_names = orig_band_names.map(lambda x: ee.String('gse_').cat(ee.String(x).toLowerCase()))
    return (sat_embed_y.rename(new_band_names)
                       .set({'temporal_res': 'annual',
                             'nominal_spatial_res': 10}))
