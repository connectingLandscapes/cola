import ee
# initialize with Earth Engine project id
gee_project_name = 'ee-gedibio'
ee.Initialize(project=gee_project_name)

from sat_ts_fusion.imagery import lsat_utils



##########
# INPUTS
##########
# a geometry for the region of interest
region_geom = (ee.Feature(ee.FeatureCollection('users/pb463/S2L/Sonoma_cty_v2_PBcleaned').first())
                 .geometry())

# a grid covering the region of interest geometry
proj_epsg = 'EPSG:3310'
grid_fc = region_geom.coveringGrid(proj=proj_epsg, scale=25000)

# CCDC run parameters
start_year = 1999
end_year = 2024
start_doy = 1
end_doy = 366
max_cloud_cover_land = 50
min_sun_elev = 40
bands = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
sensors_meta = 'landsat_45789'

# export asset folder
asset_dir = 'projects/ee-gedibio/assets/ccdc/results/python_test/'



##########
# PROCESSING
##########

# Landsat Collection 2 surface reflectance image collections
l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")

# Merge Landsat collections
l4to9 = l5.merge(l7).merge(l8).merge(l9)

ccdcParams ={
  'breakpointBands':['green','red','NIR','SWIR1','SWIR2'],
  'tmaskBands' : ['green', 'SWIR1'],
  'minObservations': 6,
  'chiSquareProbability': 0.995,
  'minNumOfYearsScaler': 1.33,
  'lambda': 0.005,
  'maxIterations' : 10000,
  'dateFormat' : 1
};

runParams = {
    "start_year": start_year,
    "end_year": end_year,
    "start_doy": start_doy,
    "end_doy": end_doy,
    "min_sun_elev": min_sun_elev,
    "max_cloud_cover_land": max_cloud_cover_land,
    "lsat_sensors": sensors_meta,
    "includeSLCOffL7": True
}

# loop over grids and run/export CCDC
n_grids=ee.Number(grid_fc.size()).getInfo()
for i in range(0, n_grids):
    g_geom = ee.Feature(grid_fc.toList(1000).get(i)).geometry()

    # Filter the merged collection spatially and temporally
    # and mask low/high SR, exclude clouds, and calculate spectral indices
    l4to9_p_g = (l4to9.filterBounds(g_geom)
               .filter(ee.Filter.calendarRange(start_year, end_year, 'year'))
               .filter(ee.Filter.calendarRange(start_doy, end_doy, 'day_of_year'))
               .filter(ee.Filter.gte('SUN_ELEVATION', min_sun_elev))
               .filter(ee.Filter.lte('CLOUD_COVER_LAND', max_cloud_cover_land))
               .map(lsat_utils.l4to9_c2_scaleoff)
               .map(lsat_utils.l4to9_c2_qa_mask_clouds)
               .map(lsat_utils.l4to9_c2_indices)
               .select(bands))

    ccdcParams['collection'] = l4to9_p_g

    ccdc_array = (ee.Algorithms.TemporalSegmentation.Ccdc(**ccdcParams)
                   .clip(g_geom)
                   .set('grid_i', i)
                   .set(runParams)
                   .set(ccdcParams))

    task_name = 'lsatc2sr_ccdc_grid_'+str(i)

    task = ee.batch.Export.image.toAsset(
        image=ccdc_array,
        description=task_name,
        assetId=f"{asset_dir}{task_name}",
        pyramidingPolicy={".default": "sample"},
        scale=30,
        crs='EPSG:4326',
        maxPixels=1e13
    )
    #task.start()
