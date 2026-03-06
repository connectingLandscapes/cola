# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Cinnects EE from command line
@author: ig299
"""

#%%
# IMPORTS
import ee
from sat_ts_fusion.imagery import lsat_utils

def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    param1 = sys.argv[1] # project name
    param2 = sys.argv[2] # ee point layer
    param3 = sys.argv[3] # EPSG
    param4 = sys.argv[4] # grid scale
    param5 = sys.argv[5] # ee folder
    
# param1 = 'gonzalezivan'
# param2 = 'projects/gonzalezivan/assets/cola/anoa1'
# param3 = 'EPSG:3395'
# param4 = 10000
# param5 = 'projects/gonzalezivan/cola/'
    
    proj_epsg = param3 # EPSG:3395
    grid_scale = param4
    asset_dir = param5    

    try:
      ee.Authenticate()
      ee.Initialize(project = param1)
    except:
      print(' ERROR: no EE initialized')
      

    # CCDC run parameters
    start_year = 1999
    end_year = 2024
    start_doy = 1
    end_doy = 366
    max_cloud_cover_land = 50
    min_sun_elev = 40
    bands = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
    sensors_meta = 'landsat_45789'
    
    region_geom = ee.Feature(ee.FeatureCollection(param2).geometry().bounds())

    # a grid covering the region of interest geometry
    grid_fc = region_geom.coveringGrid(proj=proj_epsg, scale=grid_scale)


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
    
    try:
      print(task_config)
      task = ee.batch.Export.table.toAsset(**task_config)
      task.start()
      print(' ++ CCDC submitted: ', task)
    except Exception as e:
      # Handle errors (e.g., asset doesn't exist)
      print(' ERROR, task not initialized: ')
      print(e)
      return False



# End function

if __name__ == "__main__":
    main()


