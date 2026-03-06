# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Cinnects EE from command line
@author: ig299
"""

#%%
# IMPORTS
import ee
import sys
import csv
from sat_ts_fusion.imagery import lsat_utils


# /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_lsat_ccdc.py /home/shiny/sat_ts_fusion/params_ccdc.csv

def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    paramcsv = sys.argv[1] # project name
    
    with open(paramcsv) as f:
      dic = dict(filter(None, csv.reader(f)))
    
    # print(dic)
    
    # param1 = 'gonzalezivan'
    # param2 = 'projects/gonzalezivan/assets/cola/anoamanual'
    # param3 = 10000
    # param4 = 'projects/gonzalezivan/cola/sim1'
  
    # param2 = sys.argv[2] # ee point layer
    # param3 = sys.argv[3] # grid scale
    # param4 = sys.argv[4] # ee folder
    # param5 = sys.argv[5] # EPSG:3395
    
    aoi = dic['aoi']
    eeproject = dic['eeProject']
    grid_scale = int(dic['gridScale'])
    asset_dir = dic['assetDir']
    proj_epsg = dic['projEpsg']
    
    start_year = int(dic['startYear'])
    end_year = int(dic['endYear'])
    start_doy = int(dic['startDoy'])
    end_doy = int(dic['endDoy'])
    max_cloud_cover_land = int(dic['maxCloudCover_land'])
    min_sun_elev = int(dic['minSunElev'])
    bands = dic['bands'].strip().split(',')
    sensors_meta = dic['sensorsMeta']
    break_point_bands = dic['breakPointBands'].strip().split(',')
    tmask_bands = dic['tmaskBands'].strip().split(',')
    min_observations = int(dic['minObservations'])
    chi_sq_prb = float(dic['chiSquareProbability'])
    min_numb_years_scale = float(dic['minNumOfYearsScaler'])
    lamb = float(dic['lambda0'])
    max_iter = int(dic['maxIterations'])
    date_format = int(dic['dateFormat'])
    export_scale = int(dic['exportScale'])
    out_crs = dic['outCrs']
    
    # # CCDC run parameters
    # start_year = 1999
    # end_year = 2024
    # start_doy = 1
    # end_doy = 366
    # max_cloud_cover_land = 50
    # min_sun_elev = 40
    # bands = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
    # sensors_meta = 'landsat_45789'
    
    try:
        ee.Authenticate()
        ee.Initialize(project = eeproject)
    except:
        print(' ERROR: no EE initialized')
    
    region_geom = ee.Feature(ee.FeatureCollection(aoi).geometry().bounds()).geometry()
    
    # a grid covering the region of interest geometry
    grid_fc = region_geom.coveringGrid(proj=proj_epsg, scale=grid_scale)
    
    ee.data.createFolder(asset_dir)
    #ee.data.createAsset({'type': 'Folder'}, asset_dir)

    # Landsat Collection 2 surface reflectance image collections
    l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
    l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
    l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
    
    # Merge Landsat collections
    l4to9 = l5.merge(l7).merge(l8).merge(l9)
    
    ccdcParams = {
      'breakpointBands':break_point_bands,
      'tmaskBands' : tmask_bands,
      'minObservations': min_observations,
      'chiSquareProbability': chi_sq_prb,
      'minNumOfYearsScaler': min_numb_years_scale,
      'lambda': lamb,
      'maxIterations' : max_iter,
      'dateFormat' : date_format
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
    
    n_grids=ee.Number(grid_fc.size()).getInfo()
    print ( 'Launching ' + str(n_grids) + ' tasks')
    
    for i in range(0, n_grids): # n_grids
    #if False:
        g_geom = ee.Feature(grid_fc.toList(n_grids).get(i)).geometry()
        
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
        #print(l4to9_p_g.first().getInfo())
        
        ccdcParams['collection'] = l4to9_p_g
        
        ccdc_array = (ee.Algorithms.TemporalSegmentation.Ccdc(**ccdcParams)
                     .clip(g_geom)
                     .set('grid_i', i)
                     .set(runParams)
                     .set(ccdcParams))
        
        task_name = str(asset_dir)+'__'+str(i) # 'lsatc2sr_ccdc_grid__'+
        assetID = asset_dir +'/'+str(i)
        
        task = ee.batch.Export.image.toAsset(
          image=ccdc_array, 
          description=task_name.replace('/', '_'), 
          assetId=assetID,
          pyramidingPolicy={".default": "sample"}, 
          scale=export_scale, crs=out_crs, maxPixels=1e13)
        
        #assetId=f"{asset_dir}{task_name}",
        
        try:
          task.start()
          print(' ++ CCDC submitted: ', task)
        except Exception as e:
          # Handle errors (e.g., asset doesn't exist)
          print(' ERROR, task not initialized: ')
          print(e)



# End function

if __name__ == "__main__":
    main()
