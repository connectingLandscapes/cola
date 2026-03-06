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
import geemap
import pandas as pd 
import geopandas as gpd
import json
from sat_ts_fusion.lidar import gedi_utils

# /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_uploadFeature.py gonzalezivan /home/shiny/test/anoa66sp.shp projects/gonzalezivan/assets/cola/anoa66sp
def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    project_name = sys.argv[1] # project name // gonzalezivan ///
    shapefile_path = sys.argv[2] # shapefile path
    ee_name = sys.argv[3] # ee asset
    
    

  # param2 = '/data/tempR/colaBWW2025081407481705/Anoa_present_ardianti.shp'
  
  # project_name = 'gonzalezivan'
  # shapefile_path = '/home/shiny/test/anoa66.shp'
  # ee_name = 'projects/gonzalezivan/assets/cola/anoa66'
  # project_name = 'gonzalezivan'
    #shapefile_path = '/home/shiny/test/anoa66.shp'
    #shapefile_path = 'C:/cola/Anoa/Anoa_present_ardianti.shp'
    #ee_name = 'projects/gonzalezivan/assets/cola/anoa1'

      
  
  #param = 'gonzalezivan'
  # Convert distance threshold to float or integer
    try:
      ee1 = ee.Initialize(project = project_name)
    except:
      ee.Authenticate()
      ee.Initialize(project = project_name)
    
    
    
    try:
      ee.Authenticate()
      ee.Initialize(project = project_name)
    except:
      print(' ERROR: no EE initialized')
      
      
      
    gedi_keep_cols = ['shot_num', 'beam', 'delta_time', 'date_dec', 'date', 'time', 'doy', 'sol_elev', 'sol_azim',
                  'lon_lm_a0', 'lat_lm_a0', 'elev_lm_a0', 'elev_dem_srtm', 'elev_dem_tdx', 'elev_diff_dem', 'deg_flag',
                  'l2a_qflag_a0', 'l2_hqflag', 'sens_a0', 'ls_treecov', 'ls_waterp', 'leafoff_flag', 'leafoff_doy',
                  'leafon_doy', 'pft', 'urb_prop',
                  'num_modes_a0', 'rh_25_a0', 'rh_50_a0', 'rh_75_a0', 'rh_90_a0', 'rh_95_a0', 'rh_98_a0', 'rh_99_a0',
                  'rh_100_a0',
                  'rhvdr_b', 'rhvdr_m', 'rhvdr_t', 'l2b_qflag_a0', 'fhd_pai_1m_a0', 'cover_a0', 'pai_a0', 'pai_l1',
                  'pai_l2', 'pai_l3', 'pai_l4', 'pai_l5', 'pai_l6', 'pai_l7', 'pai_l8', 'pai_l9', 'pai_l10', 'pai_l11',
                  'pai_l12', 'pai_l13', 'pai_l14', 'pai_l15', 'pai_l16', 'pai_l17', 'pai_l18', 'pai_l19', 'pai_l20',
                  'pai_l21',
                  'pavd_0_5', 'pavd_5_10', 'pavd_10_15', 'pavd_15_20', 'pavd_20_25', 'pavd_25_30', 'pavd_30_35',
                  'pavd_35_40', 'pavd_40_45', 'pavd_45_50', 'pavd_50_55', 'pavd_55_60', 'pavd_60_65', 'pavd_65_70',
                  'pavd_70_75', 'pavd_75_80', 'pavd_80_85', 'pavd_85_90', 'pavd_90_95', 'pavd_95_100',
                  'pavd_0_5_frac', 'pavd_bot_frac', 'pavd_top_frac', 'pavd_max_h', 'even_pai_1m_a0', 'fhd_pavd_5m_a0',
                  'even_pavd_5m_a0', 'l4a_qflag_a0', 'l4a_hqflag', 'agbd_a0', 'agbd_se_a0'] 
    # where to download and extract GEDI tables
    local_working_dir = 'D:/PB_working/reshape/gedi/'

    # Google Cloud Storage bucket name
    gcs_bucket_name = 'gedi_uploads'
    # Google Cloud Storage subdirectory. Leave as "" to store in root
    gcs_dir = ""
    # Google Cloud service account key associated with your project
    json_service_account_key = 'C:/Users/pb463/Downloads/ee-gedibio-c9474d011c82.json'  
    # TODO: how/where to store the Google Service Account key


    for subregion in ['mc', 'mwcf', 'scc', 'sr', 'wdm', 'wum', 'wcdm']:
        print(f'Working on subregion {subregion}...')
        # a shapefile for the area of interest
        aoi_shape_path = 'D:/PB_working/reshape/aoi/ReSHAPE_study_extent_'+subregion+'.shp'

        # an asset name for the FeatureCollection being uploaded
        fc_asset_name = 'gediv002_l2l4a_20190417to20230316_'+'reshape_'+subregion

        gedi_utils.make_gedi_l2l4a_fc(start_year=2019,
                                      end_year=2022,
                                      start_doy=1,
                                      end_doy=366,
                                      max_elev_ls=[1500,1500,1500,1500,2000,2500,3500,3500,3500,2500,1500,1500],
                                      aoi_shape=aoi_shape_path,
                                      max_shots_per_tile=100000,
                                      gedi_keep_cols=gedi_keep_cols,
                                      grid_crs='EPSG:6350',
                                      grid_width=10000,
                                      local_working_dir=local_working_dir,
                                      gcs_bucket_name=gcs_bucket_name,
                                      gcs_dir=gcs_dir,
                                      service_account_key=json_service_account_key,
                                      gee_project_name=gee_project_name,
                                      gee_fc_asset_dir='projects/'+gee_project_name+'/assets/',
                                      gee_fc_asset_name=fc_asset_name)
        print()







    ########################## OLD
    ########################## OLD
    
    



    # try:
    #     # Attempt to get asset information
    #     info = ee.data.getInfo(ee_name)
    #     return info is None
    # except Exception as e:
    #     # Handle errors (e.g., asset doesn't exist)
    #     print(' Asset exists, error:')
    #     print(e)
    #     return False
  
  
    print('reading')
    # Read the shapefile into a GeoDataFrame
    gdf = gpd.read_file(shapefile_path)
    # convert it into geo-json 
    json_df = json.loads(gdf.explode(index_parts=False).to_json())
    
        
    # create a gee object with geemap
    ee_object = geemap.geojson_to_ee(json_df)
        
    # upload this object to earthengine
    #asset = os.path.join(folder, asset_name)
                
    #create and launch the task
    task_config = {
        'collection': ee_object, 
        'description':'test',
        'assetId': ee_name.replace('.shp', '')
    }
  
    try:
      print(task_config)
      task = ee.batch.Export.table.toAsset(**task_config)
      task.start()
      print(' ++ Upload submitted: ', task)
    except Exception as e:
      # Handle errors (e.g., asset doesn't exist)
      print(' ERROR, task not initialized: ')
      print(e)
      return False
    
        # params2  = {
        #     'id': ee_name.replace('.shp', 'Ingest'),
        #     'description':'testfromTableIng',
        #     'assetId': ee_name.replace('.shp', '')
        # }
        # 
        # idx = ee.data.newTaskId()
        # ee.data.startTableIngestion(request_id = idx, params2, allow_overwrite=False)
    
    
        # manifest = {
        #   'name': 'projects/my-project/assets/asset-name',
        #   'tilesets': [
        #     {
        #       'sources': [
        #         {
        #           'uris': [
        #             'gs://my-bucket/filename.tif'
        #           ]
        #         }
        #       ]
        #     }
        #   ]
        # }
        # ee.data.startIngestion(None, manifest)
  
  
# End function

if __name__ == "__main__":
    main()
