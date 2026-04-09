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
from shapely.geometry import Point, box


# /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_uploadFeature.py gonzalezivan /home/shiny/test/anoa66sp.shp projects/gonzalezivan/assets/cola/anoa66sp
def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    project_name = sys.argv[1] # project name
    shapefile_path = sys.argv[2] # shapefile path
    ee_name = sys.argv[3] # ee asset
    

  # param2 = '/data/tempR/colaBWW2025081407481705/Anoa_present_ardianti.shp'
  
  # project_name = 'gonzalezivan'
  # shapefile_path = '/home/shiny/test/anoa66.shp'
  # ee_name = 'projects/gonzalezivan/assets/cola/anoa66'
  # project_name = 'gonzalezivan'
  # shapefile_path = '/home/shiny/test/anoa66.shp'
  # shapefile_path = 'C:/cola/Anoa/Anoa_present_ardianti.shp'
  # ee_name = 'projects/gonzalezivan/assets/cola/anoa1'

      
  
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
      print('\nEE initialized\n')
    except:
      print('\nERROR: no EE initialized\n')
      exit(1)
   
    # try:
    #     # Attempt to get asset information
    #     info = ee.data.getInfo(ee_name)
    #     return info is None
    # except Exception as e:
    #     # Handle errors (e.g., asset doesn't exist)
    #     print(' Asset exists, error:')
    #     print(e)
    #     return False
  
  
    #print('reading')
    # Read the shapefile into a GeoDataFrame
    gdf = gpd.read_file(shapefile_path)    
    minx, miny, maxx, maxy = gdf.total_bounds
    
    # gdf.assign(preabs=1)


    # Create bounding box polygon
    extent_poly = box(minx, miny, maxx, maxy)

    # Save as shapefile
    gpdbox = gpd.GeoDataFrame(geometry=[extent_poly], crs=gdf.crs)
    
    # convert it into geo-json 
    json_df = json.loads(gdf.explode(index_parts=False).to_json())
    json_box = json.loads(gpdbox.explode(index_parts=False).to_json())
    
        
    # create a gee object with geemap
    ee_object = geemap.geojson_to_ee(json_df)
    ee_box = geemap.geojson_to_ee(json_box)
    
    
    # upload this object to earthengine
    #asset = os.path.join(folder, asset_name)
                
    #create and launch the task
    task_config = {
        'collection': ee_object, 
        'description':'Uploading shapefile', #+ hapefile_path,
        'assetId': ee_name.replace('.shp', '')
    }
    
    box_tax_config = {
        'collection': ee_box, 
        'description':'Uploading box ',
        'assetId': ee_name.replace('.shp', '')+'_box'
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
    
      
    try:
      print(box_tax_config)
      taskbox = ee.batch.Export.table.toAsset(**box_tax_config)
      taskbox.start()
      print(' ++ Upload box submitted: ', box_tax_config)
      
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
