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
from sat_ts_fusion import eecolatools as ect
import os

# import os
# os.chdir('N:/My Drive/git/cola/inst/ee')
# /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_uploadFeature.py gonzalezivan /home/shiny/test/anoa66sp.shp projects/gonzalezivan/assets/cola/anoa66sp
# C:/Users/ig299/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_uploadFeature.py" gonzalezivan C:/Users/ig299/cola/ptsa.shp projects/gonzalezivan/assets/cola/aaupl20 20
# C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_uploadFeature.py" gonzalezivan C:/cola/anoa/anoa_test.shp projects/gonzalezivan/assets/cola/anoa66 20
def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    project_name = sys.argv[1] # project name
    shapefile_path = sys.argv[2] # shapefile path
    ee_name = sys.argv[3] # ee asset
    percabs = int(sys.argv[4]) # % of absences based of presences
    try:
        eelogpath = str(sys.argv[5]) # 
        if not os.path.isdir(eelogpath):
            Path(eelogpath).mkdir(parents=True, exist_ok=True)
        if not os.path.isdir(eelogpath):
            eelogpath = ""
    except IndexError:
        eelogpath = ""
        
        
    # shapefile_path = '/home/shiny/test/anoa66.shp'
    # shapefile_path = 'C:/cola/anoa/anoa66.shp'
    # shapefile_path = 'C:/Users/ig299/cola/ptsa.shp'
    # shapefile_path = 'C:/cola/anoa/anoa_test.shp'
    # eelogpath = 'C:/Users/ig299/cola/eelogpath'
    # project_name = 'gonzalezivan'
    # percabs = 60
    # ee_name = 'projects/gonzalezivan/assets/cola/a'  
  
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
  
    #%%
    #print('reading')
    # Read the shapefile into a GeoDataFrame
    gdf = gpd.read_file(shapefile_path)
    if 'MultiPoint' in list(set(gdf.geom_type)):
        print (' Exploding multipoint')
        # Convert MULTIPOINT to individual POINT features
        from shapely.geometry import MultiPoint
        #
        def explode_multipoints(geom):
            if geom.geom_type == 'MultiPoint':
                return list(geom.geoms)
            return [geom]
        #
    #
    gdf = gdf.explode(index_parts=False)
    #
    minx, miny, maxx, maxy = gdf.total_bounds
    #
    # Assign preabs column to 0 if not existing
    if 'preabs' not in gdf.columns:
        gdf = gdf.assign(preabs=1)
        print (  ' Creating the missing "preabs" column, assigning 1 to all features')
    # Create bounding box polygon
    extent_poly = box(minx, miny, maxx, maxy)
    gpdbox = gpd.GeoDataFrame(geometry=[extent_poly], crs=gdf.crs)
    #
    #
    if percabs != 0:
        # bounds = ee.Geometry.Rectangle(extents.values(["longitude_min", "latitude_min", "longitude_max", "latitude_max"]))
        # bounds = ee.Geometry.Rectangle([minx, miny, maxx, maxy])
        # refimage = ee.Image("USGS/SRTMGL1_003").clip(bounds)
        # extabs = refimage.sample(scale=30, numPixels=percabs, geometries=True)
        from shapely.geometry import Polygon
        # x = 3, 4; y = -1, 1
        # s = gpd.GeoSeries([Polygon([(3, -1), (4, -1), (4, 1), (3, 1)])])
        pol = gpd.GeoSeries([Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])])
        numabs  = int((percabs/100)*gdf.size)
        extabs = pol.sample_points(size = numabs )
        abspts = gpd.GeoDataFrame(geometry=extabs, crs=gdf.crs).explode(index_parts=True)
        abspts = abspts.assign(preabs=0)
        # abspts.plot()
        print (  ' Creating ' + str(numabs ) +' absences, assigning 0 to all new features in the "preabs" column')
        gdf = pd.concat([gdf, abspts], ignore_index=True)
        # gdf2.plot(c=gdf2['preabs']+1) 
        gdf.to_file(shapefile_path.replace(".shp", "and" + str(percabs) + "simabs.shp") )
        print (  ' Writing file ')
    #    
    #
    #
    # convert it into geo-json 
    json_df = json.loads(gdf.explode(index_parts=False).to_json())
    json_box = json.loads(gpdbox.explode(index_parts=False).to_json())
    #        
    # create a gee object with geemap
    ee_object = geemap.geojson_to_ee(json_df)
    ee_box = geemap.geojson_to_ee(json_box)
    #
    #
    # upload this object to earthengine
    #asset = os.path.join(folder, asset_name)
    #
    #%%
    #create and launch the task
    task_config = {
        'collection': ee_object, 
        'description':'Uploading shapefile', #+ hapefile_path,
        'assetId': ee_name.replace('.shp', '')
    }
    #
    box_tax_config = {
        'collection': ee_box, 
        'description':'Uploading box ',
        'assetId': ee_name.replace('.shp', '')+'_box'
    }
    #
    try:
      print(task_config)
      task = ee.batch.Export.table.toAsset(**task_config)
      task.start()
      print(' ++ Upload submitted: ', task)
      if not eelogpath == "":
          ect.writeTask(task, eelogpath)
    except Exception as e:
      # Handle errors (e.g., asset doesn't exist)
      print(' ERROR, task not initialized: ')
      print(e)
      return False
    
      
    try:
      print(box_tax_config)
      taskbox = ee.batch.Export.table.toAsset(**box_tax_config)
      taskbox.start()
      print(' ++ Upload box submitted: ', taskbox)
      if not eelogpath == "":
          ect.writeTask(taskbox, eelogpath)
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
