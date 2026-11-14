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
#from geemap import geojson_to_ee, ee_to_geojson
import pandas as pd 
import geopandas as gpd
import json
from shapely.geometry import Point, box
from sat_ts_fusion import eecolatools as ect
import os
from pathlib import Path
import argparse

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='Download GEE files to Google Drive')
    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--point_path',      type=str, required=True,
                   help = 'Path to the local point file. If there spaces "quote the whole path"')
    p.add_argument('--ee_layer',      type=str, required=True,
                   help='Path of the layer to create')
    p.add_argument('--presence_col', type=str, default='preabs',
                   help='Name of the column containing the 0/1 integer values of presence and absences'
                   'If not provided, preabs col will be created and all points filled with 1')
    p.add_argument('--number_abs', type=str, default='0',
                   help='Number of absences to simulate. Can be an integer or a % of the original presence'
                   ' e.g. 100 for ine hundred points or 100% for same amount in the points in the layer')
    p.add_argument('--abs_method',  type=str, required=False,
                   help='Method for creating absences: background, blocks, environmental'
                   'similarity')
    p.add_argument('--block_size',  type=str, required=False,
                   help='Block size to create absences. Default 10km if not provided')

    return p.parse_args()

args = parse_args()

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

# Path to file holding xy coordinates
project_name = args.ee_project # project name
shapefile_path = args.point_path # shapefile path
ee_layer = args.ee_layer # ee asset
number_abs = str(args.number_abs) # % of absences based of presences
presence_col = args.presence_col
abs_method = args.abs_method
block_size = args.block_size

# percabs 


# import os
# os.chdir('N:/My Drive/git/cola/inst/ee')
# /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_uploadFeature.py gonzalezivan /home/shiny/test/anoa66sp.shp projects/gonzalezivan/assets/cola/anoa66sp
# C:/Users/ig299/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_uploadFeature.py" gonzalezivan C:/Users/ig299/cola/ptsa.shp projects/gonzalezivan/assets/cola/aaupl20 20
# C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_uploadFeature.py" gonzalezivan C:/cola/anoa/anoa_test.shp projects/gonzalezivan/assets/cola/anoa66 20
def main() -> None:
    #%%

    # INPUTS
    # Path to file holding xy coordinates
    #project_name = sys.argv[1] # project name
    #shapefile_path = sys.argv[2] # shapefile path
    #ee_layer = sys.argv[3] # ee asset
    #percabs = int(sys.argv[4]) # % of absences based of presences
    try:
        eelogpath = str(sys.argv[5]) # 
        if not os.path.isdir(eelogpath):
            Path(eelogpath).mkdir(parents=True, exist_ok=True)
        if not os.path.isdir(eelogpath):
            eelogpath = ""
    except IndexError:
        eelogpath = ""
        
        
    # shapefile_path = '/home/shiny/test/anoa66.shp'
    # shapefile_path = 'C:/Users/ig299/cola/ptsa.shp'
    # shapefile_path = 'C:/cola/tiger40pts.shp'
    # eelogpath = 'C:/Users/ig299/cola/eelogpath'
    # project_name = 'gonzalezivan'
    # percabs = 0 60
    # ee_layer = 'projects/gonzalezivan/assets/cola2/tiger50'  
  
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
    #     info = ee.data.getInfo(ee_layer)
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
    # Estimating how many absences
    if '%' in number_abs:
        try:
            numabs = float(number_abs.replace("%", ""))
            numabs  = int((float(numabs)/100)*gdf.size)
        except:
            print( 'No valid percentage. Skipping absences')
    else:
        try:
            numabs = float(number_abs)
        except:
            print( 'No valid number. Skipping absences')
        
        
    # Simulate abs
    if numabs != '0':
        # bounds = ee.Geometry.Rectangle(extents.values(["longitude_min", "latitude_min", "longitude_max", "latitude_max"]))
        # bounds = ee.Geometry.Rectangle([minx, miny, maxx, maxy])
        # refimage = ee.Image("USGS/SRTMGL1_003").clip(bounds)
        # extabs = refimage.sample(scale=30, numPixels=percabs, geometries=True)
        from shapely.geometry import Polygon
        # x = 3, 4; y = -1, 1
        # s = gpd.GeoSeries([Polygon([(3, -1), (4, -1), (4, 1), (3, 1)])])
        
        if abs_method == 'background':
            0
        else if abs_method == 'blocks':
            1
        else if abs_method == 'environmental':
            2
        #
        pol = gpd.GeoSeries([Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])])
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
    # Convert it into geo-json 
    json_df = json.loads(gdf.explode(index_parts=False).to_json())
    json_box = json.loads(gpdbox.explode(index_parts=False).to_json())
    #
    #    
    # create a gee object with geemap
    # ee_object = geemap.geojson_to_ee(json_df) # old - no working
    ee_object = ee.FeatureCollection(json_df)
    #ee_object = geemap.shp_to_ee(shapefile_path)
    #ee_object = geemap.geopandas_to_ee(json_df)
    #ee_object = geemap.geopandas_to_ee(gdf)
    #ee_box = geemap.geojson_to_ee(json_box) # old - no working
    ee_box = ee.FeatureCollection(json_box)
    #
    #
    # upload this object to earthengine
    #asset = os.path.join(folder, asset_name)
    #
    #%%
    #create and launch the task
    task_config = {
        'collection': ee_object, 
        'description':'Uploading shapefile', # + ee_layer, #+ hapefile_path,
        'assetId': ee_layer.replace('.shp', '')
    }
    #
    box_tax_config = {
        'collection': ee_box, 
        #'description':'Uploading box '+ee_layer,
        'description':'Uploading box',
        'assetId': ee_layer.replace('.shp', '')+'_box'
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
      print(' ERROR, box task not initialized: ')
      print(e)
      return False
    
    
    # params2  = {
        #     'id': ee_layer.replace('.shp', 'Ingest'),
        #     'description':'testfromTableIng',
        #     'assetId': ee_layer.replace('.shp', '')
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
