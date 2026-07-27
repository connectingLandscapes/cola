# -*- coding: utf-8 -*-

"""
SDM GEE export to drive
==========================================

Export assets available in a EE folder

Usage:
  python gee_to_gd.py --ee_project project --ee_folder projects/PROJECT/assets/FOLDER \
      --gd_folder gee/gd_export --prefix cola2 --crs EPSG:4326 --scale 250 \
          --gee_assets projects/PROJ/assets --region projects/PROJ/assets/AOI
"""

import ee, argparse

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description='SDM MODIS Wall-to-Wall Prediction')
    # ── Core identifiers (shared across scripts) ──────────────────────────────
    p.add_argument('--ee_project',   type=str, default='geersprocessing')
    p.add_argument('--ee_folder',      type=str, required=True)
    p.add_argument('--gd_folder',      type=str, required=True,
                   help='Folder in your Google Drive account to save the results')
    p.add_argument('--prefix', type=str, default='cola2_',
                   help='Prefix on the files in your folder')
    # ── Spatial / scale (shared across scripts) ───────────────────────────────
    p.add_argument('--crs',          type=str, default='EPSG:4326')
    p.add_argument('--scale',        type=int, default=250)
    # ── Asset paths ───────────────────────────────────────────────────────────
    p.add_argument('--gee_assets',   type=str, required=True,
                   default='projects/geersprocessing/assets')
    p.add_argument('--region',  type=str, required=True,
                   help='Feature with the extent to export the layers')
    return p.parse_args()
## Output to be saved here:  desc = safe_name(f'{MODEL_ID}_suitability_{TARGET_YEAR}_{tile["name"]}')

args = parse_args()
ee.Initialize(project=args.ee_project)
# ee.Initialize(project='gonzalezivan')
# EEFOLDER = 'projects/gonzalezivan/assets/cola/m6_prediction_2025'
# REGION = 'projects/gonzalezivan/assets/cola/sulawesi'

# =============================================================================
# 1. DERIVED PARAMETERS
# =============================================================================

EEFOLDER = args.ee_project
EEFOLDER = args.ee_folder
GDFOLDER = args.gd_folder
PREFIX = args.prefix
SCALE = args.scale
CRS = args.crs
REGION = args.region

# =============================================================================
# GDFOLDER = 'gee'
# PREFIX = 'pref'
# SCALE = 250
# CRS = 'EPSG:4326'
# =============================================================================



def main() -> None:
    #%%
    reg = ee.FeatureCollection( REGION )
    region = reg.bounds() # Define the export region
    result = ee.data.listAssets({"parent": EEFOLDER})
    results = result.get('assets')
    for rr in results:
        DESC = PREFIX +'_Cola2_'+ rr.get('name').replace(EEFOLDER, '').replace('/', '')
        nname = rr.get('name') 
        try:
            task = ee.batch.Export.image.toDrive(
                image= ee.Image(rr.get('name')),
                description = DESC,
                folder = GDFOLDER,
                fileNamePrefix = DESC,
                region = region,
                scale = SCALE,
                crs = CRS,
                maxPixels=1e9)
            print(f'Task submited: {nname}')
            task.start()
        except Exception as e:
            print(f' ERROR in {nname}: {e}')
    print( '\n'  * 2)
    print( f'Finished. {len(results)} layers export tasks submmited.')
    print( '\n'  * 2)
#
#    
if __name__ == "__main__":
    main()
