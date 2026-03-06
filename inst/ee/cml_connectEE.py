# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Connects EE from command line
@author: ig299
"""

#%%
# IMPORTS
import ee
import sys 

def main() -> None:
    #%%
    # Start timer

    # INPUTS
    # Path to file holding xy coordinates
    param1 = sys.argv[1] # project name # param1 = 'gonzalezivan'
    param2 = sys.argv[2] # tempfolder were to write # param2 = 'C:/cola/colaNAR2026030412391405/'
    
    # C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe C:/cola/cola2/cml_connectEE.py gonzalezivan C:/cola/colaNAR2026030412391405/


    # param1 = 'gonzalezivan'
    # param2 = 'C:/cola/cola2/id2000.txt'
    # Convert distance threshold to float or integer
    try:
        ee.Authenticate()
        ee.Initialize(project = param1)
        print('\n Cola2: EE initialized')
    except:
        print('\n Cola2 ERROR: no EE initialized\n')
        exit(1)
        

    # Convert gaussian smoother radius to integer
    try:
        collection = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS')
        collection_id = collection.getMapId().get('mapid')
        print('\n Cola2: EE conection works!\n')
        with open(param2+'/ee_conected.txt', "w") as file:
          file.write(collection_id )
          # print(collection.getMapId())
        #
    except ValueError as e:
        print('Cola 2 ERROR: ', e)
        exit(1)


if __name__ == "__main__":
    main()
