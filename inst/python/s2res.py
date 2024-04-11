# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 14:00:33 2023
Script to rescale habitat suitability values
to resistance values.
@author: pj276
"""

#%%
# Imports
import sys
import osgeo
import cola_functions as cf
import numpy as np
import rasterio as rio
from rasterio.crs import CRS

def main() -> None:
    #%%
    # INPUTS    
    # Path to suitability grid
    sg = sys.argv[1] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example_pro.tif'

    # Output file path for resistance grid
    ofile = sys.argv[2] #r"C:/Users/pj276/Projects/UNICOR_arch/unicor/sabah_example_XY200_nik_crk.tif"

    # Suitability minimum value
    sMin = sys.argv[3]
    
    # Suitability maximum value
    sMax = sys.argv[4]

    # Maximum resistance value (min is set to 1)
    rMax = sys.argv[5]
    
    # Shape parameter
    shpParam = sys.argv[6]
    
    # No data value of suitability raster
    ndVal = sys.argv[7]
    
    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[8] # Default None
    
    # Convert suitability min to float or integer
    try:
        if cf.is_float(sMin):
            sMin = float(sMin)
        elif sMin.isdigit():
            sMin = int(sMin)
        else:
            float(sMin)
            int(sMin)
    except ValueError:
        print('Suitability min must be either a float or integer')

    # Convert suitability max to float or integer
    try:
        if cf.is_float(sMax):
            sMax = float(sMax)
        elif sMax.isdigit():
            sMax = int(sMax)
        else:
            float(sMax)
            int(sMax)
    except ValueError:
        print('Suitability max must be either a float or integer')

    # Convert resistance max to float or integer
    try:
        if cf.is_float(rMax):
            rMax = float(rMax)
        elif rMax.isdigit():
            rMax = int(rMax)
        else:
            float(rMax)
            int(rMax)
    except ValueError:
        print('Resistance max must be either a float or integer')

    # Convert shape param to float or integer
    try:
        if cf.is_float(shpParam):
            shpParam = float(shpParam)
        elif shpParam.lstrip('-').isdigit():
            shpParam = int(shpParam)
        else:
            float(shpParam)
            int(shpParam)
    except ValueError:
        print('Shape param must be either a float or integer')

    # Convert no data value to float or integer
    try:
        if cf.is_float(ndVal):
            ndVal = float(ndVal)
        elif ndVal.lstrip('-').isdigit():
            ndVal = int(ndVal)
        elif ndVal == 'nan':
            ndVal = -9999
        else:
            float(ndVal)
            int(ndVal)
    except ValueError:
        print('No data value must be either float, integer, or nan string')

    #%%
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array(upCRS, sg)
  
    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)

    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')
        
    # Rescale to 0-1 range, then to 1 - max resistance range
    resVec = cf.vecrescale(r, sMin, sMax, rMax, ndVal, -9999, True, shpParam)

    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    resVec = np.expand_dims(resVec, axis=0)
    
    # Write to file
    cf.arrayToGeoTiff(resVec, ofile, profile)

if __name__ == "__main__":
    main()