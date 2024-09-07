# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:48:54 2023
Script to create source points spatially 
distributed according to habitat suitability.
If the min or max values provided by the user
are > or < the min and max values present in the
raster, we assume the user wants points distributed
within this restricted range. We therefore mask out
values outside the provided range.
 
@author: pj276
"""
#%%
# IMPORTS
import sys
import rasterio as rio
import cola_functions as cf
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
import time

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()

    # INPUTS    
    # Path to suitability or resistance grid
    srg = sys.argv[1]

    # Output file path
    ofile = sys.argv[2]

    # Suitability grid min value
    sMin = sys.argv[3]

    # Suitability grid max value
    sMax = sys.argv[4]
    
    # Number of source points to create
    numSourcePoints = sys.argv[5]
    
    # Is provided grid suitability? If not, assume resistance. Should by 'Yes' or 'No'
    isSuitability = sys.argv[6]
    
    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[7] # Default None

    #%%
    # Convert min, max and nsources to float or integer
    try:
        if cf.is_float(sMin):
            sMin = float(sMin)
        elif sMin.isdigit():
            sMin = int(sMin)
        else:
            float(sMin)
            int(sMin)
    except ValueError:
        print('Min value must be either a float or integer')
    try:
        if cf.is_float(sMax):
            sMax = float(sMax)
        elif sMax.isdigit():
            sMax = int(sMax)
        else:
            float(sMax)
            int(sMax)
    except ValueError:
        print('Max value must be either a float or integer')   
    try:
        numSourcePoints = int(numSourcePoints)
    except ValueError:
        print('Num source points must be integer.')
    
    #%%
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array(upCRS, srg)

    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)

    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')
   
    #%%
    # Convert no data value to nan
    r[r==-9999] = np.nan
    
    # Convert raster values outside the provided range to nan
    r[r < sMin] = np.nan
    r[r > sMax] = np.nan
    
    # Rescale to 0-1
    rRsc = (r - sMin) / (sMax - sMin)
    
    # If resistance, invert
    if isSuitability == 'No':
        rRsc = 1-rRsc
    
    #%%
    # Subtract random array from rescaled suitability array
    rRsc = rRsc - np.random.rand(r.shape[0],r.shape[1])

    # Get row, column indices where values are > 0
    rcIndices = np.argwhere(rRsc > 0)
    
    # Throw an error if too many source points
    # have been requested (i.e. if there are not enough
    # habitat cells > 0 to satisfy the request)
    if numSourcePoints > rcIndices.shape[0]:
        raise Exception('There are not enough habitat cells to satisfy the requested number of source points. Try requesting fewer source points.')

    # Get random subset of indices to index into row, column object.
    # Sort from smallest to largest.
    sSet = np.random.choice(rcIndices.shape[0], numSourcePoints, replace=False)
    sSet.sort()

    # Select random row, column indices
    randIndices = rcIndices[sSet,:]
    
    # Convert row, column to coordinates
    with rio.open(srg) as src:
        xCs, yCs = rio.transform.xy(src.transform, randIndices[:,0], randIndices[:,1])
    
    # Convert to pandas df
    dfCoords = pd.DataFrame({'X': xCs, 'Y': yCs})
    
    # Write to shapefile or csv
    if Path(ofile).suffix == '.csv':
        dfCoords.to_csv(ofile, index=False)
    elif Path(ofile).suffix == '.xy':
        dfCoords.to_csv(ofile, index=False)
    elif Path(ofile).suffix == '.shp':
        # Convert to geodataframe
        gdf = gpd.GeoDataFrame(dfCoords, geometry=gpd.points_from_xy(dfCoords.X, dfCoords.Y), crs=src.crs)
        gdf.to_file(ofile)
    else:
        raise Exception('Output source point file should have a .csv, .xy, or .shp extension.')
    
    toc = time.perf_counter()
    print(f"Calculating source points took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
