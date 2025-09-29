# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 11:43:04 2025

@author: pj276
"""

#%%
# IMPORTS
import sys
import cola_functions as cf
import numpy as np
from scipy.ndimage import gaussian_filter
import time
from joblib import Parallel, delayed
import shutil, os

def main() -> None:
    #%%
   
    # Start timer
    tic = time.perf_counter()
    
    # INPUTS    
    # Path to resistance grid
    rg = sys.argv[1] 

    # Output file path
    ofile = sys.argv[2]

    # Zarr file name
    dazarr = sys.argv[3]
    
    # Radius for gaussian smoother (in number of cells)
    # The size of the kernel on each side is 2*radius + 1
    # E.g. a radius of 2 gives a 5x5 cell kernel
    gRad = sys.argv[4]

    # Amount to add to the least cost path (in cost distance units)
    # in order to generate a swath of low cost pixels,
    # termed the least cost corridor.
    # If 0, returns the least cost path.
    # If > 0, this amount is added to the least cost path
    # value so that all pixels with values <= to that value
    # will be returned. This results in a swath of pixels
    # instead of a single pixel wide path.
    corrTolerance = sys.argv[5]

    # Number of threads to use
    nThreads = sys.argv[6] # Default 1

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[7] # Default None
   
    # point pair output file name
    reOrderFile = sys.argv[8]

    # nodeids file name
    nodeidsFile = sys.argv[9]

    # Start corridor index
    # For now, these should be zero indexed python style
    # E.g. for a landscape with 10,000 corridors
    # a first batch of corridors could be 0-500
    # Python range is such that this would process
    # corridors 0-499. Then next batch would be 500-1000,
    # which would process corridors 500-999. The next
    # batch would be 1000-1500, and so on.
    sci = sys.argv[10] # Default is None
    
    # End corridor index
    eci = sys.argv[11] # Default is None
    
    # Convert start corridor to integer
    try:
        if sci != 'None':
            sci = int(sci)
    except ValueError:
        print('Start corridor index must be integer or "None"')
    
    # Convert end corridor to integer
    try:
        if eci != 'None':
            eci = int(eci)
    except ValueError:
        print('End corridor index must be integer or "None"')
    
    # Convert gaussian smoother radius to integer
    try:
        gRad = int(gRad)
    except ValueError:
        print('Gaussian smoother radius must be integer.')

    # Corridor diameter to integer
    try:
        corrTolerance = int(corrTolerance)
    except ValueError:
        print('Corridor radius must be integer.')

    # Convert number of threads to integer
    try:
        nThreads = int(nThreads)
    except ValueError:
        print('Number of processors must be integer.')
    
    #%%
    # Get reorder list
    reOrder = np.loadtxt(reOrderFile, delimiter=',').astype(np.int32)

    # Get nodeids
    nodeids = np.loadtxt(nodeidsFile, delimiter=',').astype(np.int32)
    
    # Get nodeids length
    nodeidsLen = len(nodeids)
        
    #%%
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array(upCRS, rg)
    
    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)
    
    # Check values between 0-1. If detected, this converts those values to
    # 1 and prints a message.
    r, profile = cf.checkRasterVals(r, profile)
    
    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')

    #%%
    # Use joblib to calculate corridors in parallel and sum on the fly
    if (sci != 'None') and (eci != 'None') and (eci > sci):
        startCorridor = sci
        endCorridor = eci
    else:
        startCorridor = 0
        endCorridor = len(reOrder)
        
    print("Summing corridors", flush=True)
    res = Parallel(n_jobs=nThreads, return_as="generator_unordered", prefer="threads") (
        delayed(cf.calcPathsMod)(reOrder[i], dazarr, corrTolerance) for i in range(startCorridor,endCorridor)
        )
    lccSum = np.zeros((1,nodeidsLen))
    corridorCount = 1
    tic1 = time.perf_counter()
    for result in res:
        lccSum += result
        if corridorCount % 10000 == 0:
            toc1 = time.perf_counter()
            print(f'Finished corridor {corridorCount}', flush=True)
            print(f"in {((toc1 - tic1)/60):0.4f} minutes", flush=True)
            
        corridorCount += 1
        #print(f"Finished corridor {corridorCount}", flush=True)
    
    # Write to file
    print('Writing corridors to file', flush=True)
    
    # Create 1D zeros array to hold lcc values
    dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    # Update with lcc values
    dArr[nodeids] = lccSum 
    # Reshape into rectangular array
    dArr = dArr.reshape(r.shape[0],r.shape[1])
    del lccSum
    del nodeids
    del r
    
    #%%
    # Use gaussian filter to smooth least cost values
    # Sigma is standard deviation of the kernel
    # The size of the kernel on each side is 2*radius + 1
    # E.g. a radius of 2 gives a 5x5 cell kernel
    # If radius is 0, don't smooth
    if gRad == 0:
        lccSmooth = dArr
    else:
        lccSmooth = gaussian_filter(dArr, sigma=1, radius=gRad)
    del dArr
    #%% 
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    lccSmooth = np.expand_dims(lccSmooth, axis=0)
    
    # Write to file
    cf.arrayToGeoTiff(lccSmooth, ofile, profile)
    
    # If zarr files were created, delete them
    if os.path.exists(dazarr):
        try:
            shutil.rmtree(dazarr)
            print(f"Successfully deleted Zarr store at: {dazarr}")
        except OSError as e:
            print(f"Error deleting Zarr store at {dazarr}: {e}")
    else:
        print(f"Zarr store not found at: {dazarr}")

    # Delete csv files
    if os.path.exists(reOrderFile):
        os.remove(reOrderFile)
    if os.path.exists(nodeidsFile):
        os.remove(nodeidsFile)

    toc = time.perf_counter()
    print(f"Calculating least cost corridors took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()

