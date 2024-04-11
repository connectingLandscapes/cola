# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 20:39:07 2023
Script to calculate cumulative resistant kernels using Networkit.
Runs all source points at the same time with no initial distance
threshold. This can produce large arrays so the code
determines whether the product of n points and n cells
will consume more than 10GB of RAM. If so, it breaks the
calculation into batches.
@author: pj276
"""

#%%
# IMPORTS
import sys
import osgeo
import cola_functions as cf
import networkit as nk
import rasterio as rio
from rasterio.crs import CRS
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
    # Path to file holding xy coordinates
    xyf = sys.argv[1] # 'C:/Users/pj276/Projects/CDPOP_arch/data/sabah_example1_1686679883/batchrun0mcrun0/XY200.csv'
    
    # Path to resistance grid
    rg = sys.argv[2] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example_pro.tif'

    # Output file path
    ofile = sys.argv[3] #r"C:/Users/pj276/Projects/UNICOR_arch/unicor/sabah_example_XY200_nik_crk.tif"

    # Distance threshold (in cost distance units)
    dThreshold = sys.argv[4] # 500000
    
    # Transform. Can be either linear, gaussian, inverse, or inversesquare
    tForm = sys.argv[5] # linear
    
    # Kernel volume
    kvol = sys.argv[6]
    
    # Number of threads to use. 
    nThreads = sys.argv[7] # Default 1.

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[8] # Default None

    # Memory to use when reading in sources
    gbLim = sys.argv[9]

    # Convert kernel volume to float or integer
    try:
        if cf.is_float(kvol):
            kvol = float(kvol)
        elif dThreshold.isdigit():
            kvol = int(kvol)
        else:
            float(kvol)
            int(kvol)
    except ValueError:
        print('Kernel volume value must be either a float or integer')
    
    # Convert kernel volume to float or integer
    try:
        if cf.is_float(dThreshold):
            dThreshold = float(dThreshold)
        elif dThreshold.isdigit():
            dThreshold = int(dThreshold)
        else:
            float(dThreshold)
            int(dThreshold)
    except ValueError:
        print('Threshold value must be either a float or integer')

    # Convert number of threads to integer
    try:
        nThreads = int(nThreads)
    except ValueError:
        print('Number of processors must be integer.')

    # Convert memory limit to integer
    try:
        gbLim = int(gbLim)
    except ValueError:
        print('Memory limit must be integer')
    
    #%%
    # Set number of processors
    nk.setNumberOfThreads(nThreads)
    
    #%%
    # Read xy file to dataframe
    if Path(xyf).suffix == '.csv':
        xy = pd.read_csv(xyf)
    elif Path(xyf).suffix == '.shp':
        xy = gpd.read_file(xyf)
        xy = xy.get_coordinates()
    else:
        raise Exception('Source point file should have a .csv or .shp extension.')
    
    #%%
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array(upCRS, rg)

    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)

    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')

    # Get cell size from profile
    cSize = profile['transform'][0]
    
    #%%
    # Convert resistance grid to graph
    print("Converting image to graph", flush=True)
    nkG, nodeids, idmap = cf.image_to_graph(r, cSize, -9999, 8)
    print(nk.overview(nkG))
    
    #%%
    # Get sources as networkit node ids
    with rio.open(rg) as src:
        # Convert xy coordinates to rows and columns using
        # the resistance grid as a reference.
        cinds = cf.cell_indices_from_coords(src, r, np.array(xy))
        # Raise an exception if there are no source points in the landscape
        if np.sum(np.isnan(cinds)) >= 1:
            raise Exception('Source points do not intersect resistance grid.')
        else:
            # Get sources as networkit nodeids
            sources = cf.inds2nodeids(cinds, r, idmap)
            del idmap

    #%%
    # Estimate memory required for corridor extraction by multiplying
    # the number of graph nodes by the number of sources
    gbCf = 0.000000000125 # bits to gb conversion factor
    memReq = len(nodeids)*len(sources)*64*gbCf

    # If it's more than the threshold amount, divide into batches
    gbCLimit = gbLim
    if memReq > gbCLimit:

        # Number of corridor batches
        nCBatches = int(np.ceil(memReq/gbCLimit))
        print('Memory required for kernels exceeds ' + str(gbCLimit) + 'GB')
        print('Calculating kernels in ' + str(nCBatches) + ' batches')
        # Divide sources into batches
        sBatches = np.array_split(sources, nCBatches)            
        # Create zeros array to hold corridor sums
        crkSum = np.zeros((len(nodeids),))
        
        # Loop over batches
        for s in sBatches:
            # Calculate shortest path distance from each source to every
            # cell in the landscape
            spspDist = nk.distance.SPSP(nkG, s)
            spspDist.run()
            ccArr = spspDist.getDistances(asarray=True)
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
            ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        
            # Transform distances
            if tForm == 'linear':
                # Linear transform of distances
                ccArr = 1 - (1/dThreshold) * ccArr        
                # Convert negative values  to 0. (These are beyond the threshold)
                ccArr[ccArr < 0] = 0
                # Set nan to 0
                ccArr[np.isnan(ccArr)] = 0
            elif tForm == 'inverse':
                # Inverse square transform of distances
                ccArr = 1/(ccArr + 1)
                # Convert values beyond the distance threshold to 0
                ccArr[ccArr < 1/(dThreshold + 1)] = 0
                # Set nan to 0
                ccArr[np.isnan(ccArr)] = 0
            elif tForm == 'gaussian':
                dispScale = dThreshold/4
                ccArr = np.exp(-1*((ccArr**2)/(2*(dispScale**2))))
                # Set values beyond the distance threshold to 0
                ccArr[ccArr < np.exp(-1*((dThreshold**2)/(2*(dispScale**2))))] = 0
                # Set nan to 0
                ccArr[np.isnan(ccArr)] = 0
            
            # Multiplier
            if kvol > 1:
                vc = kvol * 3/(np.pi*dThreshold**2)
                ccArr = ccArr*vc
            elif kvol < 1:
                raise Exception('Kernel volume should be >= 1.')
            
            # Sum batch of kernals and 
            # add to accumulating sum
            crkSum += np.sum(ccArr, axis=0)
            del ccArr

    # Otherwise, process as normal
    else:
        # Calculate shortest path distance from each source to every
        # cell in the landscape
        print('Calculating kernels in 1 batch', flush=True)
        spspDist = nk.distance.SPSP(nkG, sources)
        spspDist.run()
        ccArr = spspDist.getDistances(asarray=True)
        # Networkit gives inaccessible nodes the max float 64 value.
        # Set these to nan.
        ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        ccArr = ccArr.astype('float32')
        del nkG
        
        # Transform distances
        print('Applying ' + tForm + ' transformation', flush=True)
        if tForm == 'linear':
            # Linear transform of distances
            ccArr = 1 - (1/dThreshold) * ccArr        
            # Convert negative values  to 0. (These are beyond the threshold)
            ccArr[ccArr < 0] = 0
            # Set nan to 0
            ccArr[np.isnan(ccArr)] = 0
        elif tForm == 'inverse':
            # Inverse square transform of distances
            ccArr = 1/(ccArr + 1)
            # Convert values beyond the distance threshold to 0
            ccArr[ccArr < 1/(dThreshold + 1)] = 0
            # Set nan to 0
            ccArr[np.isnan(ccArr)] = 0
        elif tForm == 'inversesquare':
            # Inverse transform of distances
            ccArr = 1/(ccArr**2 + 1)
            # Convert values beyond the distance threshold to 0
            ccArr[ccArr < 1/(dThreshold + 1)] = 0
            # Set nan to 0
            ccArr[np.isnan(ccArr)] = 0
        elif tForm == 'gaussian':
            dispScale = dThreshold/4
            ccArr = np.exp(-1*((ccArr**2)/(2*(dispScale**2))))
            # Set values beyond the distance threshold to 0
            ccArr[ccArr < np.exp(-1*((dThreshold**2)/(2*(dispScale**2))))] = 0
            # Set nan to 0
            ccArr[np.isnan(ccArr)] = 0
        
        # Multiplier
        if kvol > 1:
            vc = kvol * 3/(np.pi*dThreshold**2)
            ccArr = ccArr*vc
        elif kvol < 1:
            raise Exception('Kernel volume should be >= 1.')
        
        # Sum kernels
        print('Summing kernels', flush=True)
        crkSum = np.sum(ccArr, axis=0)
        del ccArr

    #%%
    print('Writing kernels to file', flush=True)
    # Create 1D zeros array to hold kernel values
    dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    # Update with kernel values
    dArr[nodeids] = crkSum
    del crkSum
    # Reshape into rectangular array
    dArr = dArr.reshape(r.shape[0],r.shape[1])

    # Smooth (not implemented yet)
    #from scipy.ndimage import gaussian_filter
    #dArr = gaussian_filter(dArr, sigma=1, radius=2)
    
    #%% 
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    dArr = np.expand_dims(dArr, axis=0)
    
    # Write to file
    cf.arrayToGeoTiff(dArr, ofile, profile)
    
    toc = time.perf_counter()
    print(f"Calculating kernels took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
