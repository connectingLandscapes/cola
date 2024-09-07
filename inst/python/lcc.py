# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Script to calculate pairwise cost distance corridors using Networkit.
Runs all source points at the same time with no initial distance
threshold. This can produce large arrays so isn't as memory
safe but is much faster than iterating through each source
point.
@author: pj276
"""

#%%
# IMPORTS
import sys
#import osgeo
import cola_functions as cf
import networkit as nk
import rasterio as rio
#from rasterio.crs import CRS
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter
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

    # Radius for gaussian smoother (in number of cells)
    # The size of the kernel on each side is 2*radius + 1
    # E.g. a radius of 2 gives a 5x5 cell kernel
    gRad = sys.argv[5] # 2

    # Amount to add to the least cost path (in cost distance units)
    # in order to generate a swath of low cost pixels,
    # termed the least cost corridor.
    # If 0, returns the least cost path.
    # If > 0, this amount is added to the least cost path
    # value so that all pixels with values <= to that value
    # will be returned. This results in a swath of pixels
    # instead of a single pixel wide path.
    corrTolerance = sys.argv[6] # 1000

    # Number of threads to use
    nThreads = sys.argv[7] # Default 1

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[8] # Default None
    
    # Convert distance threshold to float or integer
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
    # Set number of processors
    nk.setNumberOfThreads(nThreads)
    
    #%%
    # Read xy file to dataframe
    if Path(xyf).suffix == '.csv':
        xy = pd.read_csv(xyf)
    elif Path(xyf).suffix == '.xy':
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
            # Get number of points supplied by user
            nPts = len(cinds)
            # Check for points and remove any that are out of raster bounds
            # Apparently this can happen using cell_indices_from_coords
            # Outside outer dims
            if any(cinds[:,0] >= r.shape[0]):
                cinds = cinds[cinds[:,0] < r.shape[0]]
            if any(cinds[:,1] >= r.shape[1]):
                cinds = cinds[cinds[:,1] < r.shape[1]]
            # Outside inner dims
            if any(cinds[:,0] < 0):
                cinds = cinds[cinds[:,0] >= 0]
            if any(cinds[:,1] < 0):
                cinds = cinds[cinds[:,1] >= 0]
            # Check for individual points with no data
            checkND = r[np.array(cinds[:,0]),np.array(cinds[:,1])]
            if -9999 in checkND:
                cinds = cinds[checkND != -9999,:]
            # If there are no points after removing no data vals,
            # exit script.
            if len(cinds) == 0:
                raise Exception('No valid source points. This may be because none of your source points overlap with your resistance grid.')
            # Check if any points were removed
            if nPts - len(cinds) > 0:
                print('Ignoring ' + str(nPts - len(cinds)) + " source point(s) with no data values.", flush=True)
            # Get unique indices
            cinds = np.unique(cinds,axis=0)
            # Convert to list
            sources = [(i[0]*r.shape[1])+i[1] for i in cinds]
            # Convert to networkit nodeids
            sources = [idmap[n] for n in sources]
            print("Number of valid sources = " + str(len(sources)), flush=True)

    #%%
    # Calculate shortest path distance from each source to
    # every other source
    spspDist = nk.distance.SPSP(nkG, sources, sources)
    spspDist.run()
    ccList = spspDist.getDistances()
    
    #%%
    # If threshold > 0, apply threshold
    if dThreshold > 0:
        # Get nodes that are less than the threshold distance from the target node
        thList = []
        for lst in ccList:
            thList.append([i for i,j in enumerate(lst) if (j <= dThreshold) and (j > 0)])
        # Need to iterate through sublists to check if there are any point pairs
        # within the threshold distance.
        if np.max([len(i) for i in thList]) == 0:
            raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
    else:
        thList = lst

    #%%
    # Get unique source target pairs
    reOrder = cf.sourceTargetPairs(list(range(len(thList))), thList)
    
    #%%
    # Calculate shortest path distance from each source to every
    # cell in the landscape (slightly different than above spsp calcs)
    spspDist = nk.distance.SPSP(nkG, sources)
    spspDist.run()
    ccArr = spspDist.getDistances(asarray=True)
    # Networkit gives inaccessible nodes the max float 64 value.
    # Set these to nan.
    ccArr[ccArr == np.finfo(np.float64).max] = np.nan
    
    #%%
    # Iterate through target nodes and calculate paths to nodes within the threshold distance
    # plus the tolerance
    print("Calculating corridors", flush=True)
    # Create zeros array to hold corridor sums
    lccSum = np.zeros((ccArr.shape[1],))
    # Iterate through source target pairs and calculate distances
    for i in reOrder:
        source = i[0]
        target = i[1]
        # Get conditional minimum transit cost
        lcc = ccArr[source,:]+ccArr[target,:]
        # Set any inf values to nan
        lcc[np.isinf(lcc)] = np.nan
        # Get least cost path value
        lcpVal = np.nanmin(lcc)
        # Convert to 0/1 corridor
        lcc = np.where(lcc <= lcpVal + corrTolerance, 1, 0)
        lcc[np.isnan(lcc)] = 0
        lccSum = lcc + lccSum
    del ccList
    
    # Create 1D zeros array to hold lcc values
    dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    # Update with kernel values
    dArr[nodeids] = lccSum
    #del ccArr
    # Reshape into rectangular array
    dArr = dArr.reshape(r.shape[0],r.shape[1])

    #%%
    # Use gaussian filter to smooth least cost values
    # Sigma is standard deviation of the kernel
    # The size of the kernel on each side is 2*radius + 1
    # E.g. a radius of 2 gives a 5x5 cell kernel
    # If radius is 0, don't smooth
    # Use rule of thumb for setting sigma from here:
    # https://dsp.stackexchange.com/questions/10057/gaussian-blur-standard-deviation-radius-and-kernel-size
    if gRad == 0:
        lccSmooth = dArr
    else:
        sigma = (gRad-1)/4
        # Sigma should be at least 1
        sigma = np.max([1,sigma])
        lccSmooth = gaussian_filter(dArr, sigma=sigma, radius=gRad)

    #%% 
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    lccSmooth = np.expand_dims(lccSmooth, axis=0)
    
    # Write to file
    cf.arrayToGeoTiff(lccSmooth, ofile, profile)
    
    toc = time.perf_counter()
    print(f"Calculating least cost corridors took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()