# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 16:40:13 2023
Script to calculate pairwise cost distance matrix
using Networkit and plot histogram of distances.
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
import matplotlib.pyplot as plt

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()

    # INPUTS
    # Path to file holding xy coordinates
    xyf = sys.argv[1] # 'C:/Users/pj276/Projects/CDPOP_arch/data/sabah_example1_1686679883/batchrun0mcrun0/XY200.csv'
    
    # Path to resistance grid
    rg = sys.argv[2] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example_pro.tif'
    
    # Number of processors
    nThreads = sys.argv[3]
    
    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[4] # Default None
    
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
    else:
        raise Exception('Source point file should have a .csv or .shp extension.')
    
    #%%
    # Read resistance grid to array
    # Assign crs if provided by user
    if upCRS != 'None':
        with rio.open(rg, 'r+') as src:
            # Assign projection
            src.crs = CRS.from_string(upCRS)
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype in ['int8','uint8','intp','uintp']:
                r = r.astype('int16')
                profile['dtype'] = 'int16'
    # Otherwise, read as usual
    else:            
        with rio.open(rg) as src:
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype in ['int8','uint8','intp','uintp']:
                r = r.astype('int16')
                profile['dtype'] = 'int16'
    
    # Check no data value
    if profile['nodata'] != -9999:
        r[r==profile['nodata']] = -9999
        profile.update({'nodata':-9999})
        print('Converting no data values to -9999.')
        #raise Exception('No data value is not equal to -9999. Reformat resistance grid so that no data values are set to -9999.')
    
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
    # Calculate shortest path distance from each source to
    # every other source
    spspDist = nk.distance.SPSP(nkG, sources, sources)
    spspDist.run()
    ccArr = spspDist.getDistances(asarray=True)

    # Networkit gives inaccessible nodes the max float 64 value.
    # Set these to nan.
    ccArr[ccArr == np.finfo(np.float64).max] = np.nan
    ccArr = ccArr.astype('float32')
   
    # Get upper triangle (above the diagonal to avoid zeros)
    ccArrUt = ccArr[np.triu_indices(ccArr.shape[0], k = 1)]
    
    # Set zero values to nan. Apparently networkit gives a value
    # of zero as the distance between two unreachable nodes.
    ccArrUt[ccArrUt == 0] = np.nan
    
    # Plot histogram
    # Creating histogram
    fig, ax = plt.subplots(1, 1)
    ax.hist(ccArrUt, bins='auto')
      
    # Set title
    #ax.set_title("Frequency of distances")
      
    # adding labels
    ax.set_xlabel('Cost Distance', fontsize=16)
    ax.set_ylabel('Count', fontsize=16)
    ax = plt.hist(ccArrUt, bins='auto')
    #ax.set_xlabel('Cost Distance', fontsize=16)
    #ax.set_ylabel('Count', fontsize=16)
    #ax.get_legend().remove()
    
    print('Minimum cost distance is ' + str(np.nanmin(ccArrUt)))
    print('Maximum cost distance is ' + str(np.nanmax(ccArrUt)))
    print('Mean cost distance is ' + str(np.nanmean(ccArrUt)))
    
    toc = time.perf_counter()
    print(f"Calculating distances took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()