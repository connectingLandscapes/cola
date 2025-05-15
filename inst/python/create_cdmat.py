# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:33:54 2023
Script to calculate pairwise cost distance matrix using Networkit.
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
    
    # Number of processors
    nThreads = sys.argv[5]
    
    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[6] # Default None
    
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
    
    # Check values between 0-1. If detected, this converts those values to
    # 1 and prints a message.
    r, profile = cf.checkRasterVals(r, profile)
    
    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')
    # Get cell size from profile
    cSize = profile['transform'][0]
    
    #%%
    # Create edges, nodeids, and mapping between old and new nodeids
    # from resistance grid
    edges, nodeids, idmap = cf.generate_edges(r, cSize)
    print('Generated edges')

    # Create graph (nodes only)
    nkG = nk.Graph(len(idmap), weighted=True)
    # Add edges to graph
    for i in edges:
        nkG.addEdge(i[0], i[1], w=i[2], addMissing=False, checkMultiEdge=False)
    print('Created graph')
    print('Number of nodes: ' + str(nkG.numberOfNodes()))
    print('Number of edges: ' + str(nkG.numberOfEdges()))
    #print(nk.overview(nkG))
    
    # Convert resistance grid to graph
#    print("Converting image to graph", flush=True)
#    nkG, nodeids, idmap = cf.image_to_graph(r, cSize, -9999, 8)
#    print(nk.overview(nkG))
    
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
    ccArr[ccArr == np.finfo(np.float64).max] = 99887766543211
    #ccArr = ccArr.astype('float32')
    
    # If distance threshold is > zero, apply thresholding
    if dThreshold > 0:
        # Set values beyond threshold value to max float 64 values
        ccArr[ccArr > dThreshold] = 99887766543211
        #ccArr[ccArr > dThreshold] = np.nan

    # Networkit sets the distance between 
    # two inaccessible points to zero.    
    # Find cells where off diagonal elements
    # are 0 and set to max for float64
    ccArr[np.where(np.logical_and(~np.eye(ccArr.shape[0],dtype=bool), ccArr == 0))] = 99887766543211
        
    # Save distance array to csv
    ccArr = pd.DataFrame(ccArr)
    ccArr.to_csv(ofile, header=False, index=False)
    
    toc = time.perf_counter()
    print(f"Calculating cost distance matrix took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()