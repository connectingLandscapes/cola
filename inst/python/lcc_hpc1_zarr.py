# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 20:24:14 2025
Script to calculate pairwise cost distance arrays using Networkit.
If the memory threshold is set low enough, the script calculates
distances from source points to every cell in the map in batches
and writes to an hdf file on disk. It then reads distances from
the hdf file and calculates cost distance corridors in parallel
using joblib. Best performance will be with a solid state drive.
This script only calculates the distance arrays and saves to 
a file. Subsequent corridor mapping is done by lcc_hpc2.py
@author: pj276
"""

#%%
# IMPORTS
import sys
#import osgeo
import cola_functions as cf
import cola_zarr_functions as czf
import networkit as nk
import rasterio as rio
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
import time
import psutil
import os, shutil

def main() -> None:

    #%%
    # Start timer
    tic = time.perf_counter()
    
    # INPUTS
    # Path to file holding xy coordinates
    xyf = sys.argv[1] 
    
    # Path to resistance grid
    rg = sys.argv[2] 

    # Distance threshold
    dThreshold = sys.argv[3]

    # Number of threads to use
    nThreads = sys.argv[4] # Default 1

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[5] # Default None
    
    # Output pairwise point hdf name
    # This is a temporary file and gets deleted at the end of the script,
    # unless there's a script failure.
    ppzarr = sys.argv[6]
    
    # Output distance array hdf name
    # This is a temporary file and gets deleted at the end of the script,
    # unless there's a script failure.
    dazarr = sys.argv[7]
    
    # point pair output file name
    reOrderFile = sys.argv[8]

    # Node ids output
    nodeidsFile = sys.argv[9]
    
    # Set memory size for processing corridors
    # I.e. set to 16 if you want to use 16GB of RAM
    # when processing. Make sure you have enough RAM
    # available when setting this value. Consider
    # the total amount of RAM available on your computer
    # and the amount used by other programs that may
    # be running.
    gbLim = sys.argv[10] # Default 6

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
    process = psutil.Process(os.getpid())
    max_memory_mb = 0
    # Initial memory check
    initial_memory = process.memory_info().rss / (1024 * 1024)  # in MB
    max_memory_mb = initial_memory

    #%%
    # Create edges, nodeids, and mapping between old and new nodeids
    # from resistance grid
    # Get nonzero pixel indices
    nzp = np.nonzero(r > 0)
    # Get nonzero pixel values
    nzpv = r[nzp]
    # Combine into 3 column array (row, column, value)
    iArray = np.vstack((nzp[0],nzp[1],nzpv)).T
    del nzp, nzpv
    # Get original node ids, accounting for gaps associated
    # with non-valid cells
    nCols = r.shape[1]
    nodeids = czf.idSkip(iArray, nCols)
    # Create dictionary mapping between orginal node ids
    # and continuous node ids
    idmap = dict((zip(nodeids, range(len(iArray)))))
    # Create numpy array to serve as key value pairs for 
    # node list with gaps and continuous node list
    nPix = r.shape[0]*r.shape[1]
    keyArray = np.empty((nPix, ), dtype=int)
    keyArray[nodeids] = range(len(iArray))
    # Generator for edges
    edgeGen = czf.generate_edgesMod(r, iArray, cSize, keyArray)
    # Initialize graph
    nkG = nk.Graph(len(iArray), weighted=True)
    del keyArray, iArray
    
    # Iterate over generator and add edges to graph
    print('Adding edges to graph', flush=True)
    for j in edgeGen:
        for i in j[1]:
            nkG.addEdge(i[0], i[1], w=i[2], addMissing=False, checkMultiEdge=False)
    print('Added edges to graph')
    print('Number of nodes: ' + str(nkG.numberOfNodes()))
    print('Number of edges: ' + str(nkG.numberOfEdges()))
    del edgeGen

    # Monitor memory after the function completes
    final_memory = process.memory_info().rss / (1024 * 1024)  # in MB
    max_memory_mb = max(max_memory_mb, final_memory)
    print(f"Initial RAM usage: {initial_memory:.2f} MB")
    print(f"Final RAM usage: {final_memory:.2f} MB")
    print(f"Maximum observed RAM usage: {max_memory_mb:.2f} MB")

    # Write nodeids to file
    np.savetxt(nodeidsFile, nodeids, delimiter=",", fmt="%d")
    
    # Get length of nodeids
    nodeidsLen = len(nodeids)
    
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
    # Calculate GB required by pairwise source points
    # number of nodes/750 is an empirical factor
    gbCf = 0.000000000125 # bits to gb conversion factor
    memReq = (len(sources)**2)*64*gbCf*(nkG.numberOfNodes()/500)
    
    # Use only 80% of user supplied threshold to be conservative
    gbThreshold = gbLim*0.8
    
    #%%
    # Split into batches if too many source points
    # For --pairwise distance-- calculations
    if memReq > gbThreshold:
        # Create zarr array to hold point pair distance array
        shape = (len(sources), len(sources))
        chunks = (1, len(sources))
        ppz = czf.createZarrBlosc(shape, np.float64, chunks, ppzarr)
        
        # Calculate number of batches
        nPBatches = int(np.ceil(memReq/gbThreshold))
        sLength = np.arange(len(sources))
        sBatches = np.array_split(sLength, nPBatches)
        # Remove any empty batches
        sBatches = [sub_array for sub_array in sBatches if sub_array.size > 0]
        print('Calculating pairwise distances between sources in ' + str(len(sBatches)) + ' batches', flush=True)
        counter = 1
        for s in sBatches:
            tic1 = time.perf_counter()
            # Calculate shortest path distance from each source to
            # every other source
            # Networkit outputs a 64bit array all at once and this
            # can swamp RAM, resulting in a value error.
            try:            
                spspDist = nk.distance.SPSP(nkG, sources[s[0]:s[-1]+1], sources)
                spspDist.run()
                ppArr = spspDist.getDistances(asarray=True)
            except ValueError:
                # print eror
                print("""Calculating corridors failed.
                      The process may have run out of memory.
                      Consider setting a lower memory threshold.
                      Exiting program.""")
                sys.exit(1)    
            del spspDist
            # Insert in zarr file
            ppz[np.min(s):(np.max(s)+1),:] = ppArr
            del ppArr
            print(f"finished batch {counter}", flush=True)
            toc1 = time.perf_counter()
            print(f"Batch {counter} took {toc1 - tic1:0.4f} seconds", flush=True)
            counter += 1

        #%%
        # If threshold > 0, apply threshold
        if dThreshold > 0:
            print('Thresholding pairwise distances', flush=True)
            # List to hold pairwise neighbor info
            thBigList = []
            for s in sBatches:
                # Read in distance array from zarr
                ppArr = ppz[s,:]                
                # Get nodes that are less than the threshold distance from the target node
                thList = [np.where(np.logical_and(lst<= dThreshold, lst > 0)) for lst in ppArr]
                #thList = []
                #for lst in ppArr:
                #    thList.append(np.where(np.logical_and(lst<= dThreshold, lst > 0)))
                del ppArr
                # Need to iterate through sublists to check if there are any point pairs
                # within the threshold distance.
                if np.max([len(i[0]) for i in thList]) == 0:
                    raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
                thBigList.append(thList)
            thList = [val for sublist in thBigList for val in sublist]
            del thBigList
        else:
            # List to hold pairwise neighbor info
            thBigList = []
            for s in sBatches:
                # Read in distance array from zarr
                ppArr = ppz[s,:]
                # Remove self neighbors
                thList = [np.where(lst > 0) for lst in ppArr]
                #thList = []
                #for lst in ppArr:
                #    thList.append(np.where(lst > 0))
                del ppArr
                # Need to iterate through sublists to check if there are any point pairs
                # within the threshold distance.
                if np.max([len(i[0]) for i in thList]) == 0:
                    raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
                thBigList.append(thList)
            thList = [val for sublist in thBigList for val in sublist]
            del thBigList
        del ppz
    # Otherwise, process pairwise distances in memory as normal
    else:
        #%%
        # Calculate shortest path distance from each source to
        # every other source
        print('Calculating pairwise distances between sources')
        try:
            spspDist = nk.distance.SPSP(nkG, sources, sources)
            spspDist.run()
            ppArr = spspDist.getDistances(asarray=True)
        except ValueError:
            # print eror
            print("""Calculating corridors failed.
                  The process may have run out of memory.
                  Consider setting a lower memory threshold.
                  Exiting program.""")
            sys.exit(1)
        del spspDist
        #%%
        # If threshold > 0, apply threshold
        if dThreshold > 0:
            print('Thresholding pairwise distances')
            # Get nodes that are less than the threshold distance from the target node
            thList = [np.where(np.logical_and(lst<= dThreshold, lst > 0)) for lst in ppArr]
            #thList = []
            #for lst in ppArr:
            #    thList.append(np.where(np.logical_and(lst<= dThreshold, lst > 0)))
            del ppArr
            # Need to iterate through sublists to check if there are any point pairs
            # within the threshold distance.
            if np.max([len(i[0]) for i in thList]) == 0:
                raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
        else:
            # Remove self neighbors
            thList = [np.where(lst > 0) for lst in ppArr]
            #thList = []
            #for lst in ppArr:
            #    thList.append(np.where(lst > 0))
            del ppArr
            # Need to iterate through sublists to check if there are any point pairs
            # within the threshold distance.
            if np.max([len(i[0]) for i in thList]) == 0:
                raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
    
    #%%
    # Get unique source target pairs
    reOrder = cf.newSourceTargetPairs(thList)
    del thList
    print('Processing ' + str(len(reOrder)) + ' source target pairs', flush=True)
    # Save to file
    np.savetxt(reOrderFile, reOrder, delimiter=",", fmt="%d")    

    #%%
    # Estimate memory required for distance mapping by multiplying
    # the number of graph nodes by the number of sources
    memReq = nodeidsLen*len(sources)*64*gbCf
    
    # Memory check
    process = psutil.Process(os.getpid())
    max_memory_mb = 0
    # Initial memory check
    initial_memory = process.memory_info().rss / (1024 * 1024)  # in MB
    max_memory_mb = initial_memory
    
    # Monitor memory after the function completes
    final_memory = process.memory_info().rss / (1024 * 1024)  # in MB
    max_memory_mb = max(max_memory_mb, final_memory)
    print(f"Initial RAM usage: {initial_memory:.2f} MB")
    print(f"Final RAM usage: {final_memory:.2f} MB")
    print(f"Maximum observed RAM usage: {max_memory_mb:.2f} MB")
    
    # If it's more than the threshold amount, divide into batches
    # And map the distance between each source point and every cell in the landscape
    # --distance mapping--
    if memReq > gbThreshold:
        # Create zarr array to hold distance arrays
        shape = (len(sources), nodeidsLen)
        chunks = (1, nodeidsLen)
        daz = czf.createZarrBlosc(shape, np.float64, chunks, dazarr)
   
        # Divide sources into batches
        sLength = np.arange(0,len(sources))
        # Multiply number of batches by 2 to ensure batch size is small enough
        # to fit in memory.
        nCBatches = int(np.ceil(memReq/gbThreshold))*2
        sBatches = np.array_split(sLength, nCBatches)
        # Remove any empty batches
        sBatches = [sub_array for sub_array in sBatches if sub_array.size > 0]
        
        print('Calculating cost distances in ' + str(len(sBatches)) + ' batches')
      
        # Loop over batches and calculate distance array
        for b, s in enumerate(sBatches):
            tic1 = time.perf_counter()
            print('Working on batch ' + str(b+1) + ' of ' + str(len(sBatches)))
            print('Batch size is ' + str(len(s)))
            sourceBatchNK = np.array(sources)[s]
            # Memory check
            process = psutil.Process(os.getpid())
            max_memory_mb = 0
            # Initial memory check
            initial_memory = process.memory_info().rss / (1024 * 1024)  # in MB
            max_memory_mb = initial_memory
            try:
                spspDist = nk.distance.SPSP(nkG, sourceBatchNK)
                spspDist.run()
                ccArr = spspDist.getDistances(asarray=True)
            except ValueError:
                # print error
                print("""Calculating corridors failed.
                      The process may have run out of memory.
                      Consider setting a lower memory threshold.
                      Exiting program.""")
                sys.exit(1)
            del spspDist
            # Monitor memory after the function completes
            final_memory = process.memory_info().rss / (1024 * 1024)  # in MB
            max_memory_mb = max(max_memory_mb, final_memory)
            print(f"Initial RAM usage: {initial_memory:.2f} MB")
            print(f"Final RAM usage: {final_memory:.2f} MB")
            print(f"Maximum observed RAM usage: {max_memory_mb:.2f} MB")
            
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
            ccArr[ccArr == np.finfo(np.float64).max] = np.nan
            ccArr[ccArr > dThreshold] = np.nan
            # Insert in zarr file
            daz[np.min(s):(np.max(s)+1),:] = ccArr
            del ccArr
            toc1 = time.perf_counter()
            print(f"Writing batch to file took {toc1 - tic1:0.4f} seconds")

    # If zarr files were created, delete them
    if os.path.exists(ppzarr):
        try:
            shutil.rmtree(ppzarr)
            print(f"Successfully deleted Zarr store at: {ppzarr}")
        except OSError as e:
            print(f"Error deleting Zarr store at {ppzarr}: {e}")
    else:
        print(f"Zarr store not found at: {ppzarr}")

    toc = time.perf_counter()
    print(f"Calculating distance files took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
