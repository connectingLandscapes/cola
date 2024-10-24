# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 16:47:57 2023
Script to calculate pairwise cost distance corridors using Networkit.
If the memory threshold is set low enough, the script calculates
distances from source points to every cell in the map in batches
and writes to an hdf file on disk. It then reads distances from
the hdf file and calculates cost distance corridors in parallel
using joblib. Best performance will be with a solid state drive.
@author: pj276
"""

#%%
# IMPORTS
import sys
import cola_functions as cf
import networkit as nk
import rasterio as rio
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter
import time
import tables as tb
from joblib import Parallel, delayed

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()
    
    # INPUTS
    # Path to file holding xy coordinates
    xyf = sys.argv[1] 
    
    # Path to resistance grid
    rg = sys.argv[2] 

    # Output file path
    ofile = sys.argv[3]

    # Distance threshold
    dThreshold = sys.argv[4]

    # Radius for gaussian smoother (in number of cells)
    # The size of the kernel on each side is 2*radius + 1
    # E.g. a radius of 2 gives a 5x5 cell kernel
    gRad = sys.argv[5]

    # Amount to add to the least cost path (in cost distance units)
    # in order to generate a swath of low cost pixels,
    # termed the least cost corridor.
    # If 0, returns the least cost path.
    # If > 0, this amount is added to the least cost path
    # value so that all pixels with values <= to that value
    # will be returned. This results in a swath of pixels
    # instead of a single pixel wide path.
    corrTolerance = sys.argv[6]

    # Number of threads to use
    nThreads = sys.argv[7] # Default 1

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[8] # Default None
    
    # Output pairwise point hdf name
    # This is a temporary file and gets deleted at the end of the script,
    # unless there's a script failure.
    pphdf = sys.argv[9]
    
    # Output distance array hdf name
    # This is a temporary file and gets deleted at the end of the script,
    # unless there's a script failure.
    dahdf = sys.argv[10]

    # Set memory size for processing corridors
    # I.e. set to 16 if you want to use 16GB of RAM
    # when processing. Make sure you have enough RAM
    # available when setting this value. Consider
    # the total amount of RAM available on your computer
    # and the amount used by other programs that may
    # be running.
    gbLim = sys.argv[11] # Default 6
    
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
    memReq = (len(sources)**2)*64*gbCf*(nkG.numberOfNodes()/750)
    
    # Use only 80% of user supplied threshold to be conservative
    gbThreshold = gbLim*0.8
    
    #%%
    # Split into batches if too many source points
    # For --pairwise distance-- calculations
    if memReq > gbThreshold:
        # Create hdf file to hold point pair distance array
        shape = (len(sources), len(sources))
        atom = tb.Float64Atom()
        filters = tb.Filters(complib='blosc2:lz4', shuffle=True, complevel=5, fletcher32=False)
        h5f = tb.open_file(pphdf, 'w')
        h5f.create_carray(h5f.root, 'dset', atom, shape,
                               filters=filters)
        h5f.close()
        
        # Calculate number of batches
        nPBatches = int(np.ceil(memReq/gbThreshold))
        sLength = np.arange(len(sources))
        sBatches = np.array_split(sLength, nPBatches)
        print('Calculating pairwise distances between sources in ' + str(len(sBatches)) + ' batches', flush=True)
        for s in sBatches:
            # Calculate shortest path distance from each source to
            # every other source            
            spspDist = nk.distance.SPSP(nkG, sources[s[0]:s[-1]+1], sources)
            spspDist.run()
            ppArr = spspDist.getDistances(asarray=True)
            del spspDist
            # Insert in hdf file
            h5f = tb.open_file(pphdf, 'a')
            h5f.root.dset[np.min(s):(np.max(s)+1),:] = ppArr
            h5f.close()
            
        #%%
        # If threshold > 0, apply threshold
        if dThreshold > 0:
            # List to hold pairwise neighbor info
            thBigList = []
            for s in sBatches:
                # Read in array from hdf and append to list
                # Get distance array
                h5f = tb.open_file(pphdf, 'r')
                ppArr = h5f.root.dset[s,:]
                h5f.close()
                print('Thresholding pairwise distances')
                # Get nodes that are less than the threshold distance from the target node
                thList = []
                for lst in ppArr:
                    thList.append(np.where(np.logical_and(lst<= dThreshold, lst > 0)))
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
                # Read in array from hdf and append to list
                # Get distance array
                h5f = tb.open_file(pphdf, 'r')
                ppArr = h5f.root.dset[s,:]
                h5f.close()
                # Remove self neighbors
                thList = []
                for lst in ppArr:
                    thList.append(np.where(lst > 0))
                del ppArr
                # Need to iterate through sublists to check if there are any point pairs
                # within the threshold distance.
                if np.max([len(i[0]) for i in thList]) == 0:
                    raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
                thBigList.append(thList)
            thList = [val for sublist in thBigList for val in sublist]
            del thBigList

    # Otherwise, process pairwise distances in memory as normal
    else:
        #%%
        # Calculate shortest path distance from each source to
        # every other source
        print('Calculating pairwise distances between sources')
        spspDist = nk.distance.SPSP(nkG, sources, sources)
        spspDist.run()
        ppArr = spspDist.getDistances(asarray=True)

        #%%
        # If threshold > 0, apply threshold
        if dThreshold > 0:
            print('Thresholding pairwise distances')
            # Get nodes that are less than the threshold distance from the target node
            thList = []
            for lst in ppArr:
                thList.append(np.where(np.logical_and(lst<= dThreshold, lst > 0)))
            del ppArr
            # Need to iterate through sublists to check if there are any point pairs
            # within the threshold distance.
            if np.max([len(i[0]) for i in thList]) == 0:
                raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
        else:
            # Remove self neighbors
            thList = []
            for lst in ppArr:
                thList.append(np.where(lst > 0))
            del ppArr
            # Need to iterate through sublists to check if there are any point pairs
            # within the threshold distance.
            if np.max([len(i[0]) for i in thList]) == 0:
                raise Exception('No source points are within the threshold distance. You may want to increase the distance threshold.')
      
    #%%
    # Get unique source target pairs
    reOrder = cf.newSourceTargetPairs(thList)
    print('Processing ' + str(len(reOrder)) + ' source target pairs')
    
    #%%
    # Estimate memory required for distance mapping by multiplying
    # the number of graph nodes by the number of sources
    memReq = nodeidsLen*len(sources)*64*gbCf
    
    # If it's more than the threshold amount, divide into batches
    # And map the distance between each source point and every cell in the landscape
    # --distance mapping--
    if memReq > gbThreshold:
        # Create hdf file to hold distance array
        shape = (len(sources), nodeidsLen)
        atom = tb.Float64Atom()
        filters = tb.Filters(complib='blosc2:lz4', shuffle=True, complevel=5, fletcher32=False)
        h5f = tb.open_file(dahdf, 'w')
        h5f.create_carray(h5f.root, 'dset', atom, shape,
                               filters=filters)
        h5f.close()
    
        # Divide sources into batches
        sLength = np.arange(0,len(sources))
        nCBatches = int(np.ceil(memReq/gbThreshold))
        sBatches = np.array_split(sLength, nCBatches) 
        print('Calculating cost distances in ' + str(len(sBatches)) + ' batches')
      
        # Loop over batches and calculate distance array
        for b, s in enumerate(sBatches):
            tic1 = time.perf_counter()
            print('Working on batch ' + str(b+1) + ' of ' + str(len(sBatches)))
            print('Batch size is ' + str(len(s)))
            sourceBatchNK = np.array(sources)[s]
            spspDist = nk.distance.SPSP(nkG, sourceBatchNK)
            spspDist.run()
            ccArr = spspDist.getDistances(asarray=True)
            del spspDist
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
            ccArr[ccArr == np.finfo(np.float64).max] = np.nan
            ccArr[ccArr > dThreshold] = np.nan
            # Insert in hdf file
            h5f = tb.open_file(dahdf, 'a')
            h5f.root.dset[np.min(s):(np.max(s)+1),:] = ccArr
            h5f.close()
            del ccArr
            toc1 = time.perf_counter()
            print(f"Writing batch to file took {toc1 - tic1:0.4f} seconds")
    
    #%%
    # Estimate memory required for corridor calculation by multiplying
    # the number of graph nodes by the number of sources
    memReq = nodeidsLen*len(sources)*64*gbCf

    # If it's more than the threshold amount, divide into batches
    if memReq > gbThreshold:
        #%%
        print('Calculating corridors in batches')

        # Memory per processor
        memProc = gbThreshold/nThreads

        # Arrays per processor
        arraysProc = memProc/(nodeidsLen*64*gbCf)

        # Number of batches (divide by 2 becuase each line of reOrder can
        # represent two different arrays)
        nBatches = int(np.floor(len(reOrder)/(arraysProc/2)))

        # Use 50 percent of the number of lines across all threads
        # as a way to allocate memory to the output list created
        # by joblib.
        jlibIncrement = int(np.ceil((len(reOrder)/nBatches*nThreads*0.5)/nThreads))
        
        # Split reOrder into N batches
        sBatches = np.array_split(reOrder, nBatches)
        
        # Use joblib to calculate corridors in parallel and sum on the fly
        print("Summing corridors", flush=True)
        with Parallel(n_jobs=nThreads) as parallel:
            lccSum = np.zeros((1,nodeidsLen))
            n_iter = 0
            for j in range(jlibIncrement,int((np.ceil(nBatches/jlibIncrement)+1)*jlibIncrement),jlibIncrement):
                j = min(j, nBatches)
                results = parallel(delayed(cf.calcPaths)(i, nodeidsLen, dahdf, sBatches, corrTolerance) for i in range(n_iter,j))
                lccSum += sum(results)
                n_iter = j     

    # Otherwise, process as normal
    else:
        #%%
        # Calculate shortest path distance from each source to every
        # cell in the landscape (slightly different than above spsp calcs)
        print('Calculating corridors')
        spspDist = nk.distance.SPSP(nkG, sources)
        spspDist.run()
        ccArr = spspDist.getDistances(asarray=True)
        del nkG
        # Networkit gives inaccessible nodes the max float 64 value.
        # Set these to nan.
        ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        
        #%%
        # Iterate through target nodes and calculate paths to nodes within the threshold distance
        # plus the tolerance
        print("Summing corridors", flush=True)
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
            lccSum += lcc
        del ccArr
        del lcc
        
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
    
    # If hdfs were created, delete them
    if Path(dahdf).is_file():
        Path(dahdf).unlink()
    if Path(pphdf).is_file():
        Path(pphdf).unlink()

    toc = time.perf_counter()
    print(f"Calculating least cost corridors took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()