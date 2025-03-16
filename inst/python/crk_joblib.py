# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 09:11:24 2025
Script to implement cumulative resistant kernel mapping using Networkit.
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
import tables as tb
from joblib import Parallel, delayed

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
    
    # Kernel shape transform. Can be either linear, gaussian, inverse, or inversesquare
    tForm = sys.argv[5] # linear
    
    # Should kernel volume be transformed? yes or no. Default is no.
    tkv = sys.argv[6]
    
    # Kernel volume. Ignored if tkv is no.
    # If tkv is yes, then use user input.
    kvol = sys.argv[7]
    
    # Number of threads to use. 
    nThreads = sys.argv[8] # Default 1.

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[9] # Default None
  
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
    
    # Convert threshold to float or integer
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
    # Create edges, nodeids, and mapping between old and new nodeids
    # from resistance grid
    edges, nodeids, idmap = cf.generate_edges(r, cSize)
    print('generated edges')

    # Create graph (nodes only)
    nkG = nk.Graph(len(idmap), weighted=True)
    # Add edges to graph
    for i in edges:
        nkG.addEdge(i[0], i[1], w=i[2], addMissing=False, checkMultiEdge=False)
    print('created graph')
    print(nk.overview(nkG))
    
    # Convert resistance grid to graph
#    print("Converting image to graph", flush=True)
#    nkG, nodeids, idmap = cf.image_to_graph(r, cSize, -9999, 8)
#    print(nk.overview(nkG))
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
            checkND = r[np.array(cinds[:,0], dtype=int),np.array(cinds[:,1], dtype=int)]
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
    # Calculate GB required by pairwise source points
    # number of nodes/750 is an empirical factor
    gbCf = 0.000000000125 # bits to gb conversion factor
    memReq = (len(sources)**2)*64*gbCf*(nkG.numberOfNodes()/750)
    
    # Use only 80% of user supplied threshold to be conservative
    gbThreshold = gbLim*0.8

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
    # Estimate memory required for kernel calculation by multiplying
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

        # Number of batches 
        nBatches = int(np.floor(nodeidsLen/(arraysProc)))

        # Use 50 percent of the number of lines across all threads
        # as a way to allocate memory to the output list created
        # by joblib.
        jlibIncrement = int(np.ceil((nodeidsLen/nBatches*nThreads*0.5)/nThreads))
        
        # Split nodeids into N batches
        sBatches = np.array_split(nodeids, nBatches)
        
        # Use joblib to calculate kernels in parallel and sum on the fly
        print("Summing kernels", flush=True)
        with Parallel(n_jobs=nThreads) as parallel:
            ccArr = np.zeros((1,nodeidsLen))
            n_iter = 0
            for j in range(jlibIncrement,int((np.ceil(nBatches/jlibIncrement)+1)*jlibIncrement),jlibIncrement):
                j = min(j, nBatches)
                results = parallel(delayed(cf.calcKernels)(i, nodeidsLen, dahdf, sBatches, dThreshold, tForm, tkv, kvol) for i in range(n_iter,j))
                ccArr += sum(results)
                n_iter = j  

  
    
    # Create 1D zeros array to hold kernel values
    dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    # Update with kernel values
    dArr[nodeids] = ccArr
    del ccArr
    # Reshape into rectangular array
    dArr = dArr.reshape(r.shape[0],r.shape[1])

    # Calculate quantiles of values > 0 and write to csv
    qDf = np.quantile(dArr[dArr>0],np.arange(0,1.01,0.01))
    qDf = pd.DataFrame({'q': np.arange(0,1.01,0.01), 'value': qDf})
    # Write quantiles to csv
    qDf.to_csv(ofile.replace(".tif","_quantiles.csv"),index=False)
    
    # Smooth (not implemented yet)
    #from scipy.ndimage import gaussian_filter
    #dArr = gaussian_filter(dArr, sigma=1, radius=2)
    
    #%% 
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    dArr = np.expand_dims(dArr, axis=0)
    
    # Write crk to file
    cf.arrayToGeoTiff(dArr, ofile, profile)    
    
    toc = time.perf_counter()
    print(f"Calculating kernels took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
