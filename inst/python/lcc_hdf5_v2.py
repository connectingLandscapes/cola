# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 11:43:39 2023

@author: pj276
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:35:02 2023

@author: pj276
"""

#%%
# IMPORTS
import sys
import cola_functions as cf
import networkit as nk
import rasterio as rio
from rasterio.crs import CRS
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter
import time
from itertools import combinations
import h5py
from numba import set_num_threads

def main() -> None:
    #%%
    # Takes a grouped df and an integer number of max nodes
    # per batch. Groups with size >= maxL should have already
    # been removed from the df.
    def batchesCalc(grup, maxL):
        bigList = []
        uniqueNodes = []
        ncount = 0        
        for name,group in grup:
            while ncount < maxL:
                uniqueNodes = uniqueNodes + np.unique(group.to_numpy()).tolist()
                nuns = list(set(uniqueNodes)) # unique node ids
                ncount = len(uniqueNodes)
            bigList.append(nuns)
            ncount = 0
            uniqueNodes = []
            nuns = []
        return bigList

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
    
    # Output pairwise point hdf name
    pphdf = sys.argv[9]
    
    # Output distance array hdf name
    dahdf = sys.argv[10]

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
    # Otherwise, read as usual
    else:            
        with rio.open(rg) as src:
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
    
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
    # Calculate GB required by pairwise source points
    # Not sure why but the pairwise distance calculation 
    # requires more memory than expected. Set the limit 
    # here so that problems with > ~1000 points get split up
    gbCf = 0.000000000125 # bits to gb conversion factor
    memReq = (len(sources)**2)*64*gbCf
    
    #%%
    # Split into batches if too many source points
    gbPLimit = .008
    #gbPLimit = 7e-05
    if memReq > gbPLimit:
        # Create hdf file to hold distance array
        f = h5py.File(pphdf, mode='a')
        f.create_dataset("dmat", (len(sources), len(sources)),
                         compression="lzf", shuffle=True,
                         chunks=True, dtype='float32')
        f.close()
        # Calculate number of batches
        nPBatches = np.ceil(memReq/gbPLimit)        
        sLength = np.arange(len(sources))
        sBatches = np.array_split(sLength, nPBatches)
        print('Calculating pairwise distances between sources in ' + str(len(sBatches)) + ' batches', flush=True)
        for s in sBatches:
            # Calculate shortest path distance from each source to
            # every other source            
            spspDist = nk.distance.SPSP(nkG, sources[s[0]:s[-1]+1], sources)
            spspDist.run()
            ppArr = spspDist.getDistances(asarray=True)
            ppArr[ppArr == np.finfo(np.float64).max] = np.finfo(np.float32).max
            ppArr = ppArr.astype('float32')
            del spspDist
            # Insert in hdf file
            with h5py.File(pphdf, 'a', ) as f:
                f['dmat'][np.min(s):(np.max(s)+1),:] = ppArr
            del ppArr
        #%%
        # If threshold > 0, apply threshold
        if dThreshold > 0:
            # List to hold pairwise neighbor info
            thBigList = []
            for s in sBatches:
                # Read in array from hdf and append to list
                # Get distance array
                f = h5py.File(pphdf, mode='r')
                ppArr = f['dmat'][s,:]
                f.close()
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
                f = h5py.File(pphdf, mode='r')
                ppArr = f['dmat'][s,:]
                f.close()
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
    # Otherwise, process as normal
    else:
        #%%
        # Calculate shortest path distance from each source to
        # every other source
        print('Calculating pairwise distances between sources')
        spspDist = nk.distance.SPSP(nkG, sources, sources)
        spspDist.run()
        ppArr = spspDist.getDistances(asarray=True)
        ppArr[ppArr == np.finfo(np.float64).max] = np.finfo(np.float32).max
        ppArr = ppArr.astype('float32')
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
    # Estimate memory required for corridor extraction by multiplying
    # the number of graph nodes by the number of sources
    memReq = len(nodeids)*len(sources)*32*gbCf
    # If it's more than the threshold amount, divide into batches
    gbCLimit = 8
    #gbCLimit = 1.5
    if memReq > gbCLimit:
        # Create hdf file to hold distance array
        f = h5py.File(dahdf, mode='a')
        f.create_dataset("dmat", (len(sources), len(nodeids)),
                         compression="lzf", shuffle=True,
                         chunks=True, dtype='float32')
        f.close()
        
        # Divide sources into batches
        sLength = np.arange(0,len(sources))
        nCBatches = int(np.ceil(memReq/gbCLimit))
        sBatches = np.array_split(sLength, nCBatches) 
        print('Calculating cost distances in ' + str(len(sBatches)) + ' batches')
      
        # Loop over batches and calculate distance array
        for b, s in enumerate(sBatches):
            tic1 = time.perf_counter()
            print('Working on batch ' + str(b+1) + ' of ' + str(len(sBatches)))
            sourceBatchNK = np.array(sources)[s]
            spspDist = nk.distance.SPSP(nkG, sourceBatchNK)
            spspDist.run()
            ccArr = spspDist.getDistances(asarray=True)
            del spspDist
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
            ccArr[ccArr == np.finfo(np.float64).max] = np.nan
            ccArr = ccArr.astype('float32')
            ccArr[ccArr > dThreshold] = np.nan
            # Insert in hdf file
            with h5py.File(dahdf, 'a', ) as f:
                f['dmat'][np.min(s):(np.max(s)+1),:] = ccArr
            del ccArr
            toc1 = time.perf_counter()
            print(f"Writing batch to file took {toc1 - tic1:0.4f} seconds")

        #%%
        # Recalculate 32 bit memory required (data from
        # hdf is 32 bit)
        memReq = len(nodeids)*len(sources)*32*gbCf
        if memReq >= gbCLimit:
    
            # Create zeros array to hold corridor sums
            lccSum = np.zeros((len(nodeids),))

            # Indices of completed pairs
            icp = []
            
            # Set up for matching completed pairs
            nrows, ncols = reOrder.shape
            ddtype={'names':['f{}'.format(i) for i in range(ncols)],
                   'formats':ncols * [reOrder.dtype]}
            # Set up batches
            sLength = np.arange(0,len(sources))
            nCBatches = int(np.ceil(memReq/gbCLimit))
            sBatches = np.array_split(sLength, nCBatches) 
            
            ########
            # 1st stage pair calculations
            # Loop over batches and calculate corridors
            for b, s in enumerate(sBatches):
                s = s.astype(int)
                tic1 = time.perf_counter()
                print('Working on 1st stage batch ' + str(b+1) + ' of ' + str(nCBatches))
                # Get unique values of sources 
                # in batch of source pairs
                sourceBatch = np.unique(s)
                # Get indices of completed pairs
                pncx = list(combinations(sourceBatch, 2))
                pncx = [np.array(i) for i in pncx]
                pncx = np.vstack(pncx)
                nmx = np.intersect1d(reOrder.view(ddtype), pncx.view(ddtype), return_indices=True)
                icp.append(nmx[1])
                # Get distance array
                f = h5py.File(dahdf, mode='r')
                darr = f['dmat'][sourceBatch,:]
                f.close()
                # Iterate through target nodes and calculate paths to nodes within the threshold distance
                # plus the tolerance
                print("Summing corridors", flush=True)
                set_num_threads(nThreads)
                # Get pairs
                reOss = np.vstack([np.array([t[0],t[1]]) for t in nmx[0]]).astype('int')
                lcc = cf.sumLccBatch2Numba(reOss, sourceBatch, darr, corrTolerance)
                del darr
                lccSum += lcc
                del lcc

                toc1 = time.perf_counter()
                print("Calculating batch " + str(b+1) + " took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)
       
            ########
            # 2nd stage pair calculations
            # Shift sBatches by 1/2 batch size
            # This becomes like a new sourceBatches but with indices
            # shifted to pick up new combinations of sources.
            x = np.arange(len(sources)+1)
            y = np.array_split(x, nCBatches)
            y = [(a+(len(sBatches[0])/2)).astype('int') for a in y]
            # Delete indices off the end that are > length sources
            yy = y[-1][y[-1] > len(sources)-1]
            y[-1] = np.delete(y[-1], np.nonzero(np.in1d(y[-1],yy))[0])
            for b, s in enumerate(y):
                s = s.astype(int)
                tic1 = time.perf_counter()
                print('Working on 2nd stage batch ' + str(b+1) + ' of ' + str(nCBatches))
                # Get unique values of sources 
                # in batch of source pairs
                sourceBatch = np.unique(s)
                # Get indices of pairs
                pncx = list(combinations(sourceBatch, 2))
                pncx = [np.array(i) for i in pncx]
                pncx = np.vstack(pncx)
                nmx = np.intersect1d(reOrder.view(ddtype), pncx.view(ddtype), return_indices=True)
                # Check for indices of pairs processed in stage 1 
                # and filter these out
                delInds = np.intersect1d(np.hstack(icp), nmx[1], return_indices=True)
                nmxProcFilt = np.delete(nmx[1], delInds[2])
                # Only add indices for pairs that haven't been processed
                icp.append(nmxProcFilt)
                # Filter out pairs themselves that have been processed
                nmxOFilt = np.delete(nmx[0], delInds[2])
                # Get distance array
                f = h5py.File(dahdf, mode='r')
                darr = f['dmat'][sourceBatch,:]
                f.close()
                # Iterate through target nodes and calculate paths to nodes within the threshold distance
                # plus the tolerance
                print("Summing corridors", flush=True)
                set_num_threads(nThreads)
                # Get pairs and process
                reOss = np.vstack([np.array([t[0],t[1]]) for t in nmxOFilt]).astype('int')
                lcc = cf.sumLccBatch2Numba(reOss, sourceBatch, darr, corrTolerance)
                del darr
                lccSum += lcc
                del lcc

                toc1 = time.perf_counter()
                print("Calculating batch " + str(b+1) + " took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)

            ########
            # 3rd stage pair calculations
            # Group remaining pairs by first id and
            # process groups that are larger than the
            # batch limit
            
            # Delete already processed pairs
            reOremain = reOrder
            reOremain = np.delete(reOremain, np.hstack(icp), axis=0)
            # Group by first id
            grp = pd.DataFrame(reOremain)
            grp = grp.rename(columns={0: "s1", 1: "s2"})
            # Get groups larger than the threshold number of nodes
            grp1 = grp.groupby('s1').filter(lambda x: len(x) >= len(sBatches[0])).groupby('s1')

            # Loop through groups, split into batches
            # then process
            for name, group in grp1:
                uNodes = np.unique(group.to_numpy())
                ng3Batches = int(np.ceil(len(group)/len(sBatches[0])))
                g3Batches = np.array_split(uNodes, ng3Batches)                
                for b, s in enumerate(g3Batches):
                    s = s.astype(int)
                    tic1 = time.perf_counter()
                    print('Working on 3rd stage batch ' + str(b+1) + ' of ' + str(ng3Batches))
                    # Get unique values of sources 
                    # in batch of source pairs
                    sourceBatch = np.unique(s)
                    # Get indices of completed pairs
                    pncx = list(combinations(sourceBatch, 2))
                    pncx = [np.array(i) for i in pncx]
                    pncx = np.vstack(pncx)
                    nmx = np.intersect1d(reOrder.view(ddtype), pncx.view(ddtype), return_indices=True)
                    # Check for indices of pairs processed in stage 1 or 2 
                    # and filter these out
                    delInds = np.intersect1d(np.hstack(icp), nmx[1], return_indices=True)
                    nmxProcFilt = np.delete(nmx[1], delInds[2])
                    # Only add indices for pairs that haven't been processed
                    icp.append(nmxProcFilt)
                    # Filter out pairs themselves that have been processed
                    nmxOFilt = np.delete(nmx[0], delInds[2])
                    # Get distance array
                    f = h5py.File(dahdf, mode='r')
                    darr = f['dmat'][sourceBatch,:]
                    f.close()
                    # Iterate through target nodes and calculate paths to nodes within the threshold distance
                    # plus the tolerance
                    print("Summing corridors", flush=True)
                    set_num_threads(nThreads)
                    # Get pairs and process
                    reOss = np.vstack([np.array([t[0],t[1]]) for t in nmxOFilt]).astype('int')
                    lcc = cf.sumLccBatch2Numba(reOss, sourceBatch, darr, corrTolerance)
                    del darr
                    lccSum += lcc
                    del lcc

                    toc1 = time.perf_counter()
                    print("Calculating batch " + str(b+1) + " took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)
                

            ########
            # 4th stage pair calculations
            # Group remaining pairs by first id and
            # process groups that are smaller than the 
            # batch limit
            def batchesCalc4(grup, maxL):
                bigList = []
                uniqueNodes = []
                ncount = 0  
                gNames = []
                gNameSmall = []
                for name, group in grup:
                    uniqueNodes = uniqueNodes + np.unique(group.to_numpy()).tolist()
                    nuns = list(set(uniqueNodes)) # unique node ids
                    ncount = len(nuns)
                    if ncount < maxL:
                        gNameSmall.append(name)
                    else:
                        bigList.append(nuns)
                        ncount = 0
                        uniqueNodes = []
                        nuns = []
                        gNames.append(gNameSmall)
                        gNameSmall = []
                return (bigList, gNames)
            grp2 = grp.groupby('s1').filter(lambda x: len(x) < 300).groupby('s1')
            # This returns a list of the unique nodes in each batch
            # and a list of the batch names in each batch
            aaa, bbb = batchesCalc4(grp2, 300)

            gbatch = []
            for gid in bbb[1]:
                gbatch.append(grp2.get_group(gid))
            gbatch = pd.concat(gbatch)
            gbatch = gbatch.to_numpy()

            len(np.unique(gbatch))
    # Otherwise, process as normal
    else:
        #%%
        # Calculate shortest path distance from each source to every
        # cell in the landscape (slightly different than above spsp calcs)
        print('Calculating corridors in 1 batch')
        spspDist = nk.distance.SPSP(nkG, sources)
        spspDist.run()
        ccArr = spspDist.getDistances(asarray=True)
        del nkG
        # Networkit gives inaccessible nodes the max float 64 value.
        # Set these to nan.
        ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        ccArr = ccArr.astype('float32')
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
    # Update with kernel values
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
    
    toc = time.perf_counter()
    print(f"Calculating least cost corridors took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()