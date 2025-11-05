# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 14:33:26 2025
Read arrays pairwise from zarr file and calculate corridors
@author: pj276
"""

#%%
# IMPORTS
import sys, os, shutil
import zarr
import cola_functions as cf
#import cola_zarr_functions as czf
import numpy as np
from scipy.ndimage import gaussian_filter
import time
from joblib import Parallel, delayed
from collections import OrderedDict

# Zarr file name
dazarr = sys.argv[3]
# Corridor tolerance
corrTolerance = sys.argv[5]
# Zarr file
nodeidsFile = sys.argv[9]
# Get nodeids
nodeids = np.loadtxt(nodeidsFile, delimiter=',').astype(np.int32)

# Get nodeids length
nodeidsLen = len(nodeids)

# --- LRU Array Cache (used per worker) ---
class ArrayCache:
    def __init__(self, max_size=4):
        self.cache = OrderedDict()
        self.max_size = max_size

    def get(self, idx):
        if idx in self.cache:
            self.cache.move_to_end(idx)
            return self.cache[idx]
        else:
            zf = zarr.open(dazarr, mode="r")
            arr = zf[idx, :]
            if len(self.cache) >= self.max_size:
                self.cache.popitem(last=False)
            self.cache[idx] = arr
            return arr

# --- Worker Function ---
def process_batch(batch_pairs, nNodes, corrTolerance):
    start_time = time.time()
    cache = ArrayCache(max_size=4)
    partial_sum = np.zeros((1,nodeidsLen), dtype=np.float64) #np.zeros((nNodes,), dtype=np.float64)
    counter = 1
    for i, j in batch_pairs:
        lcc = cache.get(i) + cache.get(j)
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        lcc = np.where(lcc <= lcpVal + corrTolerance + 0.001, 1, 0)
        lcc[np.isnan(lcc)] = 0
        partial_sum += lcc
        if counter == 100:
            et = time.time()
            print(f"Time taken to process 100 corridors: {et - start_time:.2f} seconds", flush=True)            
        counter += 1

    end_time = time.time()
    print(f"Time taken to process batch: {end_time - start_time:.2f} seconds", flush=True)
    
    return partial_sum

# --- Worker Function ---
# ptWeights is a dictionary where each point id has an associated weight.
# Weights are the number of individuals in the population assocaited with that
# point. Weights for connected points are multiplied. This is to ensure that
# the corridor value corresponds to the number of trips taken by individuals
# traveling between points.

# For example, if point A has a weight of 3 and point B has a weight of 5,
# the corridor between them will have a value of 15.

def process_batch_weighted(batch_pairs, nNodes, corrTolerance, ptWeights):
    start_time = time.time()
    cache = ArrayCache(max_size=4)
    partial_sum = np.zeros((1,nodeidsLen), dtype=np.float64) #np.zeros((nNodes,), dtype=np.float64)
    counter = 1
    for i, j in batch_pairs:
        lcc = cache.get(i) + cache.get(j)
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        cValue = ptWeights[i] * ptWeights[j]
        lcc = np.where(lcc <= lcpVal + corrTolerance + 0.001, cValue, 0)
        lcc[np.isnan(lcc)] = 0
        partial_sum += lcc
        if counter == 100:
            et = time.time()
            print(f"Time taken to process 100 corridors: {et - start_time:.2f} seconds", flush=True)            
        counter += 1

    end_time = time.time()
    print(f"Time taken to process batch: {end_time - start_time:.2f} seconds", flush=True)
    
    return partial_sum
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
    #dazarr = sys.argv[3]
    
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
    # Get reordered list
    reOrder = np.loadtxt(reOrderFile, delimiter=',').astype(np.int32)
    # Interleave pairs to reduce cache turnover and subsequent disk reads
    tic1 = time.time()
    reOrder = cf.interLeavePairs(reOrder)
    toc1 = time.time()
    print(f"Time taken to interleave pairs: {toc1 - tic1:.2f} seconds", flush=True)
    
    # Subset pairs list if starting and ending indices are specified
    # Use joblib to calculate corridors in parallel and sum on the fly
    if (sci != 'None') and (eci != 'None') and (eci > sci):
        reOrder = reOrder[sci:eci]

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

    # If there are many corridors, do a test run to estimate total run time
    # Run 100 corridors per thread
    if len(reOrder) >= 5000:
        print("You have a large number of corridors. Implementing a test run", flush=True)
        print("to estimate total run time.", flush=True)
        reOrder100 = reOrder[0:(100*nThreads)]
        batch_size100 = 100
        batches = [reOrder100[i:i + batch_size100] for i in range(0, len(reOrder100), batch_size100)]
        tic1 = time.time()
        results = Parallel(n_jobs=nThreads, return_as="generator_unordered", prefer="threads")(
            delayed(process_batch)(batch, nodeidsLen, corrTolerance) for batch in batches
        )
        testSum = np.zeros((1,nodeidsLen))
        # Loop over generator and sum intermediates
        for result in results:
            testSum += result
        toc1 = time.time()
        tDiff = toc1-tic1
        print(f'Processing speed is {tDiff/len(reOrder100)} seconds per corridor', flush=True)
        print(f'Estimated time to run all corridors is {tDiff/len(reOrder100)*len(reOrder)} seconds', flush=True)
        del testSum
        del result
        del reOrder100
        del batch_size100
        del batches
        #testRun = np.zeros((1,nodeidsLen))
        #tRun = []
        #for i in [(0,100),(100,200),(200,300),(300,400)]:    
        #    tic1 = time.time()
        #    aa = process_batch(reOrder[i[0]:i[1]], nodeidsLen, corrTolerance)
        #    testRun += aa
        #    toc1 = time.time()
        #    tDiff = toc1-tic1
        #    tRun.append(tDiff)
        #tRun = np.mean(tRun)
        #print(f'Average time to run 100 corridors is {tRun} seconds')
        #print(f'Estimated time to run all corridors is {tRun/100*len(reOrder)} seconds')
        #del testRun
        #del aa

    # Equally distribute pairs among workers
    batch_size = int(np.floor(len(reOrder)/nThreads))

    # --- Divide Into Batches ---
    batches = [reOrder[i:i + batch_size] for i in range(0, len(reOrder), batch_size)]
    print("----", flush=True)
    print(f"Processing {len(reOrder)} corridors in {len(batches)} batches.", flush=True)
    
    # --- Run In Parallel ---
#    results = Parallel(n_jobs=num_workers, backend='multiprocessing', prefer="processes")(
#        delayed(process_batch)(batch, nodeidsLen, corrTolerance) for batch in batches
#    )
    # Return as an unordered generator so that intermediate results
    # can be summed as they become available
    results = Parallel(n_jobs=nThreads, return_as="generator_unordered", prefer="threads")(
        delayed(process_batch)(batch, nodeidsLen, corrTolerance) for batch in batches
    )

    # --- Final Sum ---
#    lccSum = np.sum(results, axis=0)

    lccSum = np.zeros((1,nodeidsLen))
    # Loop over generator and sum intermediates
    for result in results:
        lccSum += result
    
    # Write to file
    print('Writing corridors to file', flush=True)
    
    # Create 1D zeros array to hold lcc values
    dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    # Update with lcc values
    # Zarr file
    nodeidsFile = sys.argv[9]
    # Get nodeids
    nodeids = np.loadtxt(nodeidsFile, delimiter=',').astype(np.int32)
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
#    if os.path.exists(dazarr):
#        try:
#            shutil.rmtree(dazarr)
#            print(f"Successfully deleted Zarr store at: {dazarr}")
#        except OSError as e:
#            print(f"Error deleting Zarr store at {dazarr}: {e}")
#    else:
#        print(f"Zarr store not found at: {dazarr}")

#    # Delete csv files
#    if os.path.exists(reOrderFile):
#        os.remove(reOrderFile)
#    if os.path.exists(nodeidsFile):
#        os.remove(nodeidsFile)

    toc = time.perf_counter()
    print(f"Calculating least cost corridors took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()

