# -*- coding: utf-8 -*-
"""
Created on Sun May 11 15:46:08 2025

@author: pj276
"""

#%%
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from numba import njit, prange, set_num_threads
import networkit as nk
import time
import cola_functions as cf
import rasterio as rio

# Set number of threads for numba functions
set_num_threads(5)
nk.setNumberOfThreads(5)

@njit(parallel=True)
def kCalcs(ccArr, tForm, dThreshold, tkv, kvol):
    # Networkit gives inaccessible nodes the max float 64 value.
    # Set these to nan. Need to flatten first so numba is only
    # indexing on one dimension
    orig_shape = ccArr.shape
    ccArr = ccArr.flatten()
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
        # Inverse transform of distances
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
    
    # Transform kernel volume?
    # From UNICOR help guide
    # When „const_kernel_vol‟ is False, then the „kernel_volume‟ parameter
    # is used on the transformed kernel resistant distance following
    # „kernal_volume‟ * 3/(math.pi*kernel  resistant distances^2).
    # When „const_kernel_vol‟ is True, then no volume transformation is applied.
    # **I think the UNICOR help on this might be reversed because this is the code
	#	if const_kernal_vol:		
    #        vol_const = vol_constant * 3/(math.pi*edge_dist**2)
    #    else:
    #        vol_const = 1.
    if tkv == "yes":
        vol_const = kvol * 3/(np.pi*dThreshold**2)
        ccArr = ccArr*vol_const   

    # Reshape 1D array    
    ccArr = ccArr.reshape(orig_shape) 

    # Sum crk values
    ccArr = np.sum(ccArr, axis=0)

    return(ccArr)

@njit(parallel=True)
def kCalcs2(ccArr, tForm, dThreshold, tkv, kvol):
    orig_shape = ccArr.shape
    # Empty array to accumulate sums
    acc = np.zeros((1,orig_shape[1])).flatten()
    
    # Networkit gives inaccessible nodes the max float 64 value.
    # Set these to nan. Need to flatten first so numba is only
    # indexing on one dimension
    for i in prange(ccArr.shape[0]):
        a = ccArr[i]
        a[a == np.finfo(np.float64).max] = np.nan
        
        # Transform distances
        if tForm == 'linear':
            # Linear transform of distances
            a = 1 - (1/dThreshold) * a        
            # Convert negative values  to 0. (These are beyond the threshold)
            a[a < 0] = 0
            # Set nan to 0
            a[np.isnan(a)] = 0
        elif tForm == 'inverse':
            # Inverse transform of distances
            a = 1/(a + 1)
            # Convert values beyond the distance threshold to 0
            a[a < 1/(dThreshold + 1)] = 0
            # Set nan to 0
            a[np.isnan(a)] = 0
        elif tForm == 'inversesquare':
            # Inverse transform of distances
            a = 1/(a**2 + 1)
            # Convert values beyond the distance threshold to 0
            a[a < 1/(dThreshold + 1)] = 0
            # Set nan to 0
            a[np.isnan(a)] = 0
        elif tForm == 'gaussian':
            dispScale = dThreshold/4
            a = np.exp(-1*((a**2)/(2*(dispScale**2))))
            # Set values beyond the distance threshold to 0
            a[a < np.exp(-1*((dThreshold**2)/(2*(dispScale**2))))] = 0
            # Set nan to 0
            a[np.isnan(a)] = 0
        
        # Transform kernel volume?
        # From UNICOR help guide
        # When „const_kernel_vol‟ is False, then the „kernel_volume‟ parameter
        # is used on the transformed kernel resistant distance following
        # „kernal_volume‟ * 3/(math.pi*kernel  resistant distances^2).
        # When „const_kernel_vol‟ is True, then no volume transformation is applied.
        # **I think the UNICOR help on this might be reversed because this is the code
    	#	if const_kernal_vol:		
        #        vol_const = vol_constant * 3/(math.pi*edge_dist**2)
        #    else:
        #        vol_const = 1.
        if tkv == "yes":
            vol_const = kvol * 3/(np.pi*dThreshold**2)
            a = a*vol_const   
   
        # Sum crk values
        acc += a
    
    # Test for loop summation with prange
#    acc = np.zeros((1,orig_shape[1]))
#    n = orig_shape[0]
#    for i in prange(n):
#        acc += ccArr[i].reshape((1,orig_shape[1]))
    
    #return(acc)
    return(acc)



#%%
rg = r'C:\Users\pj276\Downloads\rstudioexport\Accessibility_layer_90m_UTM46N.tif'
#xyf = r"C:\Users\pj276\Downloads\rstudioexport\source_points_scenarioB_60_redcued.shp"
xyf = r"C:\Users\pj276\Downloads\rstudioexport\source_points_scenarioB_60_redcued.shp"
ofile = 'C:/Users/pj276/Downloads/rstudioexport/crktest4.tif'
upCRS = 'None'

#%%
# Read xy file to dataframe
if Path(xyf).suffix == '.csv':
    xy = pd.read_csv(xyf)
elif Path(xyf).suffix == '.shp':
    xy = gpd.read_file(xyf)
    xy = xy.get_coordinates()
else:
    raise Exception('Source point file should have a .csv or .shp extension.')

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
print('generated edges')

# Create graph (nodes only)
nkG = nk.Graph(len(idmap), weighted=True)
# Add edges to graph
for i in edges:
    nkG.addEdge(i[0], i[1], w=i[2], addMissing=False, checkMultiEdge=False)
print('created graph')
#print(nk.overview(nkG))
    
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

tic = time.perf_counter()

sBatches = np.array_split(sources, 200)[1:2]
bigSum = np.zeros((1,len(nodeids)))

for i,b in enumerate(sBatches):
    tic1 = time.perf_counter()
    # Calculate shortest path distance from each source to every
    # cell in the landscape
    spspDist = nk.distance.SPSP(nkG, b)
    spspDist.run()
    ccArr = spspDist.getDistances(asarray=True)
    # Convert distances to probabilities and transform if specified
    ccArr = kCalcs(ccArr, 'linear', 1000000, 'yes', 100000000)  

    # Add to big sum
    bigSum = np.add(bigSum, ccArr)
    print("Finished batch " + str(i))
    toc1 = time.perf_counter()
    print(f"Calculating kernels took {toc1 - tic1:0.4f} seconds", flush=True)

toc = time.perf_counter()
print(f"Calculating kernels took {toc - tic:0.4f} seconds", flush=True)

# Create 1D zeros array to hold kernel values
dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
    
# Update with kernel values
dArr[nodeids] = bigSum

del bigSum
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
    
    
    