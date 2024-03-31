# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:23:54 2024

@author: pj276
"""

import matplotlib.pyplot as plt
import os
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
import numpy as np
import rasterio as rio
from skimage.segmentation import watershed
from skimage.filters.rank import minimum, maximum
from skimage.morphology import binary_erosion
from skimage.morphology import remove_small_objects
from skimage.measure import regionprops
import cola_functions as cf
from pathlib import Path
import pandas as pd
import geopandas as gpd
import itertools
from itertools import chain
import networkit as nk



# Pseudocode
# Actually, actually, maybe the easiest/fastest way is to 
# identify patch cells that are adjacent to a corridor cell,
# group these, constrain the cost surface to be within 
# patches and corridors, then calculate cost distance
# corridors between groups of cells.

# Actually, maybe the best/fastest way is to 
# identify the patch cells that are adjacent to a corridor
# cell, use these as the source/target nodes, constrain
# the cost surface to be within patches and corridors
# then factorial map those or some fraction of them.

# Maybe best/fastest way is to create a new cost surface
# where the patches have a value of 1, the corridors have
# the value of the original cost surface, and outside
# of corridors is no data. This would force corridors
# from patch edges to go through the existing corridors
# and would provide edge to edge distance estimates as well.
# This should be faster and less computation than factorial
# mapping between points in patches.
# Not sure how to threshold this surface though to keep
# only corridors between the target patch pair.

# Threshold CRK surface. 
# Run focal min then focal max with 5x5 cell neighborhood
# to eliminate small patches.
# Distribute N random points proportionally to patches
# based on perimeter. Distributing by area quickly results
# in small patches getting < 1 point allocated because
# if any large patches are present. Would need to allocate
# very large numbers of points in order to get enough
# for the small patches.

# Map patch pairwise factorial corridors from random points.
# Sum corridors
# Overlay each pairwise corridor system on this sum and extract
# corridor strength from it for each patch pair corridor system (PPCN)
# Add data on corridor strength, corridor length, and patch strength
# to a table. Normalize and calculate metrics for each
# PPCN.
# Rank them.
    
q = 0.5
# Footprint
k = np.ones((5,5)).astype('uint8')
nPts = 500

# Footprint
fprint = np.ones([3,3]).astype('uint8')

# Read in crk raster
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk100_test6.tif") as src:
    profile = src.profile
    im = src.read(1)
    #cinds = cf.cell_indices_from_coords(src, im, np.array(xy))

# Read in least cost corridors
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\lcc100_test14.tif") as src:
    profile = src.profile
    lcc = src.read(1)

# Create patches from high values of dispersal kernels
hvP = np.where(im >= np.quantile(im[im>0],q), 1, 0)

# Convert corridors to 0-1
lccB = np.copy(lcc)
lccB = np.where(lccB > 0, 1, 0)

# Mask out corridors in patches
lccB[(hvP==1) & (lccB == 1)] = 0

# Label discrete patches
markers, npatch1 = ndi.label(hvP, structure=fprint)
# Get patch attribute
rprops = regionprops(markers, intensity_image=im)
pLabels = [i.label for i in rprops]
pAreas = [i.area for i in rprops]
pMax = [i.intensity_max for i in rprops]
pMean = [i.intensity_mean for i in rprops]


# Check if only one patch.
if npatch1 == 1:
    raise Exception('There is only one patch. You may want to change the threshold for creating patches.')

# Get moving window max over corridors
mx = maximum(lccB.astype('uint8'), fprint)

# Get patch edges
hvPedges = hvP - binary_erosion(hvP, fprint)

# Get intersection of patch edges and moving window == 1
# These are the target pixels
tPix = mx*hvPedges

# Label discrete groups of corridor patch edge cells
cpecells, npatch2 = ndi.label(tPix, structure=fprint)

# Use watershed algorithm to determine which edge cells
# are connected to edge cells in other patches
# To do this, use edge cells as markers, flood them
# through corridors and patches, then use patches
# to mask out the result. Use the zone overlap
# function on the masked result to determine which
# groups of edge cells are connected by corridors.
# This does allow groups of edge cells to be connected
# to other groups of edge cells from the same patch
# as long as the corridor is external to the patch.
# Corridors that run from one edge of a patch, through
# the patch, to another edge, are not included.
# In the former case, we need to explicitly exclude
# patch self connections via external corridors.

# Create mask using positive values of corridors
# and 0 values of high value patches 
wMask = np.where(np.logical_and(lcc > 0, hvP < 1), True, False)
# Add cpecells back to mask
wMask = np.where(cpecells > 0, True, wMask)
#wMask = np.where(lcc > 0, True, False)

# Euclidean distance between edge cells
#edgeDist = ndi.distance_transform_edt(cpecells)

# Use watershed algo to 'flood' least cost corridors
# so that they're allocated to the closest edge cells
# Can use euclidean distance to nearest patch
# or the least cost corridor values.
# Use negative of both so that high value areas
# are flooded first
labels = watershed(-lcc, cpecells, mask=wMask, connectivity=2)

# Get pairs of edges cell groups connected by corridors
ecc = cf.zoneAdjacency(labels)
# Reorder pairs so that largest valued id is on the right
ecc = np.array([np.min(ecc,axis=1),np.max(ecc,axis=1)]).T
# Filter again for unique pairs
ecc = np.unique(ecc, axis=0)
# Add columns for patch ids
ecc = np.hstack((ecc, np.zeros((len(ecc),2),dtype=int)))

# Create dictionary for patches associated with edge cell groups
d = {}
for i in pLabels:
    msk = np.where(markers==i, 1, 0)
    ucellids = np.unique(msk*cpecells)
    for j in ucellids:
        d[j] = i
del d[0]

# Add patch id to array of edge pairs
for i,j in enumerate(ecc):
    ecc[i,2] = d[j[0]]
    ecc[i,3] = d[j[1]]

# Delete edge pairs in the same patch
ecc = ecc[ecc[:,2] != ecc[:,3],:]

#%%
# Read in original cost surface
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\size7.tif") as src:
    profilecs = src.profile
    ocs = src.read(1)
# Get cell size from profile
cSize = profilecs['transform'][0]
# Convert non-corridor, non-patch cells to no data
newcs = np.maximum(np.where(lcc > 0, 1, 0), np.where(hvP > 0, 1, 0))
newcs = np.where(newcs==0, -9999, 1)
ocs = np.where(newcs==1, ocs, -9999)
# Write to file
ocs = np.expand_dims(ocs, axis=0)
ocsname = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/size7_masked_cs.tif'
cf.arrayToGeoTiff(ocs, ocsname, profilecs)
# Remove extra dimension for further processing
ocs = np.squeeze(ocs)
# Convert resistance grid to graph
print("Converting image to graph", flush=True)
nkG, nodeids, idmap = cf.image_to_graph(ocs, cSize, -9999, 8)
print(nk.overview(nkG))

#%%
# Loop through edge pairs and calculate corridors
# plus stats
# Empty list to hold stats
attList = []
for pp in ecc:
    # Get first group of cells
    ex1 = np.where(cpecells==pp[0], 1, 0)
    # Get row, column indices where values are > 0
    ex1Inds = np.argwhere(ex1 > 0)
    # Convert row, column to coordinates
    xCs1, yCs1 = rio.transform.xy(src.transform, ex1Inds[:,0], ex1Inds[:,1])

    # Get second group of cells
    ex2 = np.where(cpecells==pp[1], 1, 0)
    # Get row, column indices where values are > 0
    ex2Inds = np.argwhere(ex2 > 0)
    # Convert row, column to coordinates
    xCs2, yCs2 = rio.transform.xy(src.transform, ex2Inds[:,0], ex2Inds[:,1])

    # Convert to pandas df
    dfCoords1 = pd.DataFrame({'X': xCs1, 'Y': yCs1})
    # Convert to pandas df
    dfCoords2 = pd.DataFrame({'X': xCs2, 'Y': yCs2})

    # Get 1st set of sources as networkit node ids
    sources1 = cf.sources2nodeids(ocsname, ocs, dfCoords1, idmap)
    # Get 2nd set of sources as networkit node ids
    sources2 = cf.sources2nodeids(ocsname, ocs, dfCoords2, idmap)
    
    # Empty list to hold pair of distance rasters
    pdrs = []
    for j in [sources1, sources2]:
        # Calculate cost distance from each source to every
        # cell in the landscape for 1st set of sources
        spspDist = nk.distance.SPSP(nkG, j)
        spspDist.run()
        ccArr = spspDist.getDistances(asarray=True)
        # Networkit gives inaccessible nodes the max float 64 value.
        # Set these to nan.
        ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        
        # Get min of distances
        ccArr = np.min(ccArr, axis=0)
        # Create 1D zeros array to hold kernel values
        dArr = np.zeros([ocs.shape[0]*ocs.shape[1]], 'float32')
        # Update with kernel values
        dArr[nodeids] = ccArr
        del ccArr
        # Reshape into rectangular array
        pdrs.append(dArr.reshape(ocs.shape[0],ocs.shape[1]))
    
    # Distance between sets of points
    dsp = pdrs[0] + pdrs[1]
    # Threshold
    dsp = np.where(dsp < np.min(dsp[dsp>0])+1000, dsp, 0)
    # Set nan to 0
    dsp[np.isnan(dsp)] = 0
    
    # Minimum cost distance value between sets of points
    mcd = np.min(dsp[dsp>0])
    meancd = np.mean(dsp[dsp>0])
    
    # Mask out corridor where it overlaps with patches
    #dsp = np.where(np.logical_and(hvP > 1, dsp > 0), 0, dsp)
    dsp = np.where(hvP == 1, 0, dsp)
    
    # Get corridor strength values overlapping with corridor
    tCorr = np.where(dsp > 0, lcc, 0)
    
    # Get row, column indices of corridor strength 
    # array where values are > 0
    tCStrengthInds = np.argwhere(tCorr > 0)
    # Convert row, column to coordinates
    xtCS, ytCS = rio.transform.xy(src.transform, tCStrengthInds[:,0], tCStrengthInds[:,1])
    # Get mean of coordinates
    meanxtCS = np.mean(xtCS)
    meanytCS = np.mean(ytCS)
    
    # Get max corridor strength value
    tCorrMax = np.max(tCorr)
    tCorrMean = np.mean(tCorr)
    
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    tCorr = np.expand_dims(tCorr, axis=0)
    
    # Write target corridor to file
    ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\tCorr_" + str(pp[0]) + "_" + str(pp[1]) + ".tif"
    cf.arrayToGeoTiff(tCorr, ofile, profilecs)
    
    # Create point shapefile to hold corridor mean location
    # and attributes.
    
    # Get index of first patch label
    pl1 = pLabels.index(pp[2])
    pl2 = pLabels.index(pp[3])
    # Use index values to retrieve patch attributes
    # mean corridor xcoord, mean corridor ycoord, patch id 1, patch id 2, edge id 1, edge id 2, area patch 1, area patch 2, max patch 1, max patch 2, mean patch 1, mean patch 2
    attList.append(np.array([meanxtCS,meanytCS,pp[2],pp[3],pp[0],pp[1],pAreas[pl1],pAreas[pl2],pMax[pl1],pMax[pl2],pMean[pl1],pMean[pl2],mcd,meancd,tCorrMax,tCorrMean]))

# Combine into single array
attArray = np.vstack(attList)
# Shapefile column names
column_names = ['xco','yco','pid1','pid2','eid1','eid2','parea1','parea2','pmax1','pmax2','pmean1','pmean2','mincost','meancost','maxstrength','meanstrength']
# Convert to dataframe
cAttDf = pd.DataFrame(attArray, columns=column_names)
# Add corridor quality metric
cAttDf['r1'] = cAttDf.pmax1 * cAttDf.pmax2 * cAttDf.maxstrength * 1/cAttDf.mincost
# COnvert to geodataframe
gdf = gpd.GeoDataFrame(cAttDf, geometry=gpd.points_from_xy(cAttDf.xco, cAttDf.yco), crs=src.crs)
# Save to file
oshp = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\tCorr_points.shp"
gdf.to_file(oshp)


#!!!!!
# Get cpecells overlapping with a patch
ex1 = np.where(markers==1, 1, 0)*cpecells
# Get one group of cells from that patch
ex1 = np.where(ex1==1, 1, 0)
# Get row, column indices where values are > 0
ex1Inds = np.argwhere(ex1 > 0)
# Convert row, column to coordinates
xCs1, yCs1 = rio.transform.xy(src.transform, ex1Inds[:,0], ex1Inds[:,1])

# Get cpecells overlapping with another patch
ex2 = np.where(markers==2, 1, 0)*cpecells
# Get one group of cells
ex2 = np.where(ex2==6, 1, 0)
# Get row, column indices where values are > 0
ex2Inds = np.argwhere(ex2 > 0)
# Convert row, column to coordinates
xCs2, yCs2 = rio.transform.xy(src.transform, ex2Inds[:,0], ex2Inds[:,1])

# Convert to pandas df
dfCoords1 = pd.DataFrame({'X': xCs1, 'Y': yCs1})
# Convert to pandas df
dfCoords2 = pd.DataFrame({'X': xCs2, 'Y': yCs2})

# Convert to geodataframe and write to file
#gdf = gpd.GeoDataFrame(dfCoords1, geometry=gpd.points_from_xy(dfCoords1.X, dfCoords1.Y), crs=src.crs)
#ofile = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/delme.shp'
#gdf.to_file(ofile)
   
#----
# Read in original cost surface
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\size7.tif") as src:
    profilecs = src.profile
    ocs = src.read(1)
# Get cell size from profile
cSize = profilecs['transform'][0]
    
# Convert non-corridor, non-patch cells to no data
newcs = np.maximum(np.where(lcc > 0, 1, 0), np.where(hvP > 0, 1, 0))
newcs = np.where(newcs==0, -9999, 1)
ocs = np.where(newcs==1, ocs, -9999)
# Write to file
ocs = np.expand_dims(ocs, axis=0)
ocsname = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/size7_masked_cs_test3.tif'
cf.arrayToGeoTiff(ocs, ocsname, profilecs)
# Remove extra dimension for further processing
ocs = np.squeeze(ocs)

# Convert resistance grid to graph
print("Converting image to graph", flush=True)
nkG, nodeids, idmap = cf.image_to_graph(ocs, cSize, -9999, 8)
print(nk.overview(nkG))

# Get 1st set of sources as networkit node ids
sources1 = cf.sources2nodeids('C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/size7_masked_cs_test.tif', ocs, dfCoords1, idmap)
# Get 2nd set of sources as networkit node ids
sources2 = cf.sources2nodeids('C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/size7_masked_cs_test.tif', ocs, dfCoords2, idmap)

# Empty list to hold pair of distance rasters
pdrs = []
for j in [sources1, sources2]:
    # Calculate cost distance from each source to every
    # cell in the landscape for 1st set of sources
    spspDist = nk.distance.SPSP(nkG, j)
    spspDist.run()
    ccArr = spspDist.getDistances(asarray=True)
    # Networkit gives inaccessible nodes the max float 64 value.
    # Set these to nan.
    ccArr[ccArr == np.finfo(np.float64).max] = np.nan
    
    # Get min of distances
    ccArr = np.min(ccArr, axis=0)
    # Create 1D zeros array to hold kernel values
    dArr = np.zeros([ocs.shape[0]*ocs.shape[1]], 'float32')
    # Update with kernel values
    dArr[nodeids] = ccArr
    del ccArr
    # Reshape into rectangular array
    pdrs.append(dArr.reshape(ocs.shape[0],ocs.shape[1]))

# Distance between sets of points
dsp = pdrs[0] + pdrs[1]
# Threshold
dsp = np.where(dsp < np.min(dsp[dsp>0])+1000, dsp, 0)
# Set nan to 0
dsp[np.isnan(dsp)] = 0

# Minimum cost distance value between sets of points
mcd = np.min(dsp[dsp>0])

# Mask out corridor where it overlaps with patches
#dsp = np.where(np.logical_and(hvP > 1, dsp > 0), 0, dsp)
dsp = np.where(hvP == 1, 0, dsp)

# Get corridor strength values overlapping with corridor
tCorr = np.where(dsp > 0, lcc, 0)

# Get max corridor strength value
tCorrMax = np.max(tCorr)

# Add a dimension to the array (the rasterio profile expects
# a dimension corresponding to number of bands)
tCorr = np.expand_dims(tCorr, axis=0)

# Write target corridor to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\tCorr_test2.tif"
cf.arrayToGeoTiff(tCorr, ofile, profilecs)




#----
# Add a dimension to the array (the rasterio profile expects
# a dimension corresponding to number of bands)
dsp = np.expand_dims(dsp, axis=0)

# Write to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\dsp_test.tif"
cf.arrayToGeoTiff(dsp, ofile, profilecs)
    
    



# !!! Working here
aa = np.array([[1,1,1,2,2,2],[1,1,1,2,2,2],[1,1,1,2,2,2]])
aa = np.array([[1,1,1,2,2,2],[1,1,1,2,2,2],[1,1,1,2,2,2],[3,3,3,4,4,4]])

# Right shift
rs = np.pad(aa,((0,0),(1,0)), mode='constant')[:, :-1]
cc = np.unique(np.vstack((aa[(aa != bb) & (aa != 0) & (bb != 0)],bb[(aa != bb) & (aa != 0) & (bb != 0)])).T)
# Left shift
ls = np.pad(aa,((0,0),(0,1)), mode='constant')[:, 1:]
# Up shift
us = np.pad(aa,((0,1),(0,0)), mode='constant')[1:,:]
# Down shift
ds = np.pad(aa,((1,0),(0,0)), mode='constant')[:-1,:]
# Diagonal upper left
dul = np.pad(aa,((0,1),(0,1)), mode='constant')[1:, 1:]
# Diagonal upper right
dur = np.pad(aa,((0,1),(1,0)), mode='constant')[1:, :-1]
# Diagonal lower left
dll = np.pad(aa,((1,0),(0,1)), mode='constant')[:-1, 1:]
# Diagonal lower right
dlr = np.pad(aa,((1,0),(1,0)), mode='constant')[:-1, :-1]



# Convert tPix to markers
# Calculate euclidean distance from markers along corridors
# Watershed algorithm along corridors
# Moving window to identify which watershed markers
# meet at least one other watershed marker.
# These are tPix that are connected to other tPix via corridors.
# Mask out tPix that are not connected.


# Get patch ids of edge markers
tPix = tPix*markers



# Remove small objects
markers = remove_small_objects(markers, 25, connectivity=1)
# Relabel
markers, _ = ndi.label(markers, structure=np.ones((3,3)))

# Check if only one patch.
if _ == 1:
    raise Exception('There is only one patch. You may want to change the threshold for creating patches.')

# Get list of patch ids
pids = list(range(1,np.max(markers)+1))

# Get patch edges
fprint = np.ones([3,3])
mBin = np.copy(markers)
mBin = np.where(mBin > 0, 1, 0)

aa = mBin - binary_erosion(mBin, fprint)

# Mask markers with edges
markers = aa*markers

# Get edges that overlap with corridors (i.e. entry-exit points)
markers[lcc <=0 ] = 0


#%%
# Loop through patch edges and get cell counts
patchCountList = []
for i in pids:
    patchCountList.append(np.sum(markers==i))
# Get proportioanal allocation
allCounts = np.sum(patchCountList)

# Proportional allocation of points to patches
patchPList = []
for i in patchCountList:
    patchPList.append(i/allCounts*nPts)
# Round to integer
patchPList = np.rint(patchPList)

# Check if each patch has at least 1 points
# If not, raise an exception suggesting user increase
# number of points
if np.min(patchPList) < 1:
    raise Exception('There are not enough points in each patch. You may want to increase the number of points.')

# Create patch pair array
ppa = np.array(list(itertools.combinations(pids, 2)))

# Reorder pairs so that largest valued id is on the right
ppa[:,] = np.array([np.min(ppa,axis=1),np.max(ppa,axis=1)]).T

# Loop through patch pairs and calculate patch pair
# corridor networks

    

    # NOTE: Need to determine how to catch cases
    # where there are no pairwise connections between
    # patches. Might have to do an initial pass with 
    # SPSP.
    print("Calculating corridors between patch " + str(pid1) + ' and patch ' + str(pid2))
    os.system("python lcc_memsafe.py 1")    

    
    
    # Get row, column indices where values are > 0
    rcIndices = np.argwhere(hvP > 0)
    
    # Get random subset of indices to index into row, column object.
    # Sort from smallest to largest.
    sSet = np.random.choice(rcIndices.shape[0], 50, replace=False)
    sSet.sort()
    
    # Select random row, column indices
    randIndices = rcIndices[sSet,:]
    
    # Convert row, column to coordinates
    xCs, yCs = rio.transform.xy(src.transform, randIndices[:,0], randIndices[:,1])

    # Convert to pandas df
    dfCoords = pd.DataFrame({'X': xCs, 'Y': yCs})
    
    # Convert to geodataframe
    ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk100_q75_sspts.shp"
    gdf = gpd.GeoDataFrame(dfCoords, geometry=gpd.points_from_xy(dfCoords.X, dfCoords.Y), crs=src.crs)
    gdf.to_file(ofile)
    


    

#%%
q = 0.75

# Read in original source points
xyf = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/sabah_100.shp'
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

# Read in crk raster
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk100_test6.tif") as src:
    profile = src.profile
    im = src.read(1)
    cinds = cf.cell_indices_from_coords(src, im, np.array(xy))

# Create patches from high values of dispersal kernels
hvP = np.where(im >= np.quantile(im[im>0],q), 1, 0)

# Find points that intersect with patches
ix = hvP[cinds[:,0],cinds[:,1]] > 0

# Subset pandas data frame, convert to geodataframe and write to file
xySS = xy.iloc[ix]

# Convert to geodataframe and write to file
gdf = gpd.GeoDataFrame(xySS, geometry=gpd.points_from_xy(xySS.x, xySS.y), crs="ESRI:102028")
ofile = 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/sabah_100_crkSS.shp'
gdf.to_file(ofile)

# Write patches to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk100_q75_patches2.tif"
hvP = np.expand_dims(hvP, axis=0)
cf.arrayToGeoTiff(hvP, ofile, profile)


#%% OLD
q = 0.75

# Read in crk raster
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk500_test2.tif") as src:
    #profile = src.profile
    im = src.read(1)

# Read in least cost corridors
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\lcc500_test10.tif") as src:
    profile = src.profile
    lcc = src.read(1)

s = ndi.generate_binary_structure(2,2)

# Patches from high values of dispersal kernels
atest = np.where(im >= np.quantile(im[im>0],q), 1, 0)
# Distance to nearest patch
distance = ndi.distance_transform_edt(atest)
# Create markers (connected components (8 neighbor))
markers, _ = ndi.label(atest, structure=s)
# Create mask using positive values of corridors
wMask = np.where(lcc > 0, True, False)

# Use watershed algo to 'flood' least cost corridors
# so that they're allocated to the closest core area
# Can use euclidean distance to nearest patch
# or the least cost corridor values.
# Use negative of both so that high value areas
# are flooded first
labels = watershed(-lcc, markers, mask=wMask)

# Use mask out patch/watershed overlap areas
pwo = np.copy(labels)
pwo[atest == 1] = 0

profile2 = profile
profile2.update({'dtype': 'int32'})

# Write to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\wshed_lccvals_500_test1.tif"
labels = np.expand_dims(labels, axis=0)
cf.arrayToGeoTiff(labels, ofile, profile2)

ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\patches8n_500_test1.tif"
markers = np.expand_dims(markers, axis=0)
cf.arrayToGeoTiff(markers, ofile, profile2)


#%%
mask = np.zeros(im.shape, dtype=bool)

mask[atest==1] = True
markers, _ = ndi.label(mask)


cds = peak_local_max(im, min_distance=50)    

mask = np.zeros(im.shape, dtype=bool)

mask[tuple(cds.T)] = True

markers, _ = ndi.label(mask)

# Read in least cost corridors
with rio.open(r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\sabah_200_lcc3.tif") as src:
    profile = src.profile
    lcc = src.read(1)

# Use watershed algo to 'flood' least cost corridors
# so that they're allocated to the closest core area
labels = watershed(-lcc, markers, mask=lcc)

# Write to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\sabah_200_lcc_labels2.tif"
labels = np.expand_dims(labels, axis=0)
cf.arrayToGeoTiff(labels, ofile, profile)

# Write peaks to file
ofile = r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\sabah_200_crk_peaks2.tif"
markers = np.expand_dims(markers, axis=0)
cf.arrayToGeoTiff(markers, ofile, profile)

# Write top 10 percent crk to file
im[im==0] = np.nan
qx = np.nanquantile(im, 0.85)
crkqx = np.zeros(im.shape)
crkqx[im >= qx] = 1
plt.imshow(crkqx)

# scipy focal median
a = np.array([[1, 2, 4, 6, 8],
              [5, 1, 1, 4, 4],
              [1, 1, 20, 7, 5],
              [9, 3, 0, 0, 3],
              [2, 3, 1, 5, 6]])
#im = img_as_float(a)
cds = peak_local_max(a, min_distance=1)


aones = np.copy(a)
aones[aones>0] = 1

k = np.array([[0,1,0],[1,1,1],[0,1,0]])
#a1 = ndimage.convolve(a, k, mode='constant', cval=0.0)
#a2 = ndimage.convolve(aones, k, mode='constant', cval=0.0)

a1 = signal.fftconvolve(a, k,mode='same')
a2 = signal.fftconvolve(aones, k, mode='same')

a3 = a1/a2
plt.imshow(a3)

a, b = 100, 100
n = 201
rad = 101

y,x = np.ogrid[-a:n-a, -b:n-b]
mask = x*x + y*y <= rad*rad

k = np.zeros((n, n))
k[mask] = 1
plt.imshow(k)


a1 = ndimage.convolve(r, k, mode='constant', cval=0.0)

from scipy import signal

with rio.open('C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/sabah_200_crk_gauss.tif') as src:
    profile = src.profile
    r = src.read(1)
    
a1 = signal.fftconvolve(r,k, mode='same')

rones = np.copy(r)
rones[rones>0] = 1
a2 = signal.fftconvolve(rones,k, mode='same')

a3 = a1/a2

rtpi = r-a3
rtpi[r==0] = 0
rtpi = np.expand_dims(rtpi, axis=0)

cf.arrayToGeoTiff(rtpi, 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/sabah_200_crk_tpi4.tif', profile)

fMinfMax = np.expand_dims(fMinfMax, axis=0)
cf.arrayToGeoTiff(fMinfMax, 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/region_group.tif', profile)

tPix = np.expand_dims(tPix, axis=0)
cf.arrayToGeoTiff(tPix, 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/tPix.tif', profile)


cf.arrayToGeoTiff(np.expand_dims(labels, axis=0), 'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/edge_cell_nbs_test4.tif', profile)
