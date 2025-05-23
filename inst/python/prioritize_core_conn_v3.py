# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 20:19:45 2025

@author: pj276
"""

#%%
# IMPORTS
import sys
from scipy import ndimage as ndi
import numpy as np
import rasterio as rio
from skimage.segmentation import watershed
from skimage.morphology import binary_erosion, binary_dilation, dilation
from skimage.measure import regionprops
import cola_functions as cf
import pandas as pd
from rasterio.features import shapes
import geopandas as gpd
import networkit as nk
import time
import itertools

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()

    # INPUTS
    # Original cost surface
    ocsFile = sys.argv[1]
    
    # Path to resistant kernel file (output of crk function)
    crkFile = sys.argv[2]
    
    # Path to least cost corridor file (output of lcc function)
    lccFile = sys.argv[3]
    
    # Output masked cost surface tiff (could be a temp file)
    maskedcsname = sys.argv[4]
    
    # Output corridor point shapefile name
    ocorrpoints = sys.argv[5]
    
    # Output corridor polygon shapefile name
    ocorrpoly = sys.argv[6]
 
    # Output high value patches shapefile name
    hvpatches = sys.argv[7]
    
    # Output high value patches tiff name
    opatchestiff = sys.argv[8]
    
    # Output sum of corridors tiff name
    osumcorr = sys.argv[9]
    
    # Quantile threshold for identifying high value patches from resistant kernel surface
    # Should be between 0 and 1.
    q = sys.argv[10] # 0.5
    q = float(q)
    
    # Corridor width parameter. As a percentage. E.g. a value of 10
    # would take the minimum cost distance + 10% and use that to determine
    # corridor width. If this is set low, the corridors from the prioritization
    # will be narrower than those in the FLCC layer. If this is set ver high,
    # the algorithm might grab additional cells nearby that are part of other
    # corridors, leading to overestimates of individual corridors. Setting it
    # to ~10-50 should yield the most interpretable results.
    corrWidthParam = sys.argv[11] # 

    # Convert tolerance to float or integer
    try:
        if cf.is_floatpy3(corrWidthParam):
            corrWidthParam = float(corrWidthParam)
        elif corrWidthParam.isdigit():
            corrWidthParam = int(corrWidthParam)
        else:
            float(corrWidthParam)
            int(corrWidthParam)
    except ValueError:
        print('Tolerance value must be either a float or integer')

    #%%    
    # Footprint to get patch edge pixels
    fprint = np.ones([3,3]).astype('uint8')

    # Read in crk raster
    with rio.open(crkFile) as src:
        profile = src.profile
        im = src.read(1)
    
    # Read in least cost corridors
    with rio.open(lccFile) as src:
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
    # Dilate markers so that corridor edge cells can be labeled
    markersDilate = dilation(markers, fprint)
    
    # Get patch attribute
    rprops = regionprops(markers, intensity_image=im)
    pLabels = [i.label for i in rprops]
    pAreas = [i.area for i in rprops]
    pMax = [i.intensity_max for i in rprops]
    pMean = [i.intensity_mean for i in rprops]
    ppCount = [i.num_pixels for i in rprops] 
    ppSum = [i[0]*i[1] for i in zip(pMean,ppCount)]
    
    # Check if only one patch.
    if npatch1 == 1:
        print('There is only one patch. You may want to change the threshold for creating patches.')
        sys.exit(0)
        #raise Exception('There is only one patch. You may want to change the threshold for creating patches.')
    
    # Convert patches to shapefile
    # Mask out zeros
    pmask = markers > 0
    # Get src for georeference
    with rio.open(ocsFile) as src:
        # Set up dictionary to hold geometries
        results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v) 
        in enumerate(
            shapes(markers, mask=pmask, transform=src.transform)))
        geoms = list(results)
        # Convert polygons to geodataframe
        patchShape = gpd.GeoDataFrame.from_features(geoms)
    # Use groupby function to convert to multipolygon    
    patchShape = cf.groupby_multipoly(patchShape, by='raster_val')
    # Append crk sums to gdf
    ddf = pd.DataFrame({'crksums': ppSum, 'ind': patchShape.index})
    ddf = ddf.set_index('ind')
    patchShape = pd.concat([patchShape, ddf], axis=1)

    # Get patch edges
    #hvPedges = hvP - binary_erosion(hvP, fprint)
    # Get cells where corridors overlap with patch edges
    #tPix = np.where((hvPedges > 0) & (lcc > 0), 1, 0)
    
    # Get corridor edges
    cedges = binary_dilation(hvP, fprint) - hvP
    # Get cells where corridors overlap with patch edges
    tPix = np.where((cedges > 0) & (lcc > 0), 1, 0)

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
    #wMask = np.where(cpecells > 0, True, wMask)
   
    # Read in original cost surface
    with rio.open(ocsFile) as src:
        profilecs = src.profile
        ocs = src.read(1)

    # Use watershed algo to 'flood' least cost corridors
    # so that they're allocated to the closest edge cells
    # Can use euclidean distance to nearest patch
    # or the least cost corridor values or the resistance surface
    # Use negative so that high value areas
    # are flooded first unless using resistance surface in which case
    # low resistance areas are flooded first.
    #labels = watershed(-lcc, cpecells, mask=wMask, connectivity=2)
    labels = watershed(ocs, cpecells, mask=wMask, connectivity=2)
    
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
        #msk = np.where(markers==i, 1, 0)
        msk = np.where(markersDilate==i, 1, 0)
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

    # Check if no connecting corridors
    if ecc.shape[0] == 0:
        print('There are no corridors connecting your patches.')
        sys.exit(0)
        #raise Exception('There are no corridors connecting your patches.')

    #%%
    # Get cell size from profile of original cost surface
    cSize = profilecs['transform'][0]
    
    # Convert non-corridor  cells to no data (this sets patch cells to nodata
    # as well so that corridors don't get routed through patches).
    # Corridors running through patches would get partially masked later
    # with patches which leads to funny results.
    # However, include patch edge pixels for now. Might remove these in 
    # a later version as sometimes corridors may run through the edge
    # of a patch and come out the other side, leading to disjunct corridors.
    ocs = wMask*ocs
    ocs[ocs == 0] = -9999
    ocs[np.isnan(ocs)] = -9999
    
    # Write to file
    ocs = np.expand_dims(ocs, axis=0)
    cf.arrayToGeoTiff(ocs, maskedcsname, profilecs)
    
    # Remove extra dimension for further processing
    ocs = np.squeeze(ocs)
    
    # Create edges, nodeids, and mapping between old and new nodeids
    # from resistance grid
    edges, nodeids, idmap = cf.generate_edges(ocs, cSize)
    print('generated edges')

    # Create graph (nodes only)
    nkG = nk.Graph(len(idmap), weighted=True)
    # Add edges to graph
    for i in edges:
        nkG.addEdge(i[0], i[1], w=i[2], addMissing=False, checkMultiEdge=False)
    print('created graph')
    print(nk.overview(nkG))
    
    # !!! Write cpecells and labels to file (temporary code) !!
    cpecells2 = np.expand_dims(cpecells, axis=0)
    cf.arrayToGeoTiff(cpecells2, 'C:/Users/pj276/Downloads/prioritization5/edge_cell_groups23.tif', profilecs)
    markers2 = np.expand_dims(markers, axis=0)
    cf.arrayToGeoTiff(markers2, 'C:/Users/pj276/Downloads/prioritization5/markers_groups23.tif', profilecs)
    labels2 = np.expand_dims(labels, axis=0)
    cf.arrayToGeoTiff(labels2, 'C:/Users/pj276/Downloads/prioritization5/labels_groups23.tif', profilecs)
    
    
    #%%
    # Loop through edge pairs and calculate corridors
    # plus stats
    # Empty array to hold running sum of binary corridors
    rsumCorr = np.zeros((ocs.shape[0],ocs.shape[1]), dtype='int32')
    
    # Empty list to hold stats
    attList = []
    # Empty list to hold polygons
    polyList = []
    for pid, pp in enumerate(ecc):
        # Get first group of cells
        ex1 = np.where(cpecells==pp[0], 1, 0)
        # Get row, column indices where values are > 0
        ex1Inds = np.argwhere(ex1 > 0)
        # Convert row, column to coordinates
        with rio.open(ocsFile) as src:
            xCs1, yCs1 = rio.transform.xy(src.transform, ex1Inds[:,0], ex1Inds[:,1])
    
        # Get second group of cells
        ex2 = np.where(cpecells==pp[1], 1, 0)
        # Get row, column indices where values are > 0
        ex2Inds = np.argwhere(ex2 > 0)
        # Convert row, column to coordinates
        with rio.open(ocsFile) as src:
            xCs2, yCs2 = rio.transform.xy(src.transform, ex2Inds[:,0], ex2Inds[:,1])
    
        # Convert to pandas df
        dfCoords1 = pd.DataFrame({'X': xCs1, 'Y': yCs1})
        # Convert to pandas df
        dfCoords2 = pd.DataFrame({'X': xCs2, 'Y': yCs2})
    
        # Get 1st set of sources as networkit node ids
        sources1 = cf.sources2nodeids(maskedcsname, ocs, dfCoords1, idmap)
        # Get 2nd set of sources as networkit node ids
        sources2 = cf.sources2nodeids(maskedcsname, ocs, dfCoords2, idmap)
        
        # Get all pairwise combinations within groups of 1st set of nodes
        if len(sources1) > 1: # no need to set edge weight if only 1 node
            s1pairs = np.array(list(itertools.combinations(sources1, 2)))
            # Make sure node ids are in ascending order
            s1pairs = np.array([np.min(s1pairs,axis=1),np.max(s1pairs,axis=1)]).T
            # Set edge weights within groups of cells to zero
            for n in s1pairs:
                if nkG.hasEdge(n[0],n[1]):
                    nkG.setWeight(n[0],n[1],0)

        # Get all pairwise combinations within groups of 2nd set of nodes
        if len(sources2) > 2: # no need to set edge weight if only 1 node
            s2pairs = np.array(list(itertools.combinations(sources2, 2)))
            # Make sure node ids are in ascending order
            s2pairs = np.array([np.min(s2pairs,axis=1),np.max(s2pairs,axis=1)]).T
            # Set edge weights within groups of cells to zero
            for n in s2pairs:
                if nkG.hasEdge(n[0],n[1]):
                    nkG.setWeight(n[0],n[1],0)
        
        # Empty list to hold pair of distance rasters
        pdrs = []
        for j in [sources1[0], sources2[0]]:
            # Calculate cost distance from a node in each group of cells to every
            # cell in the landscape. We only need one node from each group of
            # cells because the groups of cells have zero internal resistance
            spspDist = nk.distance.SPSP(nkG, [j])
            spspDist.run()
            ccArr = spspDist.getDistances(asarray=True)
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
            ccArr[ccArr == np.finfo(np.float64).max] = np.nan
            # Create 1D zeros array to hold distance values
            dArr = np.zeros([ocs.shape[0]*ocs.shape[1]], 'float32')
            # Update with distance values
            dArr[nodeids] = ccArr
            del ccArr
            # Reshape into rectangular array
            pdrs.append(dArr.reshape(ocs.shape[0],ocs.shape[1]))
        
        # Distance between sets of points
        dsp = pdrs[0] + pdrs[1]
        # Threshold
        dsp = np.where(dsp < np.min(dsp[dsp>0])*(corrWidthParam/100 + 1), dsp, 0)
        
        # Set nan to 0
        dsp[np.isnan(dsp)] = 0
        
        # Mask out corridor where it overlaps with patches
        #dsp = np.where(np.logical_and(hvP > 1, dsp > 0), 0, dsp)
        dsp = np.where(hvP == 1, 0, dsp)
        
        # Segment corridors. Retain only largest one if there are multiple.
        # This essentially erases isolated portions of corridors that are
        # an artifact of the CMTC approach
        dspBin = np.where(dsp > 0, 1, 0)
        csegs, nc1 = ndi.label(dspBin, structure=fprint)
        
        # If there is more than one corridor segment
        if nc1 > 1:            
            # Get corridor attributes
            cprops = regionprops(csegs)
            cLabels = [i.label for i in cprops]
            cAreas = [i.area for i in cprops]
            keepcid = cLabels[np.argmax(cAreas)] # get label of largest corridor segment
            dsp = np.where(csegs==keepcid, dsp, 0)
        
        # Minimum cost distance value between sets of points
        mcd = np.min(dsp[dsp>0])
        meancd = np.mean(dsp[dsp>0])
        
        # Get corridor strength values overlapping with corridor
        tCorr = np.where(dsp > 0, lcc, 0)
        
        # Corridor props
        # Label discrete corridors with points (this still doesn't work very well)
        tCorrBinary = np.where(tCorr > 0, 1, 0)
        corrMarkers, ncorr1 = ndi.label(tCorrBinary, structure=fprint)
        cprops = regionprops(corrMarkers, intensity_image=tCorr)
        centrow, centcol = [i.centroid for i in cprops][np.argmax([i.num_pixels for i in cprops])]
        with rio.open(ocsFile) as src:
            centx, centy = rio.transform.xy(src.transform, int(np.round(centrow)), int(np.round(centcol)))
        
        # Get max corridor strength value
        tCorrMax = np.max(tCorr)
        tCorrMean = np.mean(tCorr)
        tCorrSum = np.sum(tCorr)
        
        # Convert corridor to 0-1
        tCorr = np.where(tCorr > 0, 1, 0)
        tCorr = tCorr.astype('uint8')
        cmask = tCorr == 1
        cmask = cmask.astype('uint8')
        
        # Convert to shapefile
        with rio.open(ocsFile) as src:
            results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) 
            in enumerate(
                shapes(tCorr, mask=cmask, transform=src.transform)))
            geoms = list(results)
            corrShape = gpd.GeoDataFrame.from_features(geoms)
            corrShape['id'] = pid
            # Use groupby function to convert to multipolygon    
            corrShape = cf.groupby_multipoly(corrShape, by='id')
            # Append to list
            polyList.append(corrShape)

        # Add to running sum
        rsumCorr += tCorr
        
        # Get index of first patch label
        pl1 = pLabels.index(pp[2])
        pl2 = pLabels.index(pp[3])
        # Use index values to retrieve patch attributes
        # mean corridor xcoord, mean corridor ycoord, patch id 1, patch id 2, edge id 1, edge id 2, area patch 1, area patch 2, max patch 1, max patch 2, mean patch 1, mean patch 2
        # Add to list
        attList.append(np.array([centx,centy,pp[2],pp[3],pp[0],pp[1],pAreas[pl1],pAreas[pl2],pMax[pl1],pMax[pl2],pMean[pl1],pMean[pl2],ppSum[pl1],ppSum[pl2],mcd,meancd,tCorrMax,tCorrMean,tCorrSum]))
    
    # Combine attributes into single array
    attArray = np.vstack(attList)
    # Shapefile column names
    column_names = ['xco','yco','pid1','pid2','eid1','eid2','parea1','parea2','pmax1','pmax2','pmean1','pmean2','psum1','psum2','mincost','meancost','maxstrn','meanstrn','sumstrn']
    # Convert to dataframe
    cAttDf = pd.DataFrame(attArray, columns=column_names)
    # Add corridor quality metric
    cAttDf['cp1'] = np.mean([cAttDf.psum1 * cAttDf.psum2]) * cAttDf.sumstrn * 1/cAttDf.meancost
    # Add patch quality metric
    #cAttDf['pp1'] = cAttDf.psum1 * cAttDf.psum2
 
    # Convert attributes to geodataframe of points
    with rio.open(ocsFile) as src:
        gdf = gpd.GeoDataFrame(cAttDf, geometry=gpd.points_from_xy(cAttDf.xco, cAttDf.yco), crs=src.crs)
    # Save to file
    gdf.to_file(ocorrpoints)

    # Concatenate corridor polygons to a single geodataframe
    bigPoly = pd.concat(polyList)
    # Add attributes to polygons
    bigPoly = pd.concat([bigPoly, cAttDf], axis=1)
    # Set crs
    bigPoly.set_crs(src.crs, inplace=True)
    # Write corridor polygons to shapefile
    bigPoly.to_file(ocorrpoly)
    
    # Write high value patches to shapefile
    # Set crs
    patchShape.set_crs(src.crs, inplace=True)
    patchShape.to_file(hvpatches)
    
    # Write high value patches to tiff
    cf.arrayToGeoTiff(np.expand_dims(markers, axis=0), opatchestiff, profilecs)
    
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    rsumCorr = np.expand_dims(rsumCorr, axis=0)   
    # Write sum of corridors to tiff
    cf.arrayToGeoTiff(rsumCorr, osumcorr, profilecs)
        
    toc = time.perf_counter()
    print(f"Prioritizing corridors and kernels took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
