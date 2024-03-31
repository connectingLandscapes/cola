# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:08:05 2024
Script to create source points spatially 
distributed according to movement resistance.
Movement resistance is resampled to a coarser
resolution based on expected population density
for the species for the area.

The basic idea is to pre-calculate all possible
pairwise cost distances on a population grid
(with a resolution similar to that of the expected
 population density). E.g. for a species with a
population density of 4 individuals per 100km2,
the population grid resolution would be 5x5km
(1 individual every 25km2). Cost distances are, however,
calculated using the resolution of the resistance layer.
We save the cost distance matrix to an hdf5 file
which can be accessed by CDPOP as needed.

We use the points and the cdmat info to automatically
create an XY file with occupied and unoccupied cells
plus additional infomartion (e.g. sex)

@author: pj276
"""

#%%
# IMPORTS
import sys
import rasterio as rio
from rasterio.crs import CRS
import cola_functions as cf
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
import time
import subprocess
import networkit as nk
import tables as tb

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()

    # INPUTS    
    # Path to suitability grid
    sg = sys.argv[1]
    
    # Path to resistance grid
    rg = sys.argv[2]

    # Output resampled suitability grid
    sgResam = sys.argv[3]
    
    # Output occupied points file
    opoints = sys.argv[4]
    
    # Output points file all no data points of resampled grid
    opointsAll = sys.argv[5]
    
    # Output CDMAT
    ocdmat = sys.argv[6]

    # Output xy file
    oxyfile = sys.argv[7]    

    # Suitability grid min value
    sMin = sys.argv[8]

    # Suitability grid max value
    sMax = sys.argv[9]
    
    # Raster cell size of CDPOP population grid
    cpgRes = sys.argv[10]
    
    # Number of source points to create
    numSourcePoints = sys.argv[11]
    
    # Distance threshold
    dThreshold = sys.argv[12]

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[13] # Default None

    #%%
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

    # Convert min, max and nsources to float or integer
    try:
        if cf.is_float(sMin):
            sMin = float(sMin)
        elif sMin.isdigit():
            sMin = int(sMin)
        else:
            float(sMin)
            int(sMin)
    except ValueError:
        print('Min value must be either a float or integer')
    try:
        if cf.is_float(sMax):
            sMax = float(sMax)
        elif sMax.isdigit():
            sMax = int(sMax)
        else:
            float(sMax)
            int(sMax)
    except ValueError:
        print('Max value must be either a float or integer')   
    try:
        numSourcePoints = int(numSourcePoints)
    except ValueError:
        print('Num source points must be integer.')
    
    #%%
    # Read suitability grid to array
    # Assign crs if provided by user
    if upCRS != 'None':
        with rio.open(sg, 'r+') as src:
            # Assign projection
            src.crs = CRS.from_string(upCRS)
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype in ['byte','int8','uint8','intp','uintp']:
                r = r.astype('float32')
                profile['dtype'] = 'float32'
    # Otherwise, read as usual
    else:            
        with rio.open(sg) as src:
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype in ['byte','int8','uint8','intp','uintp']:
                r = r.astype('float32')
                profile['dtype'] = 'float32'

    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)

    # Write to file so that gdal can resample
    # using the correct nodata value
    tmpSg = Path(sg).with_stem(Path(sg).stem + '_tmp_prof')
    # Add extra dim for writing
    r = np.expand_dims(r, axis=0)
    cf.arrayToGeoTiff(r, tmpSg, profile)
    
    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')
   
    # Resample suitability grid to coarser resolution
    args = args = ['gdalwarp', tmpSg, sgResam, '-ot', 'float32', '-r', 'average', '-srcnodata', '-9999', '-tr', cpgRes, cpgRes, '-co', 'compress=LZW']
    result = subprocess.call(args)
    
    with rio.open(sgResam) as src:
        # Get profile
        profile = src.profile
        # Read to array
        r = src.read(1)

    # Check no data value
    r, profile = cf.checkNoData(r, profile)

    #%%
    # Get indices where resampled raster is not equal to -9999
    allNDIndices = np.argwhere(r != -9999)

    # Convert row, column to coordinates
    allxCs, allyCs = rio.transform.xy(src.transform, allNDIndices[:,0], allNDIndices[:,1])
    
    # Convert occupied points to pandas df
    dfCoordsAll = pd.DataFrame({'X': allxCs, 'Y': allyCs})
    
    # Write to shapefile or csv
    if Path(opointsAll).suffix == '.csv':
        dfCoordsAll.to_csv(opointsAll, index=False)
    elif Path(opointsAll).suffix == '.xy':
        dfCoordsAll.to_csv(opointsAll, index=False)
    elif Path(opointsAll).suffix == '.shp':
        # Convert to geodataframe
        gdf = gpd.GeoDataFrame(dfCoordsAll, geometry=gpd.points_from_xy(dfCoordsAll.X, dfCoordsAll.Y), crs=src.crs)
        gdf.to_file(opointsAll)
    else:
        raise Exception('Output source point file should have a .csv, .xy, or .shp extension.')
    
    #%%
    # Convert no data value to nan
    r[r==-9999] = np.nan
    
    # Rescale to 0-1
    rRsc = (r - sMin) / (sMax - sMin)
    
    #%%
    # Subtract random array from rescaled suitability array
    rRsc = rRsc - np.random.rand(r.shape[0],r.shape[1])

    # Get row, column indices where values are > 0
    rcIndices = np.argwhere(rRsc > 0)

    # Throw an error if too many source points
    # have been requested (i.e. if there are not enough
    # habitat cells > 0 to satisfy the request)
    if numSourcePoints > rcIndices.shape[0]:
        raise Exception('There are not enough habitat cells to satisfy the requested number of source points. Try requesting fewer source points.')

    # Get random subset of indices to index into row, column object.
    # Sort from smallest to largest.
    sSet = np.random.choice(rcIndices.shape[0], numSourcePoints, replace=False)
    sSet.sort()

    # Select random row, column indices
    randIndices = rcIndices[sSet,:]
    
    # Convert row, column to coordinates
    xCs, yCs = rio.transform.xy(src.transform, randIndices[:,0], randIndices[:,1])
    
    # Convert occupied points to pandas df
    dfCoordsOcc = pd.DataFrame({'X': xCs, 'Y': yCs})
    
    # Write to shapefile or csv
    if Path(opoints).suffix == '.csv':
        dfCoordsOcc.to_csv(opoints, index=False)
    elif Path(opoints).suffix == '.xy':
        dfCoordsOcc.to_csv(opoints, index=False)
    elif Path(opoints).suffix == '.shp':
        # Convert to geodataframe
        gdf = gpd.GeoDataFrame(dfCoordsOcc, geometry=gpd.points_from_xy(dfCoordsOcc.X, dfCoordsOcc.Y), crs=src.crs)
        gdf.to_file(opoints)
    else:
        raise Exception('Output source point file should have a .csv, .xy, or .shp extension.')

    #%%
    # Read grid to float32 array, dealing with projection as needed
    r, profile = cf.read2flt32array(upCRS, rg)
        
    # Check no data value
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

    # Get sources as networkit node ids
    sources = cf.sources2nodeids(rg, r, dfCoordsAll, idmap)

    # Create empty hdf file to hold distance matrix
    shape = (len(sources), len(sources))
    atom = tb.Float32Atom()
    filters = tb.Filters(complib='blosc2:lz4', shuffle=True, complevel=5, fletcher32=False)
    h5f = tb.open_file(ocdmat, 'w')
    h5f.create_carray(h5f.root, 'dset', atom, shape,
                           filters=filters)
    h5f.close()
        
    # Calculate number of batches
    nBatch = np.ceil(len(sources)/500)
    
    sBatches = np.array_split(sources, nBatch)
    
    # Initialize list to provide min and max 
    # row indices for hdf file. These are sequential
    # while sources are not necessarily.
    sAcc = [0,0]
    # Initialize nan counter
    nanCount = 0
    for s in sBatches:
        if sAcc[1] == 0:
            sAcc[1] = len(s)
        else:
            sAcc[0] = sAcc[1]
            sAcc[1] = sAcc[0] + len(s)
        print('Working on batch:')
        print(sAcc)
        # Calculate shortest path distance from each source to
        # every other source
        spspDist = nk.distance.SPSP(nkG, s, sources)
        spspDist.run()
        ccArr = spspDist.getDistances(asarray=True)
        if dThreshold > 0:
            ccArr[ccArr > dThreshold] = np.nan
        ccArr[ccArr == np.finfo(np.float64).max] = np.nan
        ccArr = ccArr.astype('float32')
        # Insert in hdf file
        with tb.open_file(ocdmat, 'a') as h5f:
            h5f.root.dset[sAcc[0]:sAcc[1],:] = ccArr
        nanCount += np.sum(np.isnan(ccArr))    
        
    if nanCount == len(sources)**2:
        print('Cost distance matrix is empty. Check input points and grids.')
    #!!! Need to check how source points index into the
    # array. Need to check how/where CDPOP indexes the array.
    
    toc = time.perf_counter()
    print(f"Calculating cdmat took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
