# -*- coding: utf-8 -*-
"""
Created on Tue May 20 20:24:14 2025
Script to implement cumulative resistant kernel mapping using Networkit.
If the memory threshold is set low enough, the script calculates
distances from source points to every cell in the map in batches
and writes to an hdf file on disk. It then reads distances from
the hdf file and calculates kernels in parallel
using joblib. Best performance will be with a solid state drive.
@author: pj276
"""

# Get raster bounds
with rio.open(rg) as src:
    minX = src.bounds.left
    maxX = src.bounds.right
    maxY = src.bounds.top
    minY = src.bounds.bottom
    rasterCRS = src.crs
 
# Create a fishnet
# From https://spatial-dev.guru/2022/05/22/create-fishnet-grid-using-geopandas-and-shapely/
# Creates from bottom left
x, y = (minX, minY)
geom_array = []
 
# Polygon Size
square_size = 5000*cSize
while y <= maxY:
    while x <= maxX:
        geom = geometry.Polygon([(x,y), (x, y+square_size), (x+square_size, y+square_size), (x+square_size, y), (x, y)])
        print([(x,y), (x, y+square_size), (x+square_size, y+square_size), (x+square_size, y), (x, y)])
        geom_array.append(geom)
        x += square_size
    x = minX
    y += square_size
 
fishnet = gpd.GeoDataFrame(geom_array, columns=['geometry']).set_crs(rasterCRS)
fishnet.to_file("C:/Users/pj276/Downloads/mainland_layers/fishnet_grid.shp")

#save output files as per shapefile features
for i in range(len(fishnet)):
    geom = []
    coord = shapely.geometry.mapping(fishnet)["features"][i]["geometry"]
    geom.append(coord)
    
    geomb = []
    coord = shapely.geometry.mapping(fishnet.iloc[i].geometry.buffer(50000, cap_style="square"))
    geomb.append(coord)
    
    # Unbuffered
    with rio.open(rg)as src:
        out_image, out_transform = rio.mask.mask(src,geom,crop=True)
        out_meta = src.meta
    out_meta.update({'driver':'GTiff',
                'height':out_image.shape[1],
                'width':out_image.shape[2],
                'transform':out_transform})
    # If there's data in the tile, process
    if np.sum(out_image > 0):
        with rio.open("C:/Users/pj276/Downloads/mainland_layers/testtile_" + str(i) + ".tif",'w',**out_meta) as dest:
           dest.write(out_image)
        # Buffered
        with rio.open(rg)as src:
            out_image2, out_transform2 = rio.mask.mask(src,geomb,crop=True)
            out_meta2 = src.meta
        out_meta2.update({'driver':'GTiff',
                    'height':out_image2.shape[1],
                    'width':out_image2.shape[2],
                    'transform':out_transform2})
        with rio.open("C:/Users/pj276/Downloads/mainland_layers/testtileB_" + str(i) + ".tif",'w',**out_meta2) as dest:
           dest.write(out_image2)

#%%
# IMPORTS
import sys
import osgeo
import cola_functions as cf
import networkit as nk
import rasterio as rio
from rasterio.crs import CRS
from rasterio import windows
from rasterio.windows import from_bounds
import pandas as pd
import geopandas as gpd
import numpy as np
from pathlib import Path
import time
import math
import tables as tb
from joblib import Parallel, delayed
from itertools import product
import tempfile


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

    # Set memory size for processing kernels
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
    
    # Check values between 0-1. If detected, this converts those values to
    # 1 and prints a message.
    r, profile = cf.checkRasterVals(r, profile)
    
    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat resistance grid so that cell dimensions are the same in X and Y directions.')

    # Get cell size from profile
    cSize = profile['transform'][0]
    
    # Get number of cells in the layer
    nvc = r.shape[0]*r.shape[1]
    #nvc = np.sum(r >= 1)
    
    # Get number of cells in the kernel
    nck = (math.pi*(dThreshold**2))/(cSize**2)
    
    def fdiff(lst):
        if len(lst) < 2:
            return []
        return [lst[i+1] - lst[i] for i in range(len(lst) - 1)]

    # Tile if the number of valid cells is > 20 million
    # and if the number of cells in the kernel is less
    # than 50% of the number of valid cells.
    # If the kernel is too large relative to the number
    # of valid cells, it doesn't make sense to tile
    if (nvc > 20e6) and (nck/nvc*100 < 50):
        tMult = 0
        acct = nvc
        while acct > 20e6:
            acct = acct/2
            tMult += 1
        nTiles = 2*tMult
    else:
        nTiles = 1
    
    if nTiles > 1:
        # Calculate tile overlap in number of cells using the dispersal
        # threshold
        overlap = np.ceil(dThreshold/cSize)
        int((r.shape[0]*r.shape[1]/nTiles)**0.5)
        # Get list of tile extents
        # Row offsets
        #rI = np.ceil(r.shape[0]/(nTiles/2))
        #rI = np.ceil(r.shape[0]/2)
        rI = int((r.shape[0]*r.shape[1]/nTiles)**0.5)
        acc = 0
        rspL = []
        while acc < r.shape[0]:
            rspL.append(int(acc))
            acc += rI
        # Get window heights
        wHts = fdiff(rspL + [r.shape[0]])
        wHts = wHts*int((nTiles/2)) # repeat list so there's a height for each window
        # we repeat the whole list for heights because the order of windows
        # moves from top left down, then right
        
        # Column offsets
        #cI = np.ceil(r.shape[1]/(nTiles/2))
        cI = int((r.shape[0]*r.shape[1]/nTiles)**0.5)
        acc = 0
        cspL = []
        while acc < r.shape[1]:
            cspL.append(int(acc))
            acc += cI
        # Get window widths
        wWts = fdiff(cspL + [r.shape[1]])
        wWts = list(np.repeat(wWts,(nTiles/2))) # repeat list elements so there's a width for each window
        # we repeat list elements for widths because the order of windows
        # moves from top left down, then right
        
        # This gives offsets from top left, down, then right
        offsets = product(cspL, rspL)
        big_window = windows.Window(col_off=0, row_off=0, width=r.shape[1], height=r.shape[0])
        # Window list for tiles with buffers
        wList = []
        # Window list for tiles without buffers
        wList2 = []
        for i, j in enumerate(offsets):
            col_off = j[0]
            row_off = j[1]
            # Get windows with buffers
            window = windows.Window(
                col_off=col_off - overlap,
                row_off=row_off - overlap,
                width=rI,# wWts[i] + overlap * 2,
                height=rI) #wHts[i] + overlap * 2)
            wList.append(window.intersection(big_window))
            # Get windows without buffers
            window2 = windows.Window(
                col_off,
                row_off,
                rI,#wWts[i],
                rI)#wHts[i])
            wList2.append(window2.intersection(big_window))
    else:
        wList = [big_window]

    # Loop through windows/tiles and calculate kernels
    # Empty list to hold output file names
    ofnList = []
    for wi, w in enumerate(wList):
        print("working on " + str(wi))
        
        #%%
        # Create edges, nodeids, and mapping between old and new nodeids
        # from resistance grid
        with rio.open(rg) as src:

            # Copy profile to new variable so it can be updated
            kwargs_ub = src.profile
            
            # Read in buffered window
            r_b = src.read(1, window=w).squeeze()

            # Get window specific transform
            wst = rio.windows.transform(w, src.transform)
            # Update profile
            kwargs_ub.update({
                'dtype': rio.float32,
                'height': window.height,
                'width': window.width,
                'transform': wst})
            # Write window to temp file
            wtemp_b = str(Path(ofile).parent / Path(tempfile.mktemp()).name) + '.tif'
            with rio.open(wtemp_b, 'w', **kwargs_ub) as dst:
                #dst.write_band(band_id, dArr.astype(rio.float32), window=window)
                dst.write_band(1, r_b)

            # Copy profile to new variable so it can be updated
            kwargs_b = src.profile
            
            # Read in unbuffered window
            r_ub = src.read(1, window=wList2[wi]).squeeze()

            # Get window specific transform
            wst = rio.windows.transform(wList2[wi], src.transform)
            # Update profile
            kwargs_b.update({
                'dtype': rio.float32,
                'height': window.height,
                'width': window.width,
                'transform': wst})
            # Write window to temp file
            wtemp_ub = str(Path(ofile).parent / Path(tempfile.mktemp()).name) + '.tif'
            with rio.open(wtemp_ub, 'w', **kwargs_b) as dst:
                #dst.write_band(band_id, dArr.astype(rio.float32), window=window)
                dst.write_band(1, r_ub)
        
        # Create edge list from buffered window
        edges, nodeids, idmap = cf.generate_edges(r_b, cSize)
        print('Generated edges')
    
        # Create graph (nodes only)
        nkG = nk.Graph(len(idmap), weighted=True)
        # Add edges to graph
        for e in edges:
            nkG.addEdge(e[0], e[1], w=e[2], addMissing=False, checkMultiEdge=False)
        print('Created graph')
        print('Number of nodes: ' + str(nkG.numberOfNodes()))
        print('Number of edges: ' + str(nkG.numberOfEdges()))
        del edges
    
        nodeidsLen = len(nodeids)
    
        #%%
        # Get sources as networkit node ids using non-buffered
        # windowed file as reference to use as a filter on points outside the
        # core window area
        with rio.open(wtemp_ub) as src:
            # Get values from window with no buffer
            wvals = np.array([x for x in src.sample(np.array(xy))])

        # Get sources as networkit node ids using buffered
        # windowed file as reference
        with rio.open(wtemp_b) as src:

            # Convert xy coordinates to rows and columns using
            # the resistance grid as a reference.
            cinds = cf.cell_indices_from_coords(src, r_b, np.array(xy))
            # Raise an exception if there are no source points in the landscape
            if np.sum(np.isnan(cinds)) >= 1:
                #raise Exception('Source points do not intersect resistance grid.')
                print('Source points do not intersect resistance grid. Moving on to next tile.')
                continue
            else:
                # Get number of points supplied by user
                nPts = len(cinds)
                # Check window values for points outside the windowed area
                if any(wvals < 0):
                    cinds = cinds[wvals[:,0] >= 0, :]
                # Check for points and remove any that are out of raster bounds
                # Apparently this can happen using cell_indices_from_coords
                # Outside outer dims
                if any(cinds[:,0] >= r_b.shape[0]):
                    cinds = cinds[cinds[:,0] < r_b.shape[0]]
                if any(cinds[:,1] >= r_b.shape[1]):
                    cinds = cinds[cinds[:,1] < r_b.shape[1]]
                # Outside inner dims
                if any(cinds[:,0] < 0):
                    cinds = cinds[cinds[:,0] >= 0]
                if any(cinds[:,1] < 0):
                    cinds = cinds[cinds[:,1] >= 0]
                # Check for individual points with no data
                checkND = r_b[np.array(cinds[:,0], dtype=int),np.array(cinds[:,1], dtype=int)]
                if -9999 in checkND:
                    cinds = cinds[checkND != -9999,:]

                # If there are no points after removing no data vals,
                # exit script.
                if len(cinds) == 0:
                    #raise Exception('No valid source points. This may be because none of your source points overlap with your resistance grid.')
                    print('No valid source points. This may be because none of your source points overlap with your resistance grid.')
                    continue
                # Check if any points were removed
                if nPts - len(cinds) > 0:
                    print('Ignoring ' + str(nPts - len(cinds)) + " source point(s) with no data values.", flush=True)
                # Get unique indices
                cinds = np.unique(cinds,axis=0)
                # Convert to list
                sources = [(i[0]*r_b.shape[1])+i[1] for i in cinds]
                # Convert to networkit nodeids
                sources = [idmap[n] for n in sources]
                print("Number of valid sources = " + str(len(sources)), flush=True)
    
        #%%
        # Use only 80% of user supplied threshold to be conservative
        gbThreshold = gbLim*0.8
    
        #%%
        # Estimate memory required for distance mapping by multiplying
        # the number of graph nodes by the number of sources
        gbCf = 0.000000000125 # bits to gb conversion factor
        memReq = nodeidsLen*len(sources)*64*gbCf
        
        # If it's more than the threshold amount, divide into batches
        # And map the distance between each source point and every cell in the landscape
        # --distance mapping--
        if memReq > gbThreshold:
            # Create hdf file for tile
            dahdfTile = str(Path(dahdf).parent / Path(dahdf).stem) + '_' + str(wi) + str(Path(dahdf).suffix)
            # Create hdf file to hold distance array
            shape = (len(sources), nodeidsLen)
            atom = tb.Float64Atom()
            filters = tb.Filters(complib='blosc2:lz4', shuffle=True, complevel=5, fletcher32=False)
            h5f = tb.open_file(dahdfTile, 'w')
            h5f.create_carray(h5f.root, 'dset', atom, shape,
                                   filters=filters)
            h5f.close()
            
            # Memory per processor
            memProc = gbThreshold/nThreads
    
            # Arrays per processor
            arraysProc = memProc/(nodeidsLen*64*gbCf)
    
            # Number of batches 
            nCBatches = int(np.floor(len(sources)/(arraysProc)))
            
            # Divide sources into batches
            sLength = np.arange(0,len(sources))
            #nCBatches = int(np.ceil(memReq/gbThreshold))
            
            sBatches = np.array_split(sLength, nCBatches) 
            print('Calculating cost distances in ' + str(len(sBatches)) + ' batches')
          
            # Loop over batches and calculate distance array
            for b, s in enumerate(sBatches):
                tic1 = time.perf_counter()
                print('Working on batch ' + str(b+1) + ' of ' + str(len(sBatches)))
                #print('Batch size is ' + str(len(s)))
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
                h5f = tb.open_file(dahdfTile, 'a')
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
            print('Calculating kernels in batches')
    
            # Memory per processor
            memProc = gbThreshold/nThreads
    
            # Arrays per processor
            arraysProc = memProc/(nodeidsLen*64*gbCf)
    
            # Number of batches 
            nBatches = int(np.floor(len(sources)/(arraysProc)))
    
            # Use 50 percent of the number of lines across all threads
            # as a way to allocate memory to the output list created
            # by joblib.
            jlibIncrement = int(np.ceil((len(sources)/nBatches*nThreads*0.5)/nThreads))
            
            # Main issue is that the nodeids are relative to the graph
            # while the indices need to be relative to the hdf file.
            # Dictionary to link nodeids to hdf indices
            nodehdf = dict(zip(sources, range(len(sources))))
            
            # Split nodeids into N batches
            sBatches = np.array_split(sources, nBatches)
            
            # Use joblib to calculate kernels in parallel and sum on the fly
            print("Summing kernels", flush=True)
            with Parallel(n_jobs=nThreads) as parallel:
                ccArr = np.zeros((1,nodeidsLen))
                n_iter = 0
                for j in range(jlibIncrement,int((np.ceil(nBatches/jlibIncrement)+1)*jlibIncrement),jlibIncrement):
                    j = min(j, nBatches)
                    results = parallel(delayed(cf.calcKernels)(k, nodeidsLen, nodehdf, dahdfTile, sBatches, dThreshold, tForm, tkv, kvol) for k in range(n_iter,j))
                    ccArr += sum(results)
                    n_iter = j

        # Otherwise, process as normal
        else:
            # Calculate shortest path distance from each source to every
            # cell in the landscape
            spspDist = nk.distance.SPSP(nkG, sources)
            spspDist.run()
            ccArr = spspDist.getDistances(asarray=True)
            # Networkit gives inaccessible nodes the max float 64 value.
            # Set these to nan.
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
    
            # Multiply if kernel volume not 1
            #if kvol > 1:
            #    ccArr = ccArr*(kvol * 3/(np.pi*dThreshold**2))
            #elif kvol < 1:
            #    raise Exception('Kernel volume should be >= 1.')
            
            # Sum kernels
            ccArr = np.sum(ccArr, axis=0)
        
        # Create 1D zeros array to hold kernel values
        dArr = np.zeros([r_b.shape[0]*r_b.shape[1]], 'float32')
        # Update with kernel values
        dArr[nodeids] = ccArr
        del ccArr
        # Reshape into rectangular array
        dArr = dArr.reshape(r_b.shape[0],r_b.shape[1])
    
        # Calculate quantiles of values > 0 and write to csv
        ##qDf = np.quantile(dArr[dArr>0],np.arange(0,1.01,0.01))
        ##qDf = pd.DataFrame({'q': np.arange(0,1.01,0.01), 'value': qDf})
        # Write quantiles to csv
        ##qDf.to_csv(ofile.replace(".tif","_quantiles.csv"),index=False)
        
        # Smooth (not implemented yet)
        #from scipy.ndimage import gaussian_filter
        #dArr = gaussian_filter(dArr, sigma=1, radius=2)
        
        #%% 
        # Add a dimension to the array (the rasterio profile expects
        # a dimension corresponding to number of bands)
        dArr = np.expand_dims(dArr, axis=0)
        
        # Write buffered crk to file
        ofileTile = str(Path(ofile).parent / Path(ofile).stem) + '_' + str(wi) + str(Path(ofile).suffix)
        cf.arrayToGeoTiff(dArr, ofileTile, kwargs_ub)    
        
        # If hdfs were created, delete them
        if Path(dahdf).is_file():
            Path(dahdf).unlink()

        # Delete temp files 
        if Path(wtemp_ub).is_file():
            Path(wtemp_ub).unlink()
        if Path(wtemp_b).is_file():
            Path(wtemp_b).unlink()
        
        # Add output kernel file name to list
        ofnList.append(ofileTile)
        toc = time.perf_counter()
        print(f"Calculating kernels took {toc - tic:0.4f} seconds", flush=True)

    # Check for files to mosaic
    if len(ofnList) > 1:
        # Create template for plugging in smaller arrays using resistance
        # raster profile.
        rsum = np.zeros((profile['height'],profile['width']), dtype='float32')
        # Loop through tiles and add to big raster
        for i in ofnList:
            # Get bounding box coordinates from tiles
            with rio.open(i) as src:
                rr = src.read(1)
                left = src.bounds.left
                right = src.bounds.right
                top = src.bounds.top
                bottom = src.bounds.bottom
                rr[rr < 0] = 0

            # Get smaller window row/col offsets relative to larger raster
            with rio.open(rg) as src:
                w1 = from_bounds(left, bottom, right, top, src.transform)

            # Add smaller window to existing values in larger array, then plug into larger array
            rsum[int(w1.row_off):int(w1.row_off+w1.height),int(w1.col_off):int(w1.col_off+w1.width)] = rr + rsum[int(w1.row_off):int(w1.row_off+w1.height),int(w1.col_off):int(w1.col_off+w1.width)]

        rsum = np.expand_dims(rsum, axis=0)
        cf.arrayToGeoTiff(rsum, ofile, profile)
        
# =============================================================================
#         # Create zeros image
#         crkCSum = np.zeros((r.shape[0],r.shape[1]))
#         # Mosaic by cumulative sum
#         for crkTile in ofnList:
#             # Read file subset into larger window
#             with rio.open(crkTile) as src:
#                 # Read in buffered window
#                 x = src.read(1, window=big_window).squeeze()
#                 # Add to zeros image
#                 crkCSum += x
#         # Save to file
#         # Add another dim for writing
#         crkCSum = np.expand_dims(crkCSum, axis=0)
#         # Write mosaiced crk to file
#         cf.arrayToGeoTiff(crkCSum, ofile, profile)
# =============================================================================

    else:
        # I think this branch represents what happens if there's no
        # tiling. Need to figure out how to structure this case.
        # Try just renaming the file
        Path(ofnList[0]).rename(Path(ofile))
        print('still working')

if __name__ == "__main__":
    main()

