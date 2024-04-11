# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 11:13:37 2023
Calculate kernel density estimates on point locations. When points are
clumped, ISJ will tend to provide better bandwidth estimates. When points
are spread out, the average will likely perform best. CV is robust in
most cases but can take time to run on datasets with several thousand
points (e.g. > ~5000).
Warning: this algorithm doesn't check to see if points intersect with
template raster/grid.
@author: pj276
"""
#%%
# IMPORTS
import osgeo
import rasterio as rio
from rasterio.crs import CRS
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from KDEpy.bw_selection import improved_sheather_jones
from KDEpy.bw_selection import silvermans_rule
from KDEpy.bw_selection import scotts_rule
import geopandas as gpd
import pandas as pd
import numpy as np
import cola_functions as cf
import sys, time
from pathlib import Path

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()
    
    # INPUTS
    # Path to file holding xy coordinates
    xyf = sys.argv[1] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/sabah_test_spts_5000.shp'
    
    # x coordinate column name. If shapefile, use 'none'
    xcName = sys.argv[2]
    
    # y coordinate column name. If shapefile, use 'none'
    ycName = sys.argv[3]
    
    # Path to template grid
    tg = sys.argv[4] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example_pro.tif'

    # Output file path
    ofile = sys.argv[5] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_500pts_popDens_sklearn_isj_mean.tif'
    
    # Bandwidth method - one of 'isj', 'silvermans', 'scotts', 'average', 'cv', 'user'
    bwMethod = sys.argv[6]
    
    # User supplied bandwidth. Either 'none' (default) or a float or integer
    usbw = sys.argv[7]
    
    # Return count or density
    kdeType = sys.argv[8] # either 'count' or 'density'
    
    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[9] # Default None

    # Convert bandwidth to float or integer
    try:
        if usbw != 'none':
            if cf.is_float(usbw):
                usbw = float(usbw)
            elif usbw.isdigit():
                usbw = int(usbw)
            else:
                float(usbw)
                int(usbw)
            bw = usbw
    except ValueError:
        print('Bandwidth value must be either a float or integer')

    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array(upCRS, tg)
        
    # Check no data value and convert to 
    # -9999 if necessary
    r, profile = cf.checkNoData(r, profile)

    # Check cell size
    if profile['transform'][0] != np.abs(profile['transform'][4]):
        raise Exception('X and Y cell dimensions must be equal. Reformat grid so that cell dimensions are the same in X and Y directions.')
    # Get non nodata value cells
    rows, cols = np.where(r != -9999)
    # Get coordinates of valid cells
    # using grid as template
    with rio.open(tg) as src:
        xs,ys = rio.transform.xy(src.transform, rows, cols)
    # Make into array
    xflat = np.column_stack([xs,ys])

    #%%
    # Read xy file to dataframe and convert to array with two columns (x and y coordinates)
    if Path(xyf).suffix == '.csv':
        xy = pd.read_csv(xyf, index_col=False)
        xy = np.array(xy.loc[:,[xcName,ycName]])
    elif Path(xyf).suffix == '.shp':
        xy = gpd.read_file(xyf)
        xy = np.array(xy.get_coordinates())
    else:
        raise Exception('Source point file should have a .csv or .shp extension.')

    #%%
    # Bandwidth estimates to aid in determining optimal bandwidth for application
    # Improved Sheather Jones - ISJ
    bw1isj = improved_sheather_jones(xy[:, [0]]) # Get bandwidth for x coords
    bw2isj = improved_sheather_jones(xy[:, [1]]) # Get bandwidth for y coords
    # Silverman's rule
    bw1si = silvermans_rule(xy[:, [0]]) # Get bandwidth for x coords
    bw2si = silvermans_rule(xy[:, [1]]) # Get bandwidth for y coords
    # Scott's rule
    bw1sc = scotts_rule(xy[:, [0]]) # Get bandwidth for x coords
    bw2sc = scotts_rule(xy[:, [1]]) # Get bandwidth for y coords
    
    # Print results
    print('ISJ bandwidth estimate: ' + str(np.mean([bw1isj,bw2isj])))
    print('Silverman\'s bandwidth estimate: ' + str(np.mean([bw1si,bw2si])))
    print('Scott\'s bandwidth estimate: ' + str(np.mean([bw1sc,bw2sc])))
    print('Average bandwidth: ' + str(np.mean([bw1isj,bw2isj,bw1si,bw2si,bw1sc,bw2sc])))
    
    # Get bandwidth
    if bwMethod == 'isj':
        bw = np.mean([bw1isj,bw2isj])
    if bwMethod == 'silvermans':
        bw = np.mean([bw1si,bw2si])
    if bwMethod == 'scotts':
        bw = np.mean([bw1sc,bw2sc])
    if bwMethod == 'average':
        bw = np.mean([bw1isj,bw2isj,bw1si,bw2si,bw1sc,bw2sc])
    if bwMethod == 'cv':
        # Sklearn cross-validation for empirical bandwidth estimation
        # Divide the search space (determined from rule of thumb estimate) into 100 intervals
        # and run grid search 10 fold cv to identify a reasonable bandwidth
        grid = GridSearchCV(KernelDensity(),
                            {'bandwidth': np.linspace(np.min([bw1isj,bw2isj,bw1si,bw2si,bw1sc,bw2sc]), np.max([bw1isj,bw2isj,bw1si,bw2si,bw1sc,bw2sc]), 100)},
                            cv=10) # 10-fold cross-validation
        grid.fit(xy)
        print('Cross validation bandwidth estimate: ' + str(grid.best_estimator_.bandwidth))
        bw = grid.best_estimator_.bandwidth
    if (bwMethod == 'user') and (usbw != 'none'):
        bw = usbw

    #%%        
    # Create kernel density object
    print('Bandwidth used for kernel density estimate: ' + str(bw))
    kde = KernelDensity(bandwidth=bw)
    # Fit the kernel
    kde.fit(xy)
    # Predict on grid cells. Use np.exp to transform the output from sklearn kernel estimator.
    pDensGauss = np.exp(kde.score_samples(xflat))
    # Scale so the surface integrates to the population size.
    # In this case, cell values represent the expected count.
    # Basically, sum the raw density estimates, divide into the population
    # size (assumed to be the number of rows in point file), and multiply
    # that factor against the raw density estimates.
    # E.g. if you have 500m x 500m cells, this output gives you an estimate
    # of counts per 500m x 500m cell.
    pNorm = pDensGauss*(xy.shape[0]/np.sum(pDensGauss))
    
    # If densities are requested, calculate cell area in km2 and divide into
    # expected counts so that cell values represent count per unit area
    # E.g. if you have 500m x 500m cells, this output gives you an estimate
    # of counts per km2 in each cell
    if kdeType == 'density':
        cellArea = (profile['transform'][0]**2)/1000000
        pNorm = pNorm/cellArea

    #%%    
    # Write to file
    popDens = np.copy(r)
    popDens[rows,cols] = pNorm
    popDens = np.expand_dims(popDens, axis=0)
    cf.arrayToGeoTiff(popDens, ofile, profile)
    del popDens

    toc = time.perf_counter()
    print(f"Interpolating population density took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()