# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 20:43:41 2023
Interpolate allelic diversity and heterozygosity
@author: pj276
"""

#%%
# Imports
import pandas as pd
import numpy as np
import rasterio as rio
import cola_functions as cf
from scipy.interpolate import RBFInterpolator
import sys, time
from pathlib import Path

def main() -> None:
    #%%
    # Start timer
    tic = time.perf_counter()
    
    # INPUTS
    # List of CDPOP grid.csv files containing population genetic structure
    xylist = sys.argv[1] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/sabah_test_spts_5000.shp'
    
    # Path to template grid
    tg = sys.argv[2] # 'C:/Users/pj276/Projects/UNICOR_arch/unicor/resist_sabah_example_pro.tif'
    
    # Interpolation method
    # Either 'multiquadric', 'thin_plate_spline', 'linear', or 'idw'
    iptMethod = sys.argv[3] 
    
    # Number of neighbors to use for interpolation or 'all'
    # to use all points in the dataset.
    # For idw, this argument has no effect
    rbfNbs = sys.argv[4] # 10

    # User provided CRS if using ascii or other file without projection info
    # Provide as epsg or esri string e.g. "ESRI:102028"
    upCRS = sys.argv[5] # Default None
    
    # Hard code output file names prefixes
    # Output alleles tif file name
    ofileAlleles = 'alleles'
    # Output heterozygosity tif file name
    ofileHo = 'heterozygosity'
    
    # Convert number of neighbors to integer
    try:
        if rbfNbs != 'all':
            if rbfNbs.isdigit():
                rbfNbs = int(rbfNbs)
            else:
                rbfNbs = int(rbfNbs)
        else:
            rbfNbs = None
    except ValueError:
        print('Number of neighbors must be integer')

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
    # Split string to list of files
    xylist = xylist.strip().split(',')
    
    # Loop through files and interpolate
    for xyf in xylist:
        # Read xy file to dataframe and convert to array 
        if Path(xyf).suffix == '.csv':
            xy = pd.read_csv(xyf, index_col=False)
        else:
            raise Exception('Source point file should have a .csv extension.')
        
        #%%
        # Get locations with data
        xy = xy.loc[xy['ID'] != 'OPEN',]
        # Remove unneeded columns
        xySubset = xy.iloc[:,9:]
        
        #%%
        # Get allele 1 count
        c1 = np.sum(xySubset == 1, axis=1)
        # Get allele 2 count
        c2 = np.sum(xySubset == 2, axis=1)
        # Get total allele count per location
        alleles = c1 + c2
        # Calculate observed heterozygosity
        Ho = (c1/2)/(c2+(c1/2))
        # Add allels and Ho back to dataframe
        fout = pd.concat([xy[['XCOORD','YCOORD']],alleles,Ho],axis=1)
        # Rename columns
        fout = fout.rename(columns={0: 'Alleles', 1: 'Ho'})
        
        # Check for all zeros/nan in calculated values
        # Both Ho and Alleles need to contain data, otherwise will skip
        if np.sum(np.isnan(fout.Ho)) < len(fout.Ho) and np.nanmax(fout.Alleles) > 0:
            
            #%%
            # Create interpolation model using genetic data
            # and predict for raster cell locations
            if iptMethod == 'multiquadric':
                iptAlleles = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Alleles, neighbors=rbfNbs, kernel='multiquadric', epsilon=1)(xflat)
                iptHo = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Ho, neighbors=rbfNbs, kernel='multiquadric', epsilon=1)(xflat)
                iptabb = 'mq'
            if iptMethod == 'thin_plate_spline':
                iptAlleles = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Alleles, neighbors=rbfNbs, kernel='thin_plate_spline')(xflat)
                iptHo = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Ho, neighbors=rbfNbs, kernel='thin_plate_spline')(xflat)
                iptabb = 'tps'
            if iptMethod == 'linear':
                iptAlleles = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Alleles, neighbors=rbfNbs, kernel='linear',degree=1)(xflat)
                iptHo = RBFInterpolator(fout[['XCOORD','YCOORD']], fout.Ho, neighbors=rbfNbs, kernel='linear',degree=1)(xflat)
                iptabb = 'lin'
            if iptMethod == 'idw':
                iptAlleles = cf.simple_idw(fout.XCOORD,fout.YCOORD,fout.Alleles,xflat[:,0],xflat[:,1])
                iptHo = cf.simple_idw(fout.XCOORD,fout.YCOORD,fout.Ho,xflat[:,0],xflat[:,1])
                iptabb = 'idw'
        
            #%%
            # For now, idw takes all neighbors into account so set output number
            # of neighbors to 'all' when naming files
            if iptMethod == 'idw':
                rbfNbs = 'all'
            # Write alleles to file
            rAlleles = np.copy(r)
            rAlleles[rows,cols] = iptAlleles
            rAlleles = np.expand_dims(rAlleles, axis=0)
            of1 = Path(xyf).parent / Path(ofileAlleles + '_' + iptabb + '_' + str(rbfNbs) + '_' + str((Path(xyf).stem)) + '.tif')
            print(of1)
            cf.arrayToGeoTiff(rAlleles, of1, profile)
            del rAlleles
            
            # Write Ho to file
            rHo = np.copy(r)
            rHo[rows,cols] = iptHo
            rHo = np.expand_dims(rHo, axis=0)
            # Write to file
            of2 = Path(xyf).parent / Path(ofileHo + '_' + iptabb + '_' + str(rbfNbs) + '_' + str((Path(xyf).stem)) + '.tif')
            cf.arrayToGeoTiff(rHo, of2, profile)
            del rHo

        else:
            print("In input file " + str(Path(xyf).name) + ", either heterozygosity or alleles columns are empty/contain no data. Nothing to interpolate.")

    toc = time.perf_counter()
    print(f"Interpolating population genetic structure took {toc - tic:0.4f} seconds", flush=True)

if __name__ == "__main__":
    main()
