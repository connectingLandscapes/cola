# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:41:27 2024
Script to calculate crk sums for each raster 
and graph as a barplot. If the first raster is a baseline, compare
subsequent rasters by calculating % of baseline.
@author: pj276
"""

#%%
# IMPORTS
import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import colormaps
import rasterio as rio
from rasterio import mask
import pandas as pd
import geopandas as gpd
import cola_functions as cf

#%%
# INPUTS
# Resistance file for baseline scenario (to get no data mask)
rfBase = sys.argv[1]

# Comma separated string of raster names (no spaces between commas)
# Assumes files have nodata values of -9999 and have valid projections
csn = sys.argv[2]

# Table output of absolute values (csv)
table1 = sys.argv[3]

# Table output of relative values (csv)
table2 = sys.argv[4]

# Output figure 1 name (png)
# This is barchart of the sum of crk values for each scenario
ofig1 = sys.argv [5]

# Output figure 2 name (png)
# This is the barchart comparing crk of the baseline scenario
# against all the others. Values are in % of baseline.
ofig2 = sys.argv[6]

# Output folder for writing rasters to file
# These are tiffs created by subtracting the baseline scenario
# from each of the other scenarios. Files are named using this
# pattern: s1_comp.tif, s2_comp.tif, etc.
odir = sys.argv[7]

# Shapefile for summarizing
# Needs to have an ID field for grouping polygons
# Use None to summarize over the entire raster extent
shpZones = sys.argv[8]

# ID field for grouping polygons
# Use None if shapefile arg is None
idField = sys.argv[9]

# Set style
plt.style.use('ggplot')

#%%
# Split string to list
nlist = csn.strip().split(',')

# %%
# Summarize over entire raster if no shapefile supplied
if shpZones == 'None':

    # Empty list to hold crk sums
    csumlist = []
    
    for i in nlist:
        # Read grid to array
        # Assign crs if provided by user
        # Convert profile to single band if multi
        # Convert data and profile to float32
        r, profile = cf.read2flt32array("None", i)
        # Sum raster values
        rSum = np.sum(r[r != -9999])
        # Append to list
        csumlist.append(rSum)
    
    #%%
    # Absolute values
    # Convert values to a dataframe
    csum = pd.DataFrame({'crksum': np.array(csumlist).transpose(), 'Scenario': ['S' + str(f) for f,n in enumerate(nlist)]})
    
    # Reorder columns
    csum = csum.loc[:, ['Scenario','crksum']]
    
    # Write table to file
    csum.to_csv(table1, index=False)
    
    # Barplot
    colors = ['#440154FF' if value == 'S0' else '#0072B5FF' for value in csum['Scenario']]
    ax = csum.plot.bar(x='Scenario', y='crksum', rot=0, color=colors, title="Core Movement Potential", fontsize=18)
    ax.title.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.set_ylabel('Kernel Sum', fontsize=18)
    ax.get_legend().remove()
    fig = ax.get_figure()
    plt.gcf().set_size_inches(6, 5)
    plt.tight_layout()
    fig.savefig(ofig1, dpi=300)
    
    #%%
    # Relative values
    # Get first row of dataframe
    basesum = csum['crksum'][0]
    
    # Drop first row from original dataframe
    csum = csum.drop([0])
    
    # Calculate basesum/scenario ratio
    csum['crkcomp'] = (csum['crksum']/basesum*100)-100
    
    # Drop crksum
    csum = csum.drop('crksum', axis=1)
    
    # Write table to file
    csum.to_csv(table2, index=False)
    
    # Barplot
    colors = ['#BC3C29FF' if value < 0 else '#0072B5FF' for value in csum['crkcomp']]
    ax = csum.plot.bar(x='Scenario', y='crkcomp', rot=0, color=colors, title="Core Movement Comparison", fontsize=16)
    ax.title.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.set_ylabel('% Change Relative to Baseline', fontsize=18)
    # Bar labels
    ax.bar_label(ax.containers[0], fmt=lambda x: f'{x:.2f}', fontsize=18, label_type='edge', padding=3)
    ax.margins(y=0.1)
    # Adjust margins for clarity
    ax.margins(y=0.1)
    ax.get_legend().remove()
    fig = ax.get_figure()
    plt.gcf().set_size_inches(6, 5)
    plt.tight_layout()
    fig.savefig(ofig2, dpi=300)

else:
    # Read in shapefile
    zones = gpd.read_file(shpZones)
    
    # Get unique IDs in grouping field
    uids = np.unique(zones[idField])
    
    # Empty list to hold scenario ids
    sidlist = []
    # Empty list to hold polygon ids
    pidlist = []
    # Empty list to hold crk sums
    csumlist = []

    # Split string to list
    nlist = csn.strip().split(',')

    # Loop through scenario rasters
    for s, i in enumerate(nlist):
        
        # Read grid to array
        # Assign crs if provided by user
        # Convert profile to single band if multi
        # Convert data and profile to float32
        r, profile = cf.read2flt32array("None", i)
        
        # Loop through unique polygon IDs
        for u in uids:

            # Subset shapefile to target polygon ID
            zss = zones.loc[zones[idField] == u, :]
            
            # Read in target raster and mask out values not in polygon
            # From https://gis.stackexchange.com/questions/151339/rasterize-a-shapefile-with-geopandas-or-fiona-python
            with rio.open(i) as src:
                # Template array
                tArr = src.read(1)
                 
                # Rasterized features (raster cells should have same value as id field)
                #burned = features.rasterize(shapes=shapes, fill=0, out_shape=tArr.shape, transform=src.profile['transform'])
                # Mask 
                zmask = mask.mask(src, zss.geometry, all_touched=True, invert=False, nodata=None, filled=True, crop=False, pad=False, pad_width=0.5, indexes=None)
                # Extract array
                zmask = zmask[0].squeeze()
                # Set nan values to -9999
                zmask[np.isnan(zmask)] = -9999
                
                # Check if there are any nodata values
                if np.sum(zmask > -9999) > 0:
                    # Append scenario id to list
                    sidlist.append(s)
                    # Append polygon id to list
                    pidlist.append(u)
                    # Sum raster values in valid pixels that overlap with the polygon
                    rSum = np.sum(zmask[zmask != -9999])
                    # Append to list
                    csumlist.append(rSum)
    # Convert lists to dataframe
    csum = pd.DataFrame({'crksum': np.array(csumlist).transpose(), 'Scenario': sidlist, 'PolyID': pidlist})
    csum.Scenario = 'S' + csum.Scenario.astype('str')

    #%%
    # Absolute values
    # Pivot for grouped barplot
    pivot_df = csum.pivot(index='Scenario',columns='PolyID',values='crksum')
    
    # Remove polygons with zero values across all scenarios
    pivot_df = pivot_df.loc[:,(pivot_df.sum(axis=0) != 0)]

    # Write to file
    pivot_df.to_csv(table1, index=True)

    # Barplot
    ax = pivot_df.plot.bar(rot=0, title="Core Movement Potential", fontsize=18)
    ax.title.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.set_ylabel('Kernel Sum', fontsize=18)
    #ax.get_legend().remove()
    ax.legend(loc='upper center', bbox_to_anchor=(1.05, 1.05),
          ncol=1, fancybox=True, shadow=True)
    fig = ax.get_figure()
    plt.gcf().set_size_inches(8, 5)
    plt.tight_layout()
    fig.savefig(ofig1, dpi=300)

    #%%
    # Relative values
    # Get first row of dataframe
    basesum = pivot_df.iloc[0,:]
    
    # Drop first row from original dataframe
    pivot_df = pivot_df.drop(['S0'])
    pivot_df = pivot_df.T
    pivot_df = (pivot_df.divide(basesum, axis=0)*100)-100
    pivot_df = pivot_df.T

    # Write to file
    pivot_df.to_csv(table2, index=True)
    
    # Barplot
    ax = pivot_df.plot.bar(rot=0, title="Core Movement Comparison", fontsize=18, color=colormaps['tab20'].colors)
    ax.title.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.set_ylabel('% Change Relative to Baseline', fontsize=18)
    # Bar labels
    ax.bar_label(ax.containers[0], fmt=lambda x: f'{x:.2f}', fontsize=18, label_type='edge')
    ax.margins(y=0.1)
    # Adjust margins for clarity
    ax.margins(y=0.1)
    #ax.get_legend().remove()
    ax.legend(loc='upper center', bbox_to_anchor=(1.05, 1.05),
          ncol=1, fancybox=True, shadow=True)
    fig = ax.get_figure()
    plt.gcf().set_size_inches(8, 5)
    plt.tight_layout()
    fig.savefig(ofig2, dpi=300)

#%%
# Spatially compare scenarios to baseline

# Pop first list item to new object and simultaneously remove from original list
bline = nlist.pop(0)

# Read baseline to object
bline, prof1 = cf.read2flt32array("None", bline)

# Use baseline resistance layer to get no data mask
ndmask, pro1 = cf.read2flt32array("None", rfBase)

# Loop through scenarios, subtract, and write to file
for i, j in enumerate(nlist):
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array("None", j)

    # Subtract scenario from baseline
    comp = r - bline
    
    # Apply no data mask
    comp[ndmask == -9999] = -9999
    comp[comp == 0] = -9999
    
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    comp = np.expand_dims(comp, axis=0)
    
    # Write to file
    otiff = Path(odir) / ('s' + str(i+1) + '_crk_comp.tif' )
    cf.arrayToGeoTiff(comp, otiff, profile)
