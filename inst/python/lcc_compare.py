# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:23:19 2024
Script to calculate lcc sums for each raster 
and graph as a barplot. If the first raster is a baseline, compare
subsequent rasters by calculating % of baseline.
@author: pj276
"""

#%%
# IMPORTS
import sys
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import cola_functions as cf

#%%
# INPUTS
# Resistance file for baseline scenario (to get no data mask)
rfBase = sys.argv[1]

# Comma separated string of raster names (no spaces between commas)
# Assumes files have nodata values of -9999 and have valid projections
csn = sys.argv[2]

# Output figure 1 name (png)
# This is barchart of the sum of lcc values for each scenario
ofig1 = sys.argv [3]

# Output figure 2 name (png)
# This is the barchart comparing lcc of the baseline scenario
# against all the others. Values are in % of baseline.
ofig2 = sys.argv[4]

# Output folder for writing rasters to file
# These are tiffs created by subtracting the baseline scenario
# from each of the other scenarios. Files are named using this
# pattern: s1_comp.tif, s2_comp.tif, etc.
odir = sys.argv[5]

# Set style
plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')

#%%
# Split string to list
nlist = csn.strip().split(',')

# Empty list to hold lcc sums
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
csum = pd.DataFrame({'lccsum': np.array(csumlist).transpose(), 'Scenario': ['S' + str(f) for f,n in enumerate(nlist)]})

# Barplot
ax = csum.plot.bar(x='Scenario', y='lccsum', rot=0, color="darkblue", title="Corridor Movement Potential", fontsize=16)
ax.title.set_size(16)
ax.xaxis.label.set_size(16)
ax.set_ylabel('Kernel Sum', fontsize=16)
ax.get_legend().remove()
fig = ax.get_figure()
plt.gcf().set_size_inches(6, 5)
plt.tight_layout()
fig.savefig(ofig1, dpi=300)

#%%
# Relative values
# Get first row of dataframe
basesum = csum['lccsum'][0]

# Drop first row from original dataframe
csum = csum.drop([0])

# Calculate basesum/scenario ratio
csum['lcccomp'] = csum['lccsum']/basesum*100

# Barplot
ax = csum.plot.bar(x='Scenario', y='lcccomp', rot=0, color="darkblue", title="Corridor Movement Scenario Comparison", fontsize=16)
ax.title.set_size(16)
ax.xaxis.label.set_size(16)
ax.set_ylabel('% of Baseline', fontsize=16)
#ax.text(0, 0, "{:.2f}".format(99.921600))
for p in ax.patches:
    #ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    ax.annotate("{:.2f}".format(p.get_height()), ((p.get_x() + p.get_width()/2.75), p.get_height() * 1.005))
ax.get_legend().remove()
fig = ax.get_figure()
plt.gcf().set_size_inches(6, 5)
plt.tight_layout()
fig.savefig(ofig2, dpi=300)

#%%
# Compare scenarios to baseline

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
    
    # Add a dimension to the array (the rasterio profile expects
    # a dimension corresponding to number of bands)
    comp = np.expand_dims(comp, axis=0)
    
    # Write to file
    otiff = Path(odir) / ('s' + str(i) + '_lcc_comp.tif' )
    cf.arrayToGeoTiff(comp, otiff, profile)
