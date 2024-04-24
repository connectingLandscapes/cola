# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:41:27 2024
Script to calculate crk sums for each raster (potentially within a polygon)
and graph as a barplot. If the first raster is a baseline, compare
subsequent rasters by calculating % of baseline.
@author: pj276
"""
#%%
# IMPORTS
import sys
import cola_functions as cf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
# INPUTS
# Comma separated string of raster names (no spaces between commas)
# Assumes files have nodata values of -9999 and have valid projections
csn = sys.argv[1]

# Output figure 1 name (png)
ofig1 = sys.argv [2]

# Output figure 2 name (png)
ofig2 = sys.argv[3]

#%%
# Split string to list
nlist = csn.strip().split(',')

# Empty list to hold crk sums
csumlist = []

for i in nlist:
    # Read grid to array
    # Assign crs if provided by user
    # Convert profile to single band if multi
    # Convert data and profile to float32
    r, profile = cf.read2flt32array("None", nlist[0])
    # Sum raster values
    rSum = np.sum(r[r != -9999])
    # Append to list
    csumlist.append(rSum)

#%%
# Absolute values
# Convert values to a dataframe
csum = pd.DataFrame({'crksum': np.array(csumlist).transpose(), 'Scenario': ['S' + str(f) for f,n in enumerate(nlist)]})

# Barplot
ax = csum.plot.bar(x='Scenario', y='crksum', rot=0, color="darkblue", title="Core Movement Potential", fontsize=16)
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
basesum = csum['crksum'][0]

# Drop first row from original dataframe
csum = csum.drop([0])

# Calculate basesum/scenario ratio
csum['crkcomp'] = basesum/csum['crksum']*100

# Barplot
ax = csum.plot.bar(x='Scenario', y='crkcomp', rot=0, color="darkblue", title="Core Movement Scenario Comparison", fontsize=16)
ax.title.set_size(16)
ax.xaxis.label.set_size(16)
ax.set_ylabel('% of Baseline', fontsize=16)
ax.get_legend().remove()
fig = ax.get_figure()
plt.gcf().set_size_inches(6, 5)
plt.tight_layout()
fig.savefig(ofig2, dpi=300)

