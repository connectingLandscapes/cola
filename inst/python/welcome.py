import os; 
print ("---- WD:", os.getcwd())


#from scipy.ndimage import gaussian_filter
# importance of corridors using standard approach.
from contextlib import closing
from itertools import chain
from itertools import combinations
from KDEpy.bw_selection import improved_sheather_jones
from KDEpy.bw_selection import scotts_rule
from KDEpy.bw_selection import silvermans_rule
from numba import njit, prange
from numba import set_num_threads
from pathlib import Path
from scipy import ndimage as ndi
from scipy import signal
from scipy import sparse
from scipy.interpolate import RBFInterpolator
from scipy.ndimage import gaussian_filter
from scipy.sparse import diags
from skimage.feature import peak_local_max
from skimage.filters.rank import maximum
from skimage.filters.rank import minimum, maximum
from skimage.measure import regionprops

from skimage.morphology import binary_erosion
from skimage.morphology import remove_small_objects

from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity

import base64
import gc
import glob
import h5py
import itertools
import matplotlib.pyplot as plt
import networkit as nk
import networkx as nx
import numpy as np
import numpy as np, time
import os
import pandas as pd
import random, copy
import subprocess
import sys
import sys, time
import tables
import tables as tb
import time

import skimage


from shapely.geometry import shape

# print('    .... all good until here ')

## Issues here:
#print('    .... check point 2')

import osgeo
from osgeo import gdal
import geopandas as gpd
import rasterio as rio
from rasterio.crs import CRS

#print('    .... debug ')

try: 
    from skimage.segmentation import watershed
except: 
    print ("Can't: from skimage.segmentation import watershed")

import cola_functions as cf
#print (" cola_functions loaded!")
print ("WELCOME -- libraries loaded successfully")

