import osgeo
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
from rasterio.crs import CRS
from scipy import ndimage as ndi
from scipy import signal
from scipy import sparse
from scipy.interpolate import RBFInterpolator
from scipy.ndimage import gaussian_filter
from scipy.sparse import diags
from shapely.geometry import shape
from skimage.feature import peak_local_max
from skimage.filters.rank import maximum
from skimage.filters.rank import minimum, maximum
from skimage.measure import regionprops
from skimage.morphology import binary_erosion
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
import base64
import gc
import geopandas as gpd
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
import rasterio as rio
import subprocess
import sys
import sys, time
import tables
import tables as tb
import time
print ("WELCOME -- libraries loaded succesfully")
import os; 
print ("---- WD:")
print(os.getcwd())
import cola_functions as cf
print (" cola_functions loadedS")