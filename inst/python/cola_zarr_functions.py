# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 05:12:44 2025
Functions for creating graphs and zarr files
for parallel processing least cost corridors
on hpc.
@author: pj276
"""
#%%
# Imports
import numpy as np
from numba import njit
import zarr
from zarr.codecs import BloscCodec, BloscShuffle

#%%
# Functions

# Function to generate edges based on pixel adjacency and distance
# between cells. Returns a generator. This may replace the original
# after testing.
@njit
def generate_edgesMod(pixels, iArray, cellSize, keyArray):
    """
    Creates a list of tuples representing edges between nodes (pixels)
    in a raster. Each tuple contains two node ids and an edge weight
    calculated as a function of cell size weighted by distance between cell
    centers, with cardinal cells given a weight of 1 and diagonal cells
    given a weight np.sqrt(2). Node ids are integer and correspond to 
    pixels starting at upper left corner, moving left to right, top to bottom.
    
    Parameters
    ---------- 
    pixels: numpy 2d array representing raster cell values
    cellSize: numeric representing raster cell size
    
    Returns
    ---------- 
    edges: list of tuples. each tuple has three elements, two integer node ids 
    and an edge weight
    idall: list of node ids corresponding to original pixel order. these
    do not include node ids of pixels with no data values. however, node
    ids are not renumbered so there can be gaps in numbering between ids.
    new_idmap: dictionary mapping original nodeids to node ids created
    after removing no data nodes and then renumbering remaining nodes
    consecutively        
    """
    # Convert coordinates to index function
    # Takes zero indexed cell coordinates
    # and converts to zero indexed index
    # that runs top left to right and down
    def cIndex(r, c, ncol):
        cind = ((c+1) + ncol*(r))-1
        return cind
        
    # Get number of rows and columns in array
    rows, cols = pixels.shape
    # Valid node id indexer (ids of nodes with non-nodata weights)
    vnid = 0
    for r, c, value in iArray:
        r = int(r)
        c = int(c)
        edges = []
        if value >= 1:
            #node = r * cols + c
            #node = np.where((iarray[:,0] == r) & (iarray[:,1] == c))[0][0]
            if c < cols - 1:  # Edge to the pixel to the right
                if pixels[r, c+1] >= 1: # check if neighbor is nodata
                    weight = np.mean(np.array([value, pixels[r, c+1]])) * cellSize
                    #nodenb = np.where((iarray[:,0] == r) & (iarray[:,1] == c+1))[0][0]
                    edges.append((vnid, vnid+1, weight))
            if r < rows - 1:  # Edge to the pixel below
                if pixels[r+1, c] >= 1: # check if neighbor is nodata
                    weight = np.mean(np.array([value, pixels[r+1, c]])) * cellSize
                    nodenb = keyArray[cIndex(r+1, c, cols)]
                    #nodenb = np.where((iarray[:,0] == r+1) & (iarray[:,1] == c))[0][0]
                    edges.append((vnid, nodenb, weight))
            if c > 0 and r < rows - 1: # Edge to the pixel below left
                if pixels[r+1, c-1] >= 1: # check if neighor is nodata
                    weight = np.mean(np.array([value, pixels[r+1, c-1]])) * cellSize * np.sqrt(2)
                    #nodenb = np.where((iarray[:,0] == r+1) & (iarray[:,1] == c-1))[0][0]
                    nodenb = keyArray[cIndex(r+1, c-1, cols)]
                    edges.append((vnid, nodenb, weight))
            if c < cols-1 and r < rows -1: # Edge to pixel below right
                if pixels[r+1, c+1] >= 1: # check if neighor is nodata
                    weight = np.mean(np.array([value, pixels[r+1, c+1]])) * cellSize * np.sqrt(2)
                    nodenb = keyArray[cIndex(r+1, c+1, cols)]
                    #nodenb = np.where((iarray[:,0] == r+1) & (iarray[:,1] == c+1))[0][0]
                    edges.append((vnid, nodenb, weight))

            yield vnid, edges
            # Increment nodeid
            vnid += 1

# Write data to file at index
def pwriteZarr(data, zarray, index):
    """

    Parameters
    ----------
    data : array
        numpy array to write to disk
    zarray : zarr file object
        file object to write 
    index : int
        location for writing data

    Returns
    -------
    str
        message indicating data written to file

    """
    zarray[index,:] = data
    return f"wrote index {index} to file"

# Create blosc file with compression
def createZarrBlosc(shape, dtype, chunks, store_path):
    """
    Parameters
    ----------
    shape : tuple of ints
        shape of array to create
    dtype : numpy data type object
        e.g. np.float64
    chunks : tuple of ints
        size of chunks in zarr file
    store_path : string
        file path

    Returns
    -------
    z_array : zarr file object
        zarr file object to use for reading and writing to disk
    """
    # Specify compression
    blosc_compressor = BloscCodec(cname="lz4", clevel=5, shuffle=BloscShuffle.bitshuffle)
    z_array = zarr.create_array(shape=shape, 
                         dtype=dtype, 
                         chunks=chunks, 
                         store=store_path,
                         compressors=blosc_compressor,
                         overwrite=True) # overwrite=True allows overwriting if the path exists
    return z_array

@njit
def idSkip(iArray, cols):
    """
    Parameters
    ----------
    iArray : float array
        3 columns (row, column, value) giving pixel indices
        and resistance values
    cols : int
        number of columns in resistance raster

    Returns
    -------
    list
        original nodeids accounting for gaps associated
        with non-valid cells
    """
    return [int(((x[1]+1)+(cols*((x[0]+1)-1)))-1) for x in iArray]


def calcPathsMod(x, zFileName, corrTolerance):
    """
    Function to calculate corridors/paths from one pair of distance arrays
    stored in a file on disk. Takes a node pair, sBatches, as input (as an array).
    The other inputs, nodidsLen, and corrTolerance are 
    integer, string, and integer respectively.
    Updated version of calcpaths that will probably replace the
    original.
    ----------
    x : two column array of integer node ids 
    zFileName : zarr file name on disk (string)
    corrTolerance : numeric
        
    Returns
    -------
    ccounts : numpy array
        array holding corridor counts
    """
    # Empty array to hold corridor counts
    # Open zarr file holding distances
    zfile = zarr.open(zFileName)
    lcc = zfile[x[0],:] + zfile[x[1],:]
    lcc[np.isinf(lcc)] = np.nan
    lcpVal = np.nanmin(lcc)
    lcc = np.where(lcc <= lcpVal + corrTolerance + 0.001, 1, 0)
    lcc[np.isnan(lcc)] = 0
    return lcc

