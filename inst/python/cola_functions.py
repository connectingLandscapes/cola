# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:09:43 2023
Functions for creating a least cost path density
surface from xy locations.
@author: pj276
"""

#%%
# Imports
import sys
import networkit as nk
import rasterio as rio
from rasterio.crs import CRS
import numpy as np, time
from scipy import sparse
from scipy.sparse import diags
from numba import njit, prange
from scipy.sparse import tril
from shapely.geometry import MultiPolygon
import geopandas as gpd
import tables as tb
    
#%%
def asciiToGeoTiff(infile, outfile, crs=None):
    """
    Convert UNICOR or other ascii file to geotiff,
    optionally assigning a crs.
    
    Parameters
    ----------
    infile : str
        Path to input file
    outfile : str
        Path where output file will be saved
    crs : str
        EPSG or ESRI code in string format. e.g. "ESRI:102028"
    
    Returns
    -------
    None
    """
    # Set raster projection and save to file
    # Read raster and specify projection (it's not specified in the rsg formatted file)
    with rio.open(infile, 'r+') as src:
        if crs:
            src.crs = CRS.from_string(crs)
        profile = src.profile
        data = src.read()
    # Update driver to geotiff
    profile.update(driver="GTiff")
    # Save as tiff
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(data)

def arrayToGeoTiff(inarray, outfile, profile, crs=None, driver=None, res=None, dtype=None):
    """
    Save array as geotiff using exsiting rasterio profile.
    
    Parameters
    ----------
    inarray : numpy array
        3d array to write to tiff, dim 1 is number of bands,
        dim 2 is rows, dim 3 is columns
    outfile : str
        name of file to save, with extension
    profile : rasterio.profiles.Profile
        dictionary describing the characteristics to use when writing the file
    crs : str
        epsg or esri code corresponding to the coordinate reference
        system to assign the output file. formatted as e.g.
        'ESRI:102028' or 'EPSG:4326'
    driver : str
        gdal supported driver to use when writing.
        should correspond to the file extension.
        e.g. use 'GTiff' driver when file extension is '.tif'
    res : int
        cell resolution (in meters) of output image. e.g. 30
    dtype : str
        data type of output image e.g. 'float32'
    
    Returns
    -------
    None
    """
    if crs:
        profile.update(crs=CRS.from_string(crs))
    if driver:
        profile.update(driver=driver)
    if res:
        profile['transform'][0] = res
        profile['transform'][4] = -res
    if dtype:
        profile.update(dtype=dtype)
    profile.update(tiled=False, compress='lzw')
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(inarray)

# Function to generate edges based on pixel adjacency and distance
# between cells
@njit
def generate_edges(pixels, cellSize):
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
    # Get number of rows and columns in array
    rows, cols = pixels.shape
    # Valid node id indexer (ids of nodes with non-nodata weights)
    vnid = 0
    # Node ids
    # original node ids (0 to ncells-1)
    # valid node ids (0 to ncells with no data - 1)    
    idall = []
    idvalid = []
    edges = []
    for r in range(rows):
        for c in range(cols):
            # Only process node if weight >= 1 (otherwise skip as assumed no data)
            if pixels[r, c] >= 1:
                node = r * cols + c
                if c < cols - 1:  # Edge to the pixel to the right
                    if pixels[r, c+1] >= 1: # check if neighbor is nodata
                        weight = np.mean(np.array([pixels[r, c], pixels[r, c+1]])) * cellSize
                        edges.append((node, r * cols + (c+1), weight))
                if r < rows - 1:  # Edge to the pixel below
                    if pixels[r+1, c] >= 1: # check if neighbor is nodata
                        weight = np.mean(np.array([pixels[r, c], pixels[r+1, c]])) * cellSize
                        edges.append((node, (r+1) * cols + c, weight))
                if c > 0 and r < rows - 1: # Edge to the pixel below left
                    if pixels[r+1, c-1] >= 1: # check if neighor is nodata
                        weight = np.mean(np.array([pixels[r, c], pixels[r+1, c-1]])) * cellSize * np.sqrt(2)
                        edges.append((node, (r+1) * cols + (c-1), weight))
                if c < cols-1 and r < rows -1: # Edge to pixel below right
                    if pixels[r+1, c+1] >= 1: # check if neighor is nodata
                        weight = np.mean(np.array([pixels[r, c], pixels[r+1, c+1]])) * cellSize * np.sqrt(2)
                        edges.append((node, (r+1) * cols + (c+1), weight))                    
                # Add to node map
                idall.append(node)
                idvalid.append(vnid)
                # Increment nodeid
                vnid += 1
    new_idmap = dict((zip(idall, idvalid)))
    # Renumber edges using idmap (if any nodes were invalid, this makes nodeids continuous)
    edges = [(new_idmap[i[0]], new_idmap[i[1]], i[2]) for i in edges]
    return edges, idall, new_idmap

def connected_adjacency(image, connect, patch_size=(1, 1)):
    """
    Creates an adjacency matrix from an image where nodes are considered adjacent 
    based on 4-connected or 8-connected pixel neighborhoods.
    Minor modification from 
    https://stackoverflow.com/questions/30199070/how-to-create-a-4-or-8-connected-adjacency-matrix
    Diagonal weights set to sqrt(2). Orthogonal weights set to 1.
    Calculation for conductance is then conductance value/cellres*weight

    Parameters
    ----------    
    image : 2 or 3 dim numpy array
    connect : string, either '4' or '8'
    patch_size : tuple (n,m) used if the image will be decomposed into 
                   contiguous, non-overlapping patches of size n x m. The 
                   adjacency matrix will be formed from the smaller sized array
                   e.g. original image size = 256 x 256, patch_size=(8, 8), 
                   then the image under consideration is of size 32 x 32 and 
                   the adjacency matrix will be of size 
                   32**2 x 32**2 = 1024 x 1024

    Returns
    -------
    adjacency matrix as a sparse matrix (type=scipy.sparse.csr.csr_matrix)
    """
    r, c = image.squeeze().shape
    r = int(r / patch_size[0])
    c = int(c / patch_size[1])
    if connect == '4':
        # constructed from 2 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.ones(c*(r-1))
        upper_diags = diags([d1, d2], [1, c])
        return upper_diags + upper_diags.T
    elif connect == '8':
        # constructed from 4 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.append([0], d1[:c*(r-1)])
        d3 = np.ones(c*(r-1))
        d4 = d2[1:-1]
        d4[d4==1] = 2.0**0.5
        upper_diags = diags([d1, d2, d3, d4], [1, c-1, c, c+1])
        return upper_diags + upper_diags.T
    else:
        raise ValueError('Invalid parameter \'connect\'={connect}, must be "4" or "8".'
                     .format(connect=repr(connect)))

def image_to_graph(src_data, cellSize, ndValue, pixelConnectivity=8):
    """
    Convert a raster image to a weighted Networkit graph.
    Outputs node ids from the full graph (i.e. with no data cells) and
    a dictionary mapping node ids from the full graph to the subset
    graph (i.e. without no data cells)
    
    Parameters
    ----------
    src_data : numpy array
        Array corresponding to rasterio dataset reader object.
        Should be integer data type.
    cellSize : int
        Width or length in meters of square raster cells from raster file.
    ndValue : int
        NoData value in raster file. Cells with these values are considered
        unreachable from other cells.
    pixelConnectivity: int
        Whether to use a 4 or 8 connected neighborhood when defining the connected_adjacency matrix

    Returns
    -------
    nG : networkit graph
        Network graph of connections between adjacent raster cells.
        Edge weights correspond to average of resistance values of adjacent cells.
        Diagonals are appropriately weighted.
    nid : list
        Node ids from full graph
    idmap: dictionary
        Dictionary mapping node ids from full graph to node ids from the subset graph.
        Used to map results using the subset graph back to the full graph.
    """
    tic = time.perf_counter()
    # Get array    
    A = src_data
    # Check if non-negative min value is less than 1
    if np.nanmin(A[A >= 0]) < 1:
        print('Resistance surface has values < 1. Please provide a surface where the minimum resistance value is 1.')
        sys.exit(0) 
    # Convert to 1D
    A = A.flatten()
    # Convert resistance surface to adjacency matrix
    adj8 = connected_adjacency(src_data, str(pixelConnectivity)).astype("float32").tocsr()
    # Get indices of all cells
    sinds = sparse.find(adj8)
    # Calculate cost weights as average for adjacent cells
    adj8[sinds[0],sinds[1]] = ((A[sinds[0]] + A[sinds[1]])/2)*sinds[2]*cellSize
    # Get lower triangle in triplet format
    ro, co, dat = sparse.find(tril(adj8))
    # Create graph (nodes only)
    nG = nk.Graph(src_data.size, weighted=True)
    # Add edges to graph
    for i,j in enumerate(ro):
        nG.addEdge(ro[i], co[i], w=dat[i], addMissing=False, checkMultiEdge=False)
    # Iterate through nodes and remove if nodata
    for u in nG.iterNodes():
        if A[u] < 1:
            nG.removeNode(u)
    # Get original node ids
    nid = [u for u in nG.iterNodes()]
    # Create mapping from original node ids to node ids
    # after removing no data
    idmap = nk.graphtools.getContinuousNodeIds(nG)
    # Renumber node ids to be continuous
    nG = nk.graphtools.getCompactedGraph(nG, idmap)
    # Index edges
    nG.indexEdges()
    toc = time.perf_counter()
    return nG, nid, idmap
    print(f"The process took {toc - tic:0.4f} seconds")

import skimage as ski
from scipy.sparse import coo_matrix
def sk_image_to_graph(src_data, cellSize, ndValue):
    tic = time.perf_counter()
    # Get array    
    A = src_data
    # Convert to 1D
    A = A.flatten()
    # Create mask for pixel graph
    mask = np.full(src_data.shape, True)
    #mask[src_data == ndValue] = False
    # Function to average pixel neighbors
    func = lambda x, y, z: ((x+y)/2)*z
    # Create pixel graph (adjacency matrix)
    adj = ski.graph.pixel_graph(src_data, mask=mask, edge_function = func, connectivity=2, spacing=cellSize)
    # Convert from csr to coo
    adj = adj[0].tocoo()  
    # Get lower triangle in triplet format
    ro, co, dat = sparse.find(tril(adj))
    # Remove duplicate edges
    #x = [(i,j,v) for i,j,v in zip(adj.row, adj.col, adj.data) if i < j]
    # Convert back to coo
    #adj = coo_matrix(([i[2] for i in x], ([i[0] for i in x], [i[1] for i in x])))
    adj = coo_matrix((dat, (ro, co)))
    # Convert to graph
    nG = nk.graph.GraphFromCoo(adj, weighted=True, directed=False, edgesIndexed=True)
    # Iterate through nodes and remove if nodata
    for u in nG.iterNodes():
        if A[u] == -9999:
            nG.removeNode(u)
    # Get original node ids
    nid = [u for u in nG.iterNodes()]
    # Create mapping from original node ids to node ids
    # after removing no data
    idmap = nk.graphtools.getContinuousNodeIds(nG)
    # Renumber node ids to be continuous
    nG = nk.graphtools.getCompactedGraph(nG, idmap)
    # Index edges
    nG.indexEdges()
    toc = time.perf_counter()
    return nG, nid, idmap
    print(f"The process took {toc - tic:0.4f} seconds")
    
def cell_indices_from_coords(src, src_data, coords, win=None):
    """
    Find row and column indices of spatial coordinates relative
    to a raster dataset.
    
    Parameters
    ----------
    src : rasterio DatasetReader object
    src_data : 2d numpy array,
        usually corresponding to the DatasetReader object
    coords : 2 column pandas.core.frame.DataFrame.
        First column contains coordinates from the X dimension (longitude)
        Second column contains coordinates from the Y dimension (latitude)
    win (optional): rasterio Window object

    Returns
    cinds : 2 column numpy.ndarray
        First column contains row indices
        Second column contains column indices
    ----------
    """
    # Process if there are any target values in the array
    if coords.ndim > 0: #coords.shape[0] > 0:
        # get rows and columns of coordinates
        # if a window (i.e. spatial subset) is specified
        if win:
            row, col = rio.transform.rowcol(src.window_transform(win), coords[:,0], coords[:,1])
        # otherwise, use the whole array
        else:
            row, col = rio.transform.rowcol(src.transform, coords[:,0], coords[:,1])
        cinds = np.array((row,col)).transpose()
    else:
        cinds = np.array(np.nan)
    return cinds

def sourceTargetPairs(sources, thList):
    """
    Get unique pairs of sources and targets,
    omitting self-self pairs and duplicates
    E.g. calculate only (1,2) if both (1,2) and (2,1)
    pairs are present. Don't calculate (1,1), (2,2), etc.
    
    Parameters
    ----------
    sources : List of networkit node ids corresponding to source points
    thList : List of lists. The top level list is the same length
        as sources and its index corresponds the sources
        node id at that position. The second level list corresponds
        to the nodes that are within the threshold distance of the
        corresponding source node.

    Returns
    ----------
    reOrder : 2D numpy array where the first column is the source
        point for least cost path mapping and teh second column
        is the target point
    
    """
    # Empty list to hold unfiltered source target pairs
    pairList = []
    # Iterate and remove self self pairs
    for i in range(0, len(sources)):
        source = sources[i]
        for target in thList[i]:
            if source != target:
                pairList.append((source, target))
    # Convert pair list to array
    pairArr = np.array(pairList)
    del pairList
    # Reorder pairs so that largest valued id is on the right
    reOrder = np.zeros((pairArr.shape[0],pairArr.shape[1]), dtype='int32')
    reOrder[:,0] = np.min(pairArr,axis=1)
    reOrder[:,1] = np.max(pairArr,axis=1)
    del pairArr
    # Get unique rows
    reOrder = np.unique(reOrder, axis=0)
    return reOrder

def is_float(numString):
    """
    Check string representaion of a number for a decimal point
    and if subsequent conversion to numeric is valid
    
    Parameters
    ----------
    numString : String representation of a number

    Returns
    ----------
    True or False
    
    """
    if "." in numString and numString.replace(".", "").isnumeric():
        return True
    else:
        return False

def is_floatpy3(element: any) -> bool:
    """
    From https://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python/20929881#20929881
    Parameters
    ----------
    element : any
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    #If you expect None to be passed:
    if element is None: 
        return False
    try:
        float(element)
        return True
    except ValueError:
        return False

def rescale(OldValue, OldMin, OldMax, NewRMax, OldNdval, NewNdval, DoExpTrans=False, Cshape=0.1, NewMin=0, NewMax=1, NewRMin=1):
    """
    Rescales numbers using old and new ranges as references. If old value
    is outside the range of the original range, it is set to the corresponding
    new max or min value.
    Parameters
    ----------
    OldValue : numeric
        value to be rescaled
    OldMin : numeric
        min value of original range
    OldMax : numeric
        max value of original range
    NewRMax : numeric
        new rescaled maximum value
    OldNdval: numeric
        old no data value
    NewNdval: numeric
        new no data value
    DoExpTrans: True/False
        whether to use exponential transform after rescaling
    Cshape: numeric
        controls shape of transform: > 0 resistance increases slowly with decrease in quality
        < 0 resistance increases rapidly with decrease in quality
        ~ 0 resistance increases linearly with decrease in quality
    NewMin : numeric
        min value of new range (before exponential rescale)
    NewMax : numeric
        max value of new range (before exponential rescale)
    NewRMin : numeric
        min value of new range (after exponential rescale)

    Returns
    -------
    NewValue : numeric
        rescaled value
    """
    if OldValue == OldNdval:
        NewValue = NewNdval
    elif OldValue > OldMax:
        if DoExpTrans == False:
            NewValue = NewMax
        else:
            NewValue = NewRMin
    elif OldValue < OldMin:
        if DoExpTrans == False:
            NewValue = NewMin
        else:
            NewValue = NewRMax
    else:
        OldRange = (OldMax - OldMin)
        if OldRange == 0:
            NewValue = NewMin
        else:
            NewRange = (NewMax - NewMin)  
            NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
            if DoExpTrans == True:
                NewValue = NewRMax-(NewRMax-1)*((1-np.exp(-Cshape*NewValue))/(1-np.exp(-Cshape)))
    return(NewValue)
        
# Vectorize rescale to work with numpy arrays
vecrescale = np.vectorize(rescale)

# From https://stackoverflow.com/questions/3104781/inverse-distance-weighted-idw-interpolation-with-python
def simple_idw(x, y, z, xi, yi):
    dist = distance_matrix(x,y, xi,yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / dist

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point by all observed Z-values
    zi = np.dot(weights.T, z)
    return zi

def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    # (Yay for ufuncs!)
    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])

    return np.hypot(d0, d1)

def sumLccBatch(reO, s0, s1, sB, ccA, cTol):
    """
    Loops through batches of source target pairs,
    calculates least cost corridor between each,
    converts to categorical (0,1) raster, then sums.
    
    Parameters
    ----------
    reO : 2D numpy array where the first column is the source
        point for least cost path mapping and the second column
        is the target point
    s0 : integer
        start index of batch
    s1 : integer
        end index of batch
    sB : 1D numpy array holding unique values of sources 
        within the current batch of source pairs
    ccA : 2D numpy array holding cost distances from each
        source point to every cell in the landscape
    cTol : numeric
        cost distance to add to the least cost path (lcp).
        values <= the lcp + cTol will be extracted from
        the cumulative cost surface. the lcp is returned
        if cTol is 0, otherwise a corridor around and
        including the lcp is returned

    Returns
    -------
    lccSum : 1D numpy array holding the sum of overlapping
        corridors or lcps within the current batch
    """
    lccSum = np.zeros((ccA.shape[1],))
    for i in reO[s0:s1+1]:
        source = np.where(sB==i[0])[0][0]
        target = np.where(sB==i[1])[0][0]
        lcc = ccA[source,:]+ccA[target,:]
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        lcc = np.where(lcc <= lcpVal + cTol, 1, 0)
        lcc[np.isnan(lcc)] = 0
        lccSum += lcc
    return(lccSum)

def sumLccBatch2(reO, sB, ccA, cTol):
    """
    Loops through batches of source target pairs,
    calculates least cost corridor between each,
    converts to categorical (0,1) raster, then sums.
    
    Parameters
    ----------
    reO : 2D numpy array where the first column is the source
        point for least cost path mapping and the second column
        is the target point
    sB : 1D numpy array holding unique values of sources 
        within the current batch of source pairs
    ccA : 2D numpy array holding cost distances from each
        source point to every cell in the landscape
    cTol : numeric
        cost distance to add to the least cost path (lcp).
        values <= the lcp + cTol will be extracted from
        the cumulative cost surface. the lcp is returned
        if cTol is 0, otherwise a corridor around and
        including the lcp is returned

    Returns
    -------
    lccSum : 1D numpy array holding the sum of overlapping
        corridors or lcps within the current batch
    """
    lccSum = np.zeros((ccA.shape[1],))
    for i in reO:
        source = np.where(sB==i[0])[0][0]
        target = np.where(sB==i[1])[0][0]
        lcc = ccA[source,:]+ccA[target,:]
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        lcc = np.where(lcc <= lcpVal + cTol, 1, 0)
        lcc[np.isnan(lcc)] = 0
        lccSum += lcc
    return(lccSum)

def newSourceTargetPairs(thList):
    """
    Get unique pairs of sources and targets,
    omitting self-self pairs and duplicates
    E.g. calculate only (1,2) if both (1,2) and (2,1)
    pairs are present. Don't calculate (1,1), (2,2), etc.
    Uses numpy functions instead of lists used in 
    original sourceTargetPairs function
    
    Parameters
    ----------
    thList : List of lists. The top level list is the same length
        as sources and its index corresponds the sources
        node id at that position. The second level list corresponds
        to the nodes that are within the threshold distance of the
        corresponding source node.

    Returns
    ----------
    reOrder : 2D numpy array where the first column is the source
        point for least cost path mapping and the second column
        is the target point
    
    """
    # Empty list to hold unfiltered source target pairs
    pairList = []
    for ind, a in enumerate(thList):
        if len(a[0]) > 0:
            pairList.append(np.vstack((np.repeat(ind, len(a[0])), a[0])).T)
    pairArr = np.concatenate(pairList)
    del pairList
    
    # Reorder pairs so that largest valued id is on the right
    reOrder = np.zeros((pairArr.shape[0],pairArr.shape[1]), dtype='int32')
    reOrder[:,0] = np.min(pairArr,axis=1)
    reOrder[:,1] = np.max(pairArr,axis=1)
    del pairArr
    # Get unique rows
    reOrder = np.unique(reOrder, axis=0)
    return reOrder

def inds2nodeids(cinds, r, idmap):
    """
    Check for point-raster overlap, filter out 
    nodata points, get unique indices of points
    and then map to networkit node ids
    Parameters
    ----------
    cinds : 2d numpy array with array indices
        corresponding to source point locations
    r : 2d numpy array extracted from raster dataset
        read in via rasterio
    idmap : dictionary
        Dictionary mapping node ids from full graph to node ids from the subset graph.
        Used to map results using the subset graph back to the full graph.

    Raises
    ------
    Exception
        If there are no valid source points.

    Returns
    -------
    sources : list
        List of networkit nodeids corresponding to points
    """
    # Get number of points supplied by user
    nPts = len(cinds)
    # Check for points and remove any that are out of raster bounds
    # Apparently this can happen using cell_indices_from_coords
    if any(cinds[:,0] >= r.shape[0]):
        cinds = cinds[cinds[:,0] < r.shape[0]]
    if any(cinds[:,1] >= r.shape[1]):
        cinds = cinds[cinds[:,1] < r.shape[1]]
    # Check for individual points with no data
    checkND = r[np.array(cinds[:,0]),np.array(cinds[:,1])]
    if -9999 in checkND:
        cinds = cinds[checkND != -9999,:]
    # If there are no points after removing no data vals,
    # exit script.
    if len(cinds) == 0:
        raise Exception('No valid source points. This may be because none of your source points overlap with your resistance grid.')
    # Check if any points were removed
    if nPts - len(cinds) > 0:
        print('Ignoring ' + str(nPts - len(cinds)) + " source point(s) with no data values.", flush=True)
    # Get unique indices
    cinds = np.unique(cinds,axis=0)
    # Convert to list
    sources = [(i[0]*r.shape[1])+i[1] for i in cinds]
    # Convert to networkit nodeids
    sources = [idmap[n] for n in sources]
    print("Number of valid sources = " + str(len(sources)), flush=True)
    return sources

@njit(parallel=True)
def sumLccBatch2Numba(reO, sB, ccA, cTol):
    """
    Loops through batches of source target pairs,
    calculates least cost corridor between each,
    converts to categorical (0,1) raster, then sums.
    
    Parameters
    ----------
    reO : 2D numpy array where the first column is the source
        point for least cost path mapping and the second column
        is the target point
    sB : 1D numpy array holding unique values of sources 
        within the current batch of source pairs
    ccA : 2D numpy array holding cost distances from each
        source point to every cell in the landscape
    cTol : numeric
        cost distance to add to the least cost path (lcp).
        values <= the lcp + cTol will be extracted from
        the cumulative cost surface. the lcp is returned
        if cTol is 0, otherwise a corridor around and
        including the lcp is returned

    Returns
    -------
    lccSum : 1D numpy array holding the sum of overlapping
        corridors or lcps within the current batch
    """
    lccSum = np.zeros((ccA.shape[1],))
    for i in prange(len(reO)):
        source = np.where(sB==reO[i][0])[0][0]
        target = np.where(sB==reO[i][1])[0][0]
        lcc = ccA[source,:]+ccA[target,:]
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        lcc = np.where(lcc <= lcpVal + cTol, 1, 0)
        lcc[np.isnan(lcc)] = 0
        lccSum += lcc
    return(lccSum)

def sources2nodeids(rg, r, xy, idmap):
    """
    Convert xy points to Networkit node ids
    Outputs a list of node ids
    
    Parameters
    ----------
    rg : string
        Name of raster used to build the Networkit graph for which node ids are requested
    r : numpy array
        Array corresponding to raster used to build the Networkit graph
    xy : pandas dataframe
        Two column dataframe with xy coordinates of points
    idmap : dictionary
        Dictionary mapping node ids from full graph to node ids from the subset graph.
        Used to map results using the subset graph back to the full graph.
        Outupt by the image_to_graph function
    Returns
    -------
    sources : list
        List of Networkit node ids corresponding to points
    """
    with rio.open(rg) as src:
        # Convert xy coordinates to rows and columns using
        # the resistance grid as a reference.
        cinds = cell_indices_from_coords(src, r, np.array(xy))
        # Raise an exception if there are no source points in the landscape
        if np.sum(np.isnan(cinds)) >= 1:
            raise Exception('Source points do not intersect resistance grid.')
        else:
            # Get number of points supplied by user
            nPts = len(cinds)
            # Check for points and remove any that are out of raster bounds
            # Apparently this can happen using cell_indices_from_coords
            if any(cinds[:,0] >= r.shape[0]):
                cinds = cinds[cinds[:,0] < r.shape[0]]
            if any(cinds[:,1] >= r.shape[1]):
                cinds = cinds[cinds[:,1] < r.shape[1]]
            # Check for individual points with no data
            checkND = r[np.array(cinds[:,0]),np.array(cinds[:,1])]
            if -9999 in checkND:
                cinds = cinds[checkND != -9999,:]
            # If there are no points after removing no data vals,
            # exit script.
            if len(cinds) == 0:
                raise Exception('No valid source points. This may be because none of your source points overlap with your resistance grid.')
            # Check if any points were removed
            if nPts - len(cinds) > 0:
                print('Ignoring ' + str(nPts - len(cinds)) + " source point(s) with no data values.", flush=True)
            # Get unique indices
            cinds = np.unique(cinds,axis=0)
            # Convert to list
            sources = [(i[0]*r.shape[1])+i[1] for i in cinds]
            # Convert to networkit nodeids
            sources = [idmap[n] for n in sources]
            print("Number of valid sources = " + str(len(sources)), flush=True)
            return sources

def checkNoData(r, profile):
    """
    Check no data value of a raster and change
    if not -9999
    Parameters
    ----------
    r : numpy array
        raster array
    profile : dictionary
        rasterio profile dictionary
    Returns
    -------
    r : numpy array
        raster array
    profile : dictionary
        rasterio profile dictionary
    """
    # Check no data value
    if profile['nodata'] != -9999:
        if profile['nodata'] == None:
            r[r==0] = -9999
            profile.update({'nodata':-9999})
            print('No data value was given as None')
            print('Assuming 0 corresponds to no data')
            print('Converting no data values to -9999.')
        elif np.isnan(profile['nodata']):
            r[np.isnan(r)] = -9999
            profile.update({'nodata':-9999})
            print('Converting no data values to -9999.')
        else:
            r[r==profile['nodata']] = -9999
            profile.update({'nodata':-9999})
            print('Converting no data values to -9999.')
        return([r, profile])
    else:
        return([r, profile])
            #raise Exception('No data value is not equal to -9999. Reformat resistance grid so that no data values are set to -9999.')   

def checkRasterVals(r, profile):
    """
    Check whether there are resistance values between 0-1
    Parameters
    ----------
    r : numpy array
        raster array
    profile : dictionary
        rasterio profile dictionary
    Returns
    -------
    r : numpy array
        raster array
    profile : dictionary
        rasterio profile dictionary
    """
    # Check resistance values
    if np.sum((r > 0) & (r < 1)) > 0:
        r[((r > 0) & (r < 1))] = 1
        print('Warning: Resistance values between 0 and 1 detected. Converting these values to 1.')
        print('If this is not what you want, please load a resistance layer with no resistance values less than 1.')
        return([r, profile])
    else:
        return([r, profile])
            #raise Exception('No data value is not equal to -9999. Reformat resistance grid so that no data values are set to -9999.') 

def read2flt32array(upCRS, rg):
    """
    Read resistance grid to array
    Assign crs if provided by user
    Output array as float32, update profile if needed
    Parameters
    ----------
    upCRS : string
        coordinate reference system, default is None
    rg : string
        grid file name to read to array
    Returns
    -------
    r : numpy array
        raster array
    profile : dictionary
        rasterio profile dictionary
    """
    if upCRS != 'None':
        with rio.open(rg, 'r+') as src:
            # Assign projection
            src.crs = CRS.from_string(upCRS)
            # Check if projection is geographic
            if src.crs.is_geographic:
                raise Exception('Coordinate reference system of input file is geographic. CoLa requires input files to have a projected, equal area coordinate reference system.')
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype != 'float32': #['byte','int8','uint8','intp','uintp']:
                r = r.astype('float32')
                profile['dtype'] = 'float32'
            # Convert profile to single band
            if profile['count'] > 1:
                profile['count'] = 1
    # Otherwise, read as usual
    else:            
        with rio.open(rg) as src:
            # Check if projection is geographic
            if src.crs.is_geographic:
                raise Exception('Coordinate reference system of input file is geographic. CoLa requires input files to have a projected, equal area coordinate reference system.')
            # Get profile
            profile = src.profile
            # Read to array
            r = src.read(1)
            # Convert to larger int if needed
            if r.dtype != 'float32': #['byte','int8','uint8','intp','uintp']:
                r = r.astype('float32')
                profile['dtype'] = 'float32'
            # Convert profile to single band
            if profile['count'] > 1:
                profile['count'] = 1
    return(r, profile)

def calcPaths(x, nodeidsLen, dahdf, sBatches, corrTolerance):
    """
    Function to calculate corridors/paths from pairs of distance arrays
    stored in an hdf file. Takes a batch of node pairs, sBatches, as input (as an array).
    For use in joblib parallel processing script lcc_joblib.py.
    The other inputs, nodidsLen, dahdf, and corrTolerance are 
    integer, string, and integer respectively.
    ----------
    x : integer used to index into 
        two column array of integer node ids
    Returns
    -------
    ccounts : numpy array
        array holding corridor counts
    """
    # Empty array to hold corridor counts
    ccounts = np.zeros((1,nodeidsLen))
    # Open hdf5 file holding distances
    h5f = tb.open_file(dahdf, 'r')
    # Loop through pairs and calculate corridors
    for i in sBatches[x]:
        lcc1 = h5f.root.dset[i[0],:].astype(np.float64)
        lcc1[lcc1 >= np.iinfo(np.uint64).max] = np.nan
        lcc2 = h5f.root.dset[i[1],:].astype(np.float64)
        lcc2[lcc2 >= np.iinfo(np.uint64).max] = np.nan
        lcc = (lcc1 + lcc2)/10000
        del lcc1, lcc2
        lcc[np.isinf(lcc)] = np.nan
        lcpVal = np.nanmin(lcc)
        lcc = np.where(lcc <= lcpVal + corrTolerance + 0.0001, 1, 0)
        lcc[np.isnan(lcc)] = 0
        ccounts += lcc
    h5f.close()
    return(ccounts)

def calcKernels(x, nodeidsLen, nodehdf, dahdf, sBatches, dThreshold, tForm, tkv, kvol):
    """
    Function to calculate kernels from distance arrays
    stored in an hdf file. Takes a batch of nodes as input (as an array).
    For use in joblib parallel processing script crk_joblib.py.
    ----------
    x : integer used to index into 
        two column array of integer node ids
    Returns
    -------
    ksums : numpy array
        array holding kernel sums
    """
    # Empty array to hold kernel sums
    ksums = np.zeros((1,nodeidsLen))
    # Open hdf5 file holding distances
    h5f = tb.open_file(dahdf, 'r')
    # Loop through cost distance arrays and calculate crks
    for i in sBatches[x]:
        ccArr = h5f.root.dset[nodehdf[i],:]
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
        # Add kernel to running sum
        ksums += ccArr
    h5f.close()
    return(ksums)
    
def groupby_multipoly(df, by, aggfunc="first"):
    data = df.drop(labels=df.geometry.name, axis=1)
    aggregated_data = data.groupby(by=by).agg(aggfunc)

    # Process spatial component
    def merge_geometries(block):
        return MultiPolygon(block.values)

    g = df.groupby(by=by, group_keys=False)[df.geometry.name].agg(
        merge_geometries
    )

    # Aggregate
    aggregated_geometry = gpd.GeoDataFrame(g, geometry=df.geometry.name, crs=df.crs)
    # Recombine
    aggregated = aggregated_geometry.join(aggregated_data)
    return aggregated

def zoneAdjacency(zones):
    """
    Get zone adjacency from a watershed segmentation
    Outputs array of zone neighbors. Neighbor combinations
    are not unique (i.e. there can be repeats) so further
    filtering is required if unique combinations are desired.
    Parameters
    ----------
    r : numpy array
        Array output by skimage watershed function
    Returns
    -------
    zOs : numpy array
        Two column array of zone adjacencies
    """
    # Right shift
    rs = np.pad(zones,((0,0),(1,0)), mode='constant')[:, :-1]    
    # Left shift
    ls = np.pad(zones,((0,0),(0,1)), mode='constant')[:, 1:]
    # Up shift
    us = np.pad(zones,((0,1),(0,0)), mode='constant')[1:,:]
    # Down shift
    ds = np.pad(zones,((1,0),(0,0)), mode='constant')[:-1,:]
    # Diagonal upper left shift
    dul = np.pad(zones,((0,1),(0,1)), mode='constant')[1:, 1:]
    # Diagonal upper right shift
    dur = np.pad(zones,((0,1),(1,0)), mode='constant')[1:, :-1]
    # Diagonal lower left shift
    dll = np.pad(zones,((1,0),(0,1)), mode='constant')[:-1, 1:]
    # Diagonal lower right shift
    dlr = np.pad(zones,((1,0),(1,0)), mode='constant')[:-1, :-1]
    # Overlay with original zones and get
    # unique combinations (excluding zero)
    aList = []
    for i in [rs,ls,us,ds,dul,dur,dll,dlr]:
        zO = np.unique(np.vstack((zones[(zones != i) & (zones != 0) & (i != 0)], i[(zones != i) & (zones != 0) & (i != 0)])).T, axis=0)
        if len(zO) > 0:
            aList.append(zO)
    zOs = np.vstack(aList)
    return zOs

