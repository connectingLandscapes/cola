# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:09:43 2023
Functions for creating a least cost path density
surface from xy locations.
@author: pj276
"""

#%%
# Imports
import networkit as nk
import networkx as nx
import rasterio as rio
from rasterio.crs import CRS
import numpy as np, time
from scipy import sparse
from scipy.sparse import diags
from numba import njit, prange
    
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
    profile.update(tiled=True, compress='lzw')
    with rio.open(outfile, 'w', **profile) as dst:
        dst.write(inarray)

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
    Uses a Networkx graph as an intermediate.
    Outputs node ids from the full graph (i.e. with no data cells) and
    a dictionary mapping node ids from the full graph to the subset
    graph (i.e. without no data cells)
    
    Parameters
    ----------
    src : rasterio dataset reader object
        Rasterio dataset object representing a 2D raster file.
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
    # Get indices of source points
    ##cinds = [ds.index(i[1].X,i[1].Y) for i in xys.iterrows()]
    # Convert to 1D for networkit
    ##sources = [(i[0]*ds.shape[1])+i[1] for i in cinds]
    # Get array    
    A = src_data
    # Convert to 1D
    A = A.flatten()
    # Convert resistance surface to adjacency matrix
    adj8 = connected_adjacency(src_data, str(pixelConnectivity)).astype("float32").tocsr()
    # Get indices of all cells
    sinds = sparse.find(adj8)
    # Calculate cost weights as average for adjacent cells
    adj8[sinds[0],sinds[1]] = ((A[sinds[0]] + A[sinds[1]])/2)*sinds[2]*cellSize
    # Convert to graph using networkx utility
    G8 = nx.from_scipy_sparse_array(adj8)
    # Set node attributes using cell values 
    nx.set_node_attributes(G8, dict(zip(range(0,G8.number_of_nodes()), [{'value': i} for i in A])))
    # Get node values
    node_values = nx.get_node_attributes(G8, 'value')
    # Remove no data nodes
    G8.remove_nodes_from((n for n, w in node_values.items() if w == ndValue))
    # Get list of original node ids
    nid = list(G8.nodes)
    # Convert networkx graph to networkit graph
    nG = nk.nxadapter.nx2nk(G8, weightAttr='weight')
    # Index edges
    nG.indexEdges()
    # Get mapping from node ids in G8 to node ids in nG
    # This relies on the fact that networkx preserves nodeids while 
    # networkit node ids range simply from 0 to N-1
    idmap = dict((id, u) for (id, u) in zip(G8.nodes(), range(G8.number_of_nodes())))
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

def rescale(OldValue, OldMin, OldMax, NewRMax, OldNdval, NewNdval, DoExpTrans=False, Cshape=0.1, NewMin=0, NewMax=1):
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
        min value of new range
    NewMax : numeric
        max value of new range

    Returns
    -------
    NewValue : numeric
        rescaled value
    """
    if OldValue == OldNdval:
        NewValue = NewNdval
    elif OldValue > OldMax:
        NewValue = NewMax
    elif OldValue < OldMin:
        NewValue = NewMin
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

