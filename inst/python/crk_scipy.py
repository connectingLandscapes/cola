# -*- coding: utf-8 -*-
"""
Cumulative resistant kernel mapping using scipy sparse graphs.

Replaces NetworKit SPSP with scipy.sparse.csgraph.dijkstra, which
accepts a distance limit parameter and avoids computing distances
beyond the threshold. Uses ThreadPoolExecutor with a semaphore and
sliding window to give independent control over CPU parallelism
(nThreads) and peak memory (maxLive arrays simultaneously in memory).

Peak memory ~ graph_memory + maxLive * M * 12 bytes
where M is the number of valid (non-nodata) cells and 12 bytes
accounts for the float64 Dijkstra vector + float64 kernel row
per active thread.

Usage:
    python crk_scipy.py <xyf> <rg> <ofile> <dThreshold> <tForm>
                        <tkv> <kvol> <nThreads> <maxLive> <upCRS>

Arguments:
    xyf        : path to CSV or shapefile of source point coordinates
    rg         : path to resistance raster (.tif)
    ofile      : path for output kernel raster (.tif)
    dThreshold : cost distance cutoff (float or int)
    tForm      : kernel transform: linear, inverse, inversesquare, gaussian
    tkv        : apply kernel volume transform: yes or no
    kvol       : kernel volume (float or int, ignored if tkv is no)
    nThreads   : number of parallel threads (controls CPU use)
    maxLive    : max arrays in memory at once (controls RAM use)
                 set lower than nThreads to reduce peak memory
                 e.g. nThreads=6, maxLive=2 uses 6 cores but caps
                 concurrent allocations at 2
    upCRS      : user-supplied CRS string e.g. ESRI:102028, or None
"""

import sys
import time
import threading
import itertools
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio as rio
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra

import cola_functions as cf


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def raster_to_csr(r, cSize):
    """
    Convert a 2D resistance raster to a CSR sparse graph.

    Valid cells are those with resistance >= 1 (after checkRasterVals)
    and not -9999 nodata. Edge weight between adjacent cells is the
    mean resistance multiplied by the Euclidean distance between cell
    centres (cSize for cardinal neighbours, cSize*sqrt(2) for diagonal).

    Only 4 forward-looking edges are generated per cell (E, S, SW, SE)
    matching the convention in generate_edgesMod. scipy dijkstra with
    directed=False treats these as undirected automatically.

    Returns
    -------
    graph : scipy.sparse.csr_matrix
    nodeids : np.ndarray int
        Flat indices of valid cells, used to map results back to raster.
    old_to_new : np.ndarray int, shape (nrows*ncols,)
        Maps original flat index to compact graph node id (-1 if invalid).
    """
    nrows, ncols = r.shape
    valid = (~np.isnan(r)) & (r != -9999) & (r >= 1)
    valid_flat = valid.ravel()
    nodeids = np.where(valid_flat)[0].astype(np.int32)
    M = len(nodeids)

    old_to_new = np.full(nrows * ncols, -1, dtype=np.int32)
    old_to_new[nodeids] = np.arange(M, dtype=np.int32)

    # float64 to match NetworKit internal precision
    flat = r.ravel().astype(np.float64)
    DIAG = np.sqrt(2) * float(cSize)
    CARD = float(cSize)

    # 4 forward-looking directions only: E, S, SW, SE
    # Matches generate_edgesMod; dijkstra(directed=False) handles symmetry
    offsets = [
        ( 0,  1, False),   # E
        ( 1,  0, False),   # S
        ( 1, -1, True),    # SW
        ( 1,  1, True),    # SE
    ]

    row_list, col_list, dat_list = [], [], []

    for dr, dc, is_diag in offsets:
        R, C = np.meshgrid(np.arange(nrows), np.arange(ncols), indexing='ij')
        r_dst = R + dr
        c_dst = C + dc

        in_bounds = (
            (r_dst >= 0) & (r_dst < nrows) &
            (c_dst >= 0) & (c_dst < ncols)
        )

        rs = R[in_bounds];  cs = C[in_bounds]
        rd = r_dst[in_bounds]; cd = c_dst[in_bounds]

        i_old = rs * ncols + cs
        j_old = rd * ncols + cd

        both = valid_flat[i_old] & valid_flat[j_old]
        i_old, j_old = i_old[both], j_old[both]

        i_new = old_to_new[i_old]
        j_new = old_to_new[j_old]

        w = (flat[i_old] + flat[j_old]) * 0.5
        w *= DIAG if is_diag else CARD

        row_list.append(i_new)
        col_list.append(j_new)
        dat_list.append(w)

    rows = np.concatenate(row_list)
    cols = np.concatenate(col_list)
    data = np.concatenate(dat_list).astype(np.float32)

    graph = sp.csr_matrix((data, (rows, cols)), shape=(M, M))
    return graph, nodeids, old_to_new


# ---------------------------------------------------------------------------
# Kernel transform
# ---------------------------------------------------------------------------

def apply_kernel_transform(d, dThreshold, tForm, tkv, kvol):
    """
    Transform a 1D distance array into kernel weights.

    Parameters
    ----------
    d : np.ndarray float64
        Distance values. NaN and values > dThreshold become 0.
    dThreshold : float
    tForm : str
        One of 'linear', 'inverse', 'inversesquare', 'gaussian'.
    tkv : str
        'yes' to apply volume transform.
    kvol : float
        Kernel volume multiplier (used when tkv == 'yes').

    Returns
    -------
    k : np.ndarray float64, same shape as d
    """
    k = np.zeros_like(d)

    if tForm == 'linear':
        k = 1.0 - (1.0 / dThreshold) * d
        k[k < 0] = 0.0

    elif tForm == 'inverse':
        with np.errstate(divide='ignore', invalid='ignore'):
            k = 1.0 / (d + 1.0)
        k[k < 1.0 / (dThreshold + 1.0)] = 0.0

    elif tForm == 'inversesquare':
        with np.errstate(divide='ignore', invalid='ignore'):
            k = 1.0 / (d ** 2 + 1.0)
        k[k < 1.0 / (dThreshold ** 2 + 1.0)] = 0.0

    elif tForm == 'gaussian':
        dispScale = dThreshold / 4.0
        k = np.exp(-((d ** 2) / (2.0 * dispScale ** 2)))
        k[k < np.exp(-((dThreshold ** 2) / (2.0 * dispScale ** 2)))] = 0.0

    else:
        raise ValueError(
            f"Unknown transform '{tForm}'. "
            "Choose linear, inverse, inversesquare, or gaussian."
        )

    k[np.isnan(d)] = 0.0

    if tkv == 'yes':
        k *= kvol * 3.0 / (np.pi * dThreshold ** 2)

    return k


# ---------------------------------------------------------------------------
# Per-source worker
# ---------------------------------------------------------------------------

def _single_kernel(graph, src, dThreshold, tForm, tkv, kvol):
    """
    Run Dijkstra from a single source and return the kernel contribution.
    Graph is read-only shared memory across threads -- no copying occurs.
    """
    d = dijkstra(graph, indices=int(src), directed=False,
                 limit=dThreshold, return_predecessors=False)
    d[d == np.inf] = np.nan
    return apply_kernel_transform(d.astype(np.float64), dThreshold,
                                  tForm, tkv, kvol)


# ---------------------------------------------------------------------------
# Parallel kernel accumulation
# ---------------------------------------------------------------------------

def _accumulate_kernels(graph, sources, dThreshold, tForm, tkv, kvol,
                        nThreads, maxLive):
    """
    Accumulate kernels across all sources with independent control
    over CPU parallelism and peak memory.

    Parameters
    ----------
    nThreads : int
        Number of concurrent threads (controls CPU utilisation).
    maxLive : int
        Maximum number of Dijkstra+kernel arrays in memory at once.
        Set lower than nThreads to reduce peak RAM. Threads block
        at the semaphore until a slot is free.

    Uses a sliding window (one submission per completion) so that
    exactly maxLive futures exist at any time and no sources are
    skipped regardless of how nThreads and maxLive relate to the
    total number of sources.
    """
    M = graph.shape[0]
    kernel = np.zeros(M, dtype=np.float64)
    semaphore = threading.Semaphore(maxLive)
    lock = threading.Lock()
    count = 0

    def _worker(src):
        semaphore.acquire()
        try:
            return _single_kernel(graph, src, dThreshold, tForm, tkv, kvol)
        finally:
            semaphore.release()

    with ThreadPoolExecutor(max_workers=nThreads) as executor:
        source_iter = iter(sources)
        pending = set()

        # Seed with exactly maxLive tasks using islice to avoid
        # consuming an extra item from the iterator
        for src in itertools.islice(source_iter, maxLive):
            pending.add(executor.submit(_worker, src))

        # Sliding window: submit one new task per completion
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                with lock:
                    kernel += future.result()
                    count += 1
                    if count % 500 == 0 or count == len(sources):
                        print(f'  {count}/{len(sources)} sources done',
                              flush=True)
                # Refill one slot for each completed future
                try:
                    src = next(source_iter)
                    pending.add(executor.submit(_worker, src))
                except StopIteration:
                    pass

    return kernel


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    tic = time.perf_counter()

    # ------------------------------------------------------------------
    # Parse arguments
    # ------------------------------------------------------------------
    xyf        = sys.argv[1]
    rg         = sys.argv[2]
    ofile      = sys.argv[3]
    dThreshold = sys.argv[4]
    tForm      = sys.argv[5]
    tkv        = sys.argv[6]
    kvol       = sys.argv[7]
    nThreads   = sys.argv[8]
    maxLive    = sys.argv[9]
    upCRS      = sys.argv[10]

    try:
        kvol = float(kvol) if cf.is_float(kvol) else int(kvol)
    except ValueError:
        sys.exit('Kernel volume must be a float or integer.')

    try:
        dThreshold = float(dThreshold) if cf.is_float(dThreshold) else int(dThreshold)
    except ValueError:
        sys.exit('Threshold must be a float or integer.')

    try:
        nThreads = int(nThreads)
    except ValueError:
        sys.exit('Number of threads must be an integer.')

    try:
        maxLive = int(maxLive)
    except ValueError:
        sys.exit('maxLive must be an integer.')

    if maxLive > nThreads:
        print(f'Warning: maxLive ({maxLive}) > nThreads ({nThreads}), '
              f'capping maxLive at nThreads.', flush=True)
        maxLive = nThreads

    # ------------------------------------------------------------------
    # Read resistance raster
    # ------------------------------------------------------------------
    r, profile = cf.read2flt32array(upCRS, rg)
    r, profile = cf.checkNoData(r, profile)
    r, profile = cf.checkRasterVals(r, profile)

    if profile['transform'][0] != np.abs(profile['transform'][4]):
        sys.exit('X and Y cell dimensions must be equal.')

    cSize = profile['transform'][0]

    # ------------------------------------------------------------------
    # Build CSR graph
    # ------------------------------------------------------------------
    print('Building graph...', flush=True)
    graph, nodeids, old_to_new = raster_to_csr(r, cSize)
    M = graph.shape[0]
    graph_gb = (graph.data.nbytes + graph.indices.nbytes +
                graph.indptr.nbytes) / 1e9
    peak_gb = graph_gb + maxLive * M * 12 / 1e9
    print(f'Nodes:  {M:,}', flush=True)
    print(f'Edges:  {graph.nnz:,}', flush=True)
    print(f'Graph:  {graph_gb:.2f} GB', flush=True)
    print(f'Per-thread array cost: {M * 12 / 1e6:.0f} MB', flush=True)
    print(f'Estimated peak memory: {peak_gb:.2f} GB '
          f'(graph + {maxLive} live arrays)', flush=True)

    # ------------------------------------------------------------------
    # Read source points
    # ------------------------------------------------------------------
    if Path(xyf).suffix == '.csv':
        xy = pd.read_csv(xyf)
    elif Path(xyf).suffix == '.shp':
        xy = gpd.read_file(xyf)
        xy = xy.get_coordinates()
    else:
        sys.exit('Source file must be .csv or .shp.')

    with rio.open(rg) as src:
        cinds = cf.cell_indices_from_coords(src, r, np.array(xy))
        if np.sum(np.isnan(cinds)) >= 1:
            sys.exit('Source points do not intersect resistance grid.')

        nPts = len(cinds)

        # Bounds checks
        cinds = cinds[cinds[:, 0] < r.shape[0]]
        cinds = cinds[cinds[:, 1] < r.shape[1]]
        cinds = cinds[cinds[:, 0] >= 0]
        cinds = cinds[cinds[:, 1] >= 0]

        # Remove no-data locations
        checkND = r[cinds[:, 0].astype(int), cinds[:, 1].astype(int)]
        cinds = cinds[checkND != -9999]

        if len(cinds) == 0:
            sys.exit('No valid source points after removing no-data values.')

        if nPts - len(cinds) > 0:
            print(f'Ignoring {nPts - len(cinds)} source(s) on no-data.',
                  flush=True)

        cinds = np.unique(cinds, axis=0)

        # Raster (row, col) -> flat index -> compact graph node id
        flat_ids = cinds[:, 0] * r.shape[1] + cinds[:, 1]
        sources  = old_to_new[flat_ids.astype(int)]

        valid_src = sources >= 0
        if not np.all(valid_src):
            print(f'Dropping {(~valid_src).sum()} source(s) on invalid cells.',
                  flush=True)
            sources = sources[valid_src]

        sources = sources.astype(np.int32)
        print(f'Valid sources: {len(sources):,}', flush=True)

    # ------------------------------------------------------------------
    # Accumulate kernels
    # ------------------------------------------------------------------
    print(f'Calculating kernels '
          f'(nThreads={nThreads}, maxLive={maxLive})...', flush=True)

    kernel = _accumulate_kernels(
        graph, sources, dThreshold, tForm, tkv, kvol, nThreads, maxLive)

    # ------------------------------------------------------------------
    # Map compact node ids back to full raster
    # ------------------------------------------------------------------
    dArr = np.zeros(r.shape[0] * r.shape[1], dtype=np.float32)
    dArr[nodeids] = kernel.astype(np.float32)
    dArr = dArr.reshape(r.shape[0], r.shape[1])

    # ------------------------------------------------------------------
    # Quantile summary
    # ------------------------------------------------------------------
    pos = dArr[dArr > 0]
    if len(pos) > 0:
        qDf = pd.DataFrame({
            'q':     np.arange(0, 1.01, 0.01),
            'value': np.quantile(pos, np.arange(0, 1.01, 0.01))
        })
        qDf.to_csv(ofile.replace('.tif', '_quantiles.csv'), index=False)

    # ------------------------------------------------------------------
    # Write output raster
    # ------------------------------------------------------------------
    cf.arrayToGeoTiff(np.expand_dims(dArr, axis=0), ofile, profile)
    print(f'Written: {ofile}', flush=True)

    toc = time.perf_counter()
    print(f'Total time: {(toc - tic) / 60:.2f} minutes', flush=True)


if __name__ == '__main__':
    main()
