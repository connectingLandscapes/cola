# In Progress.
# Used networkit notebook for visualisation.
# User can upload a folder of outputs from cdpop or if cdpop is ran before this then
# those output will be populated in this page. Based on the selection of csv file, visualisation is being done

import streamlit as st
import pandas as pd
import rasterio as rio
import networkit as nk
import networkx as nx
from collections import Counter
import sys
import zipfile
import os
import numpy as np, time
from numpy import copy
from itertools import zip_longest, chain
from matplotlib import pyplot
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

# Sidebar for file uploads
csv_file = st.sidebar.file_uploader("Upload CSV file", type=["csv"])
tif_file = st.sidebar.file_uploader("Upload TIF file", type=["tif", "tiff"])



# Sidebar for folder (zip) upload
folder_zip = st.sidebar.file_uploader("Upload ZIP containing folders with CSV files", type=["zip"])

all_csv_files = []
if folder_zip:
    # Extract the zip file
    with zipfile.ZipFile(folder_zip, 'r') as zip_ref:
        zip_dir = os.path.join("uploaded_folder", folder_zip.name.replace(".zip", ""))
        zip_ref.extractall(zip_dir)
    
    # List all CSV files within subdirectories
    for subdir, _, files in os.walk(zip_dir):
        for file in files:
            if file.endswith(".csv"):
                all_csv_files.append(os.path.join(subdir, file))
    
    # Display the CSV files in a dropdown for the user to select
    selected_csv = st.sidebar.selectbox("Select a CSV file:", all_csv_files)
        # Placeholder for visualization code

    if csv_file and tif_file:
        def cell_indices_from_coords_nowindow(src, src_data, coords):
            # proceed if there are any target values in the array
            if coords.ndim > 0: #coords.shape[0] > 0:
                # get rows and columns of coordinates
                row, col = rio.transform.rowcol(src.transform, coords[:,0], coords[:,1])
                cinds = np.array((row,col)).transpose()
            else:
                cinds = np.array(np.nan)
            return cinds

        # Coordinate reference system
        # crsString = 'ESRI:102028'
        # Cell size
        cSize = 500
        # Distance threshold
        dThreshold = 2000000

        # xy coordinate file
        xyf = selected_csv
        xy = pd.read_csv(selected_csv)

        # resistance grid file
        rg = tif_file
        with rio.open(tif_file) as src:
            r = src.read(1)
            # Convert to int16
            r = copy(r).astype(np.int16)

        print("Converting image to graph", flush=True)
        tic = time.perf_counter()
        nkG, nodeids, idmap = rkf.image_to_graph(r, cSize, -9999)
        toc = time.perf_counter()
        print(f"Converting image to graph took {toc - tic:0.4f} seconds", flush=True)
        print("number of edges = ", str(nkG.numberOfEdges()), flush=True)
        print(nk.overview(nkG))


        # Get sources as rows and columns
        with rio.open(tif_file) as src:
            cinds = cell_indices_from_coords_nowindow(src, r, np.array(xy))
            # Convert to 1D
            sources = [(i[0]*r.shape[1])+i[1] for i in cinds]
            # Convert to networkit nodeids
            sources = [idmap[n] for n in sources]
            print("number of sources = " + str(len(sources)), flush=True)


        #nk.setNumberOfThreads(3)

        # Get distances from source nodes to every other node in the graph
        # Note: future versions of networkit should allow target nodes to be specified
        tic = time.perf_counter()
        aa = nk.distance.SPSP(nkG, sources)
        aa.run()
        bb = aa.getDistances()
        toc = time.perf_counter()
        print(f"Calculating shortest paths took {toc - tic:0.4f} seconds", flush=True)


        # Get distances between source nodes
        ccList = []
        for i in bb:
            ccList.append([i[s] for s in sources])

        # Get nodes that are less than the threshold distance from the target node
        thList = []
        distList = []
        for lst in ccList:
            thList.append([sources[i] for i,j in enumerate(lst) if j < dThreshold])
            distList.append([j for i,j in enumerate(lst) if j < dThreshold])


        # Iterate through target nodes and calculate paths to nodes within the threshold distance
        tic = time.perf_counter()
        pDict = {}
        for i in range(0,244):
            source = sources[i]   
            for target in thList[i]:
                if source != target:
                    dij = nk.distance.BidirectionalDijkstra(nkG, source, target, storePred=True)
                    dij.run()
                    sPs = dij.getPath()
                    # Add source and target nodes to path
                    sPs = sPs + [source] + [target]
                    tempDict = Counter(sPs)
                    # From https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries
                    pDict = {k: tempDict.get(k, 0) + pDict.get(k, 0) for k in set(tempDict) | set(pDict)}
                    del dij
        toc = time.perf_counter()
        print(f"Calculating shortest paths took {toc - tic:0.4f} seconds", flush=True)


        # The keys in the idmap dict are the original array indices, values are subset array indices
        print(np.max(list(idmap.keys())))
        # Invert dict mapping
        inv_map = {v: k for k, v in idmap.items()}
        # The keys in the inverted dict are the subset array indices, values are original array indices
        print(dict(list(inv_map.items())[0:2]))
        print(np.max(list(inv_map.keys())))
        # The key/valuse in pDict are the subset array indices and least cost path counts
        print(dict(list(pDict.items())[0:2]))
        print(np.max(list(pDict.keys())))
        # Iterate through dicts to get lists of original array indices and least cost path counts
        l12 = [(inv_map[k], v) for k, v in pDict.items()]
        l1, l2 = zip(*l12)

        # Create 1D zeros array to hold least cost path values
        dArr = np.zeros([r.shape[0]*r.shape[1]], 'float32')
        # Assign least cost path values
        dArr[l1,] = l2
        # Make into rectangular array
        dArr = dArr.reshape(r.shape[0],r.shape[1])

        # Use gaussian filter to smooth least cost paths
        # Sigma is standard deviation of the kernel
        # The size of the kernel on each side is 2*radius + 1
        lcpSmooth = gaussian_filter(dArr, sigma=1, radius=2)
        if selected_csv:
            # Visualization code
            plt.imshow(lcpSmooth)
            plt.colorbar()
            plt.title('Least Cost Paths Visualization')
            st.pyplot(plt)

            nk.getCurrentNumberOfThreads()

