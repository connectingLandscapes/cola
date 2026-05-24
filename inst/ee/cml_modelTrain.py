# -*- coding: utf-8 -*-
'''
Created on Tue Apr 21 09:23:51 2026
@author: ig299
'''
# =============================================================================
# from: https://developers.google.com/earth-engine/tutorials/community/species-d
# istribution-modeling#model_fitting_and_prediction
# Author(s): osgeokr
# =============================================================================

# IMPORTS
import ee
import sys
import csv
import geemap
import geemap.colormaps as cm
import pandas as pd, geopandas as gpd
import numpy as np, matplotlib.pyplot as plt
import requests, math, random
from ipyleaflet import TileLayer
from statsmodels.stats.outliers_influence import variance_inflation_factor


import os

# os.chdir('N:/My Drive/git/cola/inst/ee')
# C:/Users/ig299/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_covsExtraction.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv" 
# C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_covsExtraction.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv"
# C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_covsExtraction.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv" 


def main() -> None:    
    #%% Load EE
    def asset_exists(asset_id):
        try:
            ee.data.getAsset(asset_id)
            return True
        except ee.EEException:
            return False

    # INPUTS
    paramcsv = sys.argv[1] # project name
    # paramcsv = "N:/My Drive/git/cola/inst/ee/params_train.csv"
    
    try:
        eelogpath = str(sys.argv[2]) # 
        if not os.path.isdir(eelogpath):
            os.makedirs(eelogpath, exist_ok=True)
        if not os.path.isdir(eelogpath):
            eelogpath = ""
    except IndexError:
        eelogpath = ""
    # eelogpath = 'C:/Users/ig299/cola'
    #    
    #
    # LOAD PARAMS CSV
    with open(paramcsv) as f:
        dic = dict(filter(None, csv.reader(f)))
    dic
    #
    eeproject = dic.get('eeProject')
    
    
    #%% Start EE
    try:
        ee.Authenticate( )
        ee.Initialize(project = eeproject)
        print('  EE initialized')
    except:
        print(' ERROR: no EE initialized')
        exit(1)
    #
    points_path = dic.get('ptsWithCovs',False)
    if points_path is False:
        print(' ERROR: ptsWithCovs not provided')
        exit(1)
    #
    #
    # check if dataset+covs are EE or local file
    if points_path.find('projects/') == 0:
        try:
            ee.data.getAsset(points_path)
            print(' Points file was found in EE')
            covs_local = False
        except ee.EEException:
            print(' ERROR: points path was not found in EE')
            exit(1)
    elif os.path.splitdrive(points_path)[0] != '':
         if not os.path.isfile(points_path):
             print(' ERROR: local file was not found in EE')
             exit(1)
         else:
             covs_local = True
    #
    ###########################################
    #%% Load modules and create local variables 
    from sat_ts_fusion.imagery import lsat_utils
    from sat_ts_fusion.ccdc import ccdc_utils
    from sat_ts_fusion.fusion import ancillary_covariates as anc
    from sat_ts_fusion.fusion import covariates as covs
    from sat_ts_fusion import eecolatools as ct
    
    # from sat_ts_fusion.fusion.models import models

    out_ee_path = dic.get('outEEPath', '') # C:/Users/ig299/cola/
    out_local_path = dic.get('outLocalPath', '') # C:/Users/ig299/cola/
    ee_log_path = dic.get('eeLogPath','')
    showplots = dic.get('showPlots',False) #False
    npres = dic.get('nPres',False) #100
    nabs = dic.get('nAbs',False) #100
    nclus = dic.get('nClus',False) #2
    disfun = dic.get('distFun',False) #Euclidean
    rand_grid_size = dic.get('randGridSize',False) #50000
    spat_dup_pixsize = int(dic.get('spatDupPixSize',200)) #50000
    split = dic.get('splot', 0.7)
    numiter = dic.get('numiter', 1)
    stack_type = int(dic.get('stackType', 1))
    points_buff_meters = int(dic.get('pointsBuffMeters', 0))
    pixel_size = int(dic.get('pixelSize',False)) # ,100
    remove_spat_dup = dic.get('removeSpatDup',False)
    area_of_interest = dic.get('areaOfInterest', '')
    
    num_folds = int(dic.get('numFolds', 0)) # 10
    num_subforests=int(dic.get('numSubForests', 0)) #10
    nodata_value = int(dic.get('noDataValue', -9999)) # -9999
    run_local = dic.get('runLocal', False) # -9999
    run_local = False
    
    region_name = dic.get('regionOfInterest', None)
    gcs_mod_path  = dic.get('gcsModPath') #     gcs_mod_path = 'gs://rf_mods/'+region_name+'/'
    
    export_crs = dic.get('exportCrs', 'error')# 'EPSG:6350'
    export_scale = int(dic.get('exportScale', 0)) # 30
    export_descr = dic.get('descExp', '') 
    save_each_results = dic.get('saveIndividualResults', None)
    export_size = dic.get('exportSize', None)
    seed = dic.get('seed', 'noseed')
    
    
    # Create outdir if not available
    if not os.path.isdir(out_local_path):
        os.makedirs(out_local_path, exist_ok=True)
        
    #%% Load EE
# =============================================================================
#     if seed == 'noseed' or int(seed) is not int:
#         import random
#         seed = random.randint(1, 100000)
#         print ( 'Seed used is ', seed)
# =============================================================================
    #
    ## Provide local files
    if covs_local:
        print(' Using Local dataset')
        data = gpd.read_file( points_path ) # species  year  month 
        data_ee = geemap.geopandas_to_ee(data)
        # Have local dataset
        ## Run the things in Python+R
        #if run_local:
        #    data = df # Keep local files
        #else:
            # Convert local to ee dataset
            #print('    Converting local to EE')
            #
        #
    ## Provide EE dataset
    else: 
        print(' Using EE dataset')
        data_ee = ee.FeatureCollection( points_path )
        data = geemap.ee_to_gdf(data_ee)
        # Convert EE to local
        if run_local:
            #data = geemap.ee_to_gdf(df)
            print('    Converting EE to local')
        else:
            # Keep EE
            1
            #data_ee = df
        #
    #
    sim_abs = False
    if len(data[data["preabs"] == 0]):
        sim_abs = True
        #s.lower() in ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']

        
    #
    if remove_spat_dup:
        def remove_duplicates(data, pixel_size):
            # Select one occurrence record per pixel at the chosen spatial resolution
            random_raster = ee.Image.random().reproject('EPSG:4326', None, spat_dup_pixsize)
            rand_point_vals = random_raster.sampleRegions(
                collection=ee.FeatureCollection(data), geometries=True
            )
            return rand_point_vals.distinct('random')
        #
        print('  Removing spatial duplicates')
        print('   --  Original data size:', data_ee.size().getInfo())
        data_ee = remove_duplicates(data_ee, pixel_size)
        # Before selection and after selection
        print('   -- Final data size:', data_ee.size().getInfo())
    
    # =============================================================================
    # Definition of the Area of Interest
    # Defining the Area of Interest (AOI below) refers to the term used by researche
    # rs to denote the geographical area they want to analyze. It has a similar mean
    # ing to the term Study Area.
    # 
    # In this context, we obtained the bounding box of the occurrence point layer ge
    # ometry and created a 50-kilometer buffer around it (with a maximum tolerance o
    # f 1,000 meters) to define the AOI.
    # =============================================================================
    
    
    # Define the AOI
    if area_of_interest != "":
        # load aoi
        if asset_exists(area_of_interest):
            aoi = ee.FeatureCollection(area_of_interest).first()
        else:
            print(' AOI not found ', area_of_interest)
    else:
        # Create aoi
        aoi = data_ee.geometry().bounds().buffer(distance = points_buff_meters, maxError=1000)
        if  not run_local:
            1
  
    # Converting predictor values from Earth Engine to a DataFrame
    
    # =============================================================================
    # Calculating Spearman correlation coefficients between the given predictor vari
    # ables and visualizing them in a heatmap.
    # =============================================================================
    def plot_correlation_heatmap(dataframe, h_size=20, show_labels=False, path = '', showplots = False, idd = '_'):
        # Calculate Spearman correlation coefficients
        correlation_matrix = dataframe.corr(method='spearman')
        # Create a heatmap
        plt.figure(figsize=(h_size, h_size-2))
        plt.imshow(correlation_matrix, cmap='coolwarm', interpolation='nearest')
        # Optionally display values on the heatmap
        if show_labels:
            for i in range(correlation_matrix.shape[0]):
                for j in range(correlation_matrix.shape[1]):
                    plt.text(j, i, f'{correlation_matrix.iloc[i, j]:.2f}',
                             ha='center', va='center', color='white', fontsize=8)
        columns = dataframe.columns.tolist()
        plt.xticks(range(len(columns)), columns, rotation=90)
        plt.yticks(range(len(columns)), columns)
        plt.title('Variables Correlation Matrix')
        plt.colorbar(label='Spearman Correlation')
        if path != '':
            try:
                plt.savefig(path + '/correlation_heatmap_plot_' + idd + '.png')
            except:
                print(" Can't save the correlation heatmap plot")  
        if not showplots:
            plt.close()
            plt.ioff()
        else:
            plt.show()
    #
    #
    # Plot the correlation heatmap of variables
    #if run_local:
    plt_cor = plot_correlation_heatmap(dataframe=data.select_dtypes(exclude=['object','geometry', 'string']), 
                                       path = out_local_path)
    
    
    # =============================================================================
    # Spearman correlation coefficient is useful for understanding the general assoc
    # iations among predictor variables but does not directly assess how multiple va
    # riables interact, specifically detecting multicollinearity.
    # 
    # The Variance Inflation Factor (VIF below) is a statistical metric used to eval
    # uate multicollinearity and guide variable selection. It indicates the degree o
    # f linear relationship of each independent variable with the other independent 
    # variables, and high VIF values can be evidence of multicollinearity.
    # 
    # Typically, when VIF values exceed 5 or 10, it suggests that the variable has a
    #  strong correlation with other variables, potentially compromising the stabili
    # ty and interpretability of the model. In this tutorial, a criterion of VIF val
    # ues less than 10 was used for variable selection. The following 6 variables we
    # re selected based on VIF.
    # 
    # =============================================================================
    # Filter variables based on Variance Inflation Factor (VIF)
    def filter_variables_by_vif(dataframe, threshold=10):
        # dataframe = data
        ignore_columns = ['geometry', 'preabs', 'lc', 'Latitude', 'Longitude', 'date', 'wrong']
        df = dataframe.select_dtypes(exclude=['object', 'string'])
        df = df.drop(columns=ignore_columns, errors='ignore')
        original_columns = df.columns.tolist()
        remaining_columns = original_columns[:]
        while True:
            vif_data = dataframe[remaining_columns]
            vif_values = [
                variance_inflation_factor(vif_data.values, i)
                for i in range(vif_data.shape[1])
            ]
            max_vif_index = vif_values.index(max(vif_values))
            max_vif = max(vif_values)
            if max_vif < threshold:
                break
            print(f"Removing '{remaining_columns[max_vif_index]}' with VIF {max_vif:.2f}")
            del remaining_columns[max_vif_index]
        removed_columns = list(set(original_columns) - set(remaining_columns))
        #filtered_data = dataframe[remaining_columns]
        bands = remaining_columns # filtered_data.columns.tolist()
        print('Bands:', bands)
        return bands, removed_columns, ignore_columns
    
    #if run_local:
    # Variable Selection Based on VIF
    selected_columns, removed_columns, ignore_columns = filter_variables_by_vif(data)
    
    # Plot the correlation heatmap of variables
    p_corrh = plot_correlation_heatmap(data[selected_columns], h_size=6, show_labels=True)
    
    
    ## Create predictors stack
    anc_static_cov_img = anc.precooked_mosaic(stack = stack_type ) 
    # anc_static_cov_img.bandNames().getInfo()
    predictors = anc_static_cov_img.select( selected_columns )
    # predictors.bandNames().getInfo()
    # 
    # =============================================================================
    # Environmental Classification Using k-means Clustering: The k-means clustering 
    # algorithm, based on Euclidean distance, will be used to divide the pixels with
    # in the study area into two clusters. One cluster will represent areas with sim
    # ilar environmental characteristics to randomly selected 100 presence locations
    # , while the other cluster will represent areas with different characteristics.
    # 
    # Generation of Pseudo-Absence Data within Dissimilar Clusters: Within the secon
    # d cluster identified in the first step (which has different environmental char
    # acteristics from the presence data), randomly generated pseudo-absence points 
    # will be created. These pseudo-absence points will represent locations where th
    # e species is not expected to exist.
    # =============================================================================
    
    # Randomly select 100 locations for occurrence
    # Presence 
    #pvals = predictors.sampleRegions(
    #    collection=data.randomColumn().sort('random').limit(nabs),
    #    properties=[],
     #   scale=pixel_size
    #)
    sim_pseabs = False
    if sim_pseabs:
        pvals = data_ee.filterMetadata('preabs', 'equals', 1).limit(int(npres)).select(selected_columns)
        
        # Perform k-means clustering
        clusterer = ee.Clusterer.wekaKMeans(
            nClusters=int(nclus),
            distanceFunction=str(disfun)
        ).train(pvals)
        
        cl_result = predictors.cluster(clusterer)
        # Get cluster ID for locations similar to occurrence
        cl_id = cl_result.sampleRegions(
            collection=data_ee.randomColumn().sort('random').limit(200),
            properties=[],
            scale=finalmap 
        )
        # Define non-occurrence areas in dissimilar clusters
        cl_id = ee.FeatureCollection(cl_id).reduceColumns(ee.Reducer.mode(),['cluster'])
        cl_id = ee.Number(cl_id.get('mode')).subtract(1).abs()
        cl_mask = cl_result.select(['cluster']).eq(cl_id)
        # Presence location mask
        presence_mask = data_ee.reduceToImage(properties=['random'],
                                              reducer=ee.Reducer.first()
                                              ).reproject('EPSG:4326', None,
                                                          pixel_size).mask().neq(1).selfMask()
        # Masking presence locations in non-occurrence areas and clipping to AOI
        area_for_pa = presence_mask.updateMask(cl_mask).clip(aoi)

    # =============================================================================
    # Model fitting and prediction
    # We will now divide the data into training data and test data. The training dat
    # a will be used to find the optimal parameters by training the model, while the
    #  test data will be used to evaluate the model trained beforehand. An important
    #  concept to consider in this context is spatial autocorrelation.
    # 
    # Spatial autocorrelation is an essential element in SDM, associated with Tobler
    # 's law. It embodies the concept that 'everything is related to everything else
    # , but near things are more related than distant things'. Spatial autocorrelati
    # on represents the significant relationship between the location of species and
    #  environmental variables. However, if spatial autocorrelation exists between t
    # he training and test data, the independence between the two data sets can be c
    # ompromised. This significantly impacts the evaluation of the model's generaliz
    # ation ability.
    # 
    # One method to address this issue is the spatial block cross-validation techniq
    # ue, which involves dividing the data into training and testing datasets. This 
    # technique involves dividing the data into multiple blocks, using each block in
    # dependently as training and test datasets to reduce the impact of spatial auto
    # correlation. This enhances the independence between datasets, allowing for a m
    # ore accurate evaluation of the model's generalization ability.
    # 
    # The specific procedure is as follows:
    # 
    # Creation of spatial blocks: Divide the entire dataset into spatial blocks of e
    # qual size (e.g., 50x50 km).
    # Assignment of training and testing sets: Each spatial block is randomly assign
    # ed to either the training set (70%) or the test set (30%). This prevents the m
    # odel from overfitting to data from specific areas and aims to achieve more gen
    # eralized results.
    # Iterative cross-validation: The entire process is repeated n times (e.g., 10 t
    # imes). In each iteration, the blocks are randomly divided into training and te
    # st sets again, which is intended to improve the model's stability and reliabil
    # ity.
    # Generation of pseudo-absence data: In each iteration, pseudo-absence data are 
    # randomly generated to evaluate the model's performance.
    # 
    # =============================================================================
    #rand_grid_size = 100
    terrain = ee.Algorithms.Terrain(ee.Image('USGS/SRTMGL1_003'))
    watermask = terrain.select('elevation').gt(0)
    grid = watermask.reduceRegions(
        collection=aoi.coveringGrid(scale=int(rand_grid_size), proj='EPSG:4326'),
        reducer=ee.Reducer.mean()).filter(ee.Filter.neq('mean', None))
    #grid.size().getInfo()
    
      # =============================================================================
    # Now we can fit the model. Fitting a model involves understanding the patterns 
    # in the data and adjusting the model's parameters (weights and biases) accordin
    # gly. This process enables the model to make more accurate predictions when pre
    # sented with new data. For this purpose, we have defined a function called SDM(
    # ) to fit the model.
    # 
    # We will use the Random Forest algorithm.
    # =============================================================================
    
    presence_points = data_ee.filterMetadata('preabs', 'equals', 1).limit(int(npres)).select(selected_columns + ['preabs'])
    absence_points = data_ee.filterMetadata('preabs', 'equals', 0).limit(int(nabs)).select(selected_columns + ['preabs'] )
    #presence_points.first().getInfo() 
    # data_ee.size().getInfo() 
    # presence_points.size().getInfo() 
    # absence_points.size().getInfo()
    
    def sdm(x, save_ind = False, eelogpath = ''):
        # x = 305
        seed = ee.Number(x)
        # Random block division for training and validation
        rand_blk = ee.FeatureCollection(grid).randomColumn(seed=seed).sort('random')
        training_grid = rand_blk.filter(ee.Filter.lt('random', split))  # Grid for training
        testing_grid = rand_blk.filter(ee.Filter.gte('random', split))  # Grid for testing
        # Presence points
        #presence_points = ee.FeatureCollection(data)
        #presence_points = presence_points.map(lambda feature: feature.set('PresAbs', 1))

        tr_presence_points = presence_points.filter( ee.Filter.bounds(training_grid) )  # Presence points for training
        te_presence_points = presence_points.filter( ee.Filter.bounds(testing_grid) )  # Presence points for testing
        # tr_presence_points.size().getInfo() 
        # te_presence_points.size().getInfo()
        
        # Pseudo-absence points for training
        # tr_pseudo_abs_points = area_for_pa.sample(
        #    region=training_grid, scale=pixel_size,
        #    numPixels=tr_presence_points.size().add(300),
        #    seed=seed, geometries=True)
        tr_pseudo_abs_points = absence_points.filter( ee.Filter.bounds(training_grid) )
        te_pseudo_abs_points = absence_points.filter( ee.Filter.bounds(testing_grid) )
        # tr_pseudo_abs_points.size().getInfo() 
        # te_pseudo_abs_points.size().getInfo()
        
        # Same number of pseudo-absence points as presence points for training
        #tr_pseudo_abs_points = (tr_pseudo_abs_points.randomColumn().sort('random').limit(ee.Number(tr_presence_points.size()))        )
        #tr_pseudo_abs_points = tr_pseudo_abs_points.map(lambda feature: feature.set('PresAbs', 0))
        
        #te_pseudo_abs_points = area_for_pa.sample(region=testing_grid, scale=pixel_size,
         #  numPixels=te_presence_points.size().add(100),seed=seed,geometries=True )
        # Same number of pseudo-absence points as presence points for testing
        #te_pseudo_abs_points = ( te_pseudo_abs_points.randomColumn()
         #   .sort('random').limit(ee.Number(te_presence_points.size())))
        #te_pseudo_abs_points = te_pseudo_abs_points.map(lambda feature: feature.set('PresAbs', 0))
        # Merge training and pseudo-absence points
        
        training_partition = tr_presence_points.merge(tr_pseudo_abs_points)
        testing_partition = te_presence_points.merge(te_pseudo_abs_points)
        # training_partition.size().getInfo()
        # testing_partition.size().getInfo()
        
        # Extract predictor variable values at training points
        #train_pvals = predictors.sampleRegions(
        #    collection=training_partition, properties=['PresAbs'],
        #    scale=pixel_size, geometries=True,
        #)
        
        # Random Forest classifier
        classifier = ee.Classifier.smileRandomForest(
            numberOfTrees=500,
            variablesPerSplit=None,
            minLeafPopulation=3,
            bagFraction=0.5,
            maxNodes=None,
            seed=seed,
        )
        # Presence probability: Habitat suitability map
        classifier_pr = classifier.setOutputMode('PROBABILITY').train(
            training_partition, 'preabs', selected_columns
        )
        classified_img_pr = predictors.classify(classifier_pr)
        
        if save_ind:
            taskA = ee.batch.Export.image.toAsset(
                image = classified_img_pr ,
                assetId = str(out_ee_path + '/export_classified_img_pr'+str(x)).replace('//', '/'),
                description = 'classified_img_pr'+str(x) ,
                scale = export_scale,  region = grid.geometry(), #crs=export_crs,
                maxPixels=1e13)
            taskA.start()
        
        # Binary presence/absence map: Potential distribution map
        classifier_bin = classifier.setOutputMode('CLASSIFICATION').train(
            training_partition, 'preabs', selected_columns
        )
        classified_img_bin = predictors.classify(classifier_bin)
        
        if save_ind:        
            taskB = ee.batch.Export.image.toAsset(
                image = classified_img_bin ,
                assetId = str(out_ee_path + '/export_classified_img_bin_'+str(x)).replace('//', '/'),
                description = 'classified_img_bin_'+str(x) ,
                scale = export_scale,  region = grid.geometry(), #crs=export_crs,
                maxPixels=1e13)
            taskB.start()
        
        answer = [classified_img_pr, classified_img_bin, training_partition, testing_partition], classifier_pr
        
        return answer
    #
    #
    #
    #
    
    # =============================================================================
    # Spatial blocks are divided into 70% for model training and 30% for model testi
    # ng, respectively. Pseudo-absence data are randomly generated within each train
    # ing and testing set in every iteration. As a result, each execution yields dif
    # ferent sets of presence and pseudo-absence data for model training and testing
    # .
    # =============================================================================
    
    # Random Seed
    runif = lambda length: [random.randint(1, 1000) for _ in range(length)]
    items = runif(numiter)
    # Fixed seed
    # items = [287, 288, 553, 226, 151, 255, 902, 267, 419, 538]
    results_list = [] # Initialize SDM results list
    importances_list = [] # Initialize variable importance list
    #
    for item in items:
        print( item )
        result, partition, trained = sdm(item)
        # result, partition, trained = answer
        ## Accumulate SDM results into the list
        results_list.extend(result)
        # Accumulate variable importance into the list
        importance = ee.Dictionary(trained.explain()).get('importance')
        importances_list.extend(importance.getInfo().items())
    
    # Flatten the SDM results list
    results = ee.List(results_list).flatten()
    # gir = results.getInfo()
    # =============================================================================
    # Now we can visualize the habitat suitability map and potential distribution ma
    # p for the Fairy pitta. In this case, the habitat suitability map is created by
    #  using the mean() function to calculate the average for each pixel location ac
    # ross all images, and the potential distribution map is generated by using the 
    # mode() function to determine the most frequently occurring value at each pixel
    #  location across all images.
    # =============================================================================
    # =============================================================================
    # SDM results
    # =============================================================================
    # Habitat suitability map
    images = ee.List.sequence(
        0, ee.Number(numiter).multiply(4).subtract(1), 4).map(
        lambda x: results.get(x))
    
    model_average = ee.ImageCollection.fromImages(images).mean()
    

    # Potential distribution map
    images2 = ee.List.sequence(1, ee.Number(numiter).multiply(4).subtract(1), 4).map(
        lambda x: results.get(x)
    )
    distribution_map = ee.ImageCollection.fromImages(images2).mode().rename('binary')
    distribution_map.id
    
    finalmap = model_average.rename('average') 
    finalmap = finalmap.addBands(distribution_map )
    # xgi = finalmap.getInfo() 
    
    taskk = ee.batch.Export.image.toAsset(image = distribution_map ,
                                          assetId = str(out_ee_path + '/export_finalmap').replace('//', '/'),
                                          description = export_descr ,
                                          scale = export_scale, 
                                          region = grid.geometry(),
                                          #crs=export_crs,
                                          maxPixels=1e13)
    taskk.start()
     
   
    # =============================================================================
    # Variable importance and accuracy assessment
    # Random Forest (ee.Classifier.smileRandomForest) is one of the ensemble learnin
    # g methods, which operates by constructing multiple decision trees to make pred
    # ictions. Each decision tree independently learns from different subsets of the
    #  data, and their results are aggregated to enable more accurate and stable pre
    # dictions.
    #
    # Variable importance is a measure that evaluates the impact of each variable on
    #  the predictions within the Random Forest model. We will use the previously de
    # fined importances_list to calculate and print the average variable importance.
    # =============================================================================
    def plot_variable_importance(importances_list, path = '', showplots = False):
        # Extract each variable importance value into a list
        variables = [item[0] for item in importances_list]
        importances = [item[1] for item in importances_list]
        # Calculate the average importance for each variable
        average_importances = {}
        for variable in set(variables):
            indices = [i for i, var in enumerate(variables) if var == variable]
            average_importance = np.mean([importances[i] for i in indices])
            average_importances[variable] = average_importance
        # Sort the importances in descending order of importance
        sorted_importances = sorted(average_importances.items(),
                                    key=lambda x: x[1], reverse=False)
        variables = [item[0] for item in sorted_importances]
        avg_importances = [item[1] for item in sorted_importances]
        # Adjust the graph size
        plt.figure(figsize=(8, 4))
        # Plot the average importance as a horizontal bar chart
        plt.barh(variables, avg_importances)
        plt.xlabel('Importance')
        plt.ylabel('Variables')
        plt.title('Average Variable Importance')
        # Display values above the bars
        for i, v in enumerate(avg_importances):
            plt.text(v + 0.02, i, f'{v:.2f}', va='center')
        # Adjust the x-axis range
        plt.xlim(0, max(avg_importances) + 5)  # Adjust to the desired range
        plt.tight_layout()
        plt.savefig(path + '/variable_importance.png')
        if not showplots:
            plt.close()
            plt.ioff()
        else:
            plt.show()
    #    
    plot_variable_importance(importances_list, path=out_local_path, showplots=showplots)
    # =============================================================================
    # Using the Testing Datasets, we calculate AUC-ROC and AUC-PR for each run. Then
    # , we compute the average AUC-ROC and AUC-PR over n iterations.
    # 
    # AUC-ROC represents the area under the curve of the 'Sensitivity (Recall) vs. 1
    # -Specificity' graph, illustrating the relationship between sensitivity and spe
    # cificity as the threshold changes. Specificity is based on all observed non-oc
    # currences. Therefore, AUC-ROC encompasses all quadrants of the confusion matri
    # x.
    # 
    # AUC-PR represents the area under the curve of the 'Precision vs. Recall (Sensi
    # tivity)' graph, showing the relationship between precision and recall as the t
    # hreshold varies. Precision is based on all predicted occurrences. Hence, AUC-P
    # R does not include the true negatives (TN).
    # 
    # Note: It's important to ensure that each run has a sufficient number of points
    #  for model validation. The final number of points may vary due to the random p
    # artitioning of spatial blocks, so it's crucial to verify if there are enough p
    # resence and pseudo-absence points for model validation. In the case of endange
    # red or rare species, there might be a shortage of occurrence data, leading to 
    # an insufficient test dataset. In such cases, alternatives may include addition
    # al data collection based on expert knowledge and experience or utilizing relev
    # ant alternative data sources.
    # =============================================================================
    def print_pres_abs_sizes(TestingDatasets, numiter):
        # Check and print the sizes of presence and pseudo-absence coordinates
        def get_pres_abs_size(x):
            fc = ee.FeatureCollection(TestingDatasets.get(x))
            presence_size = fc.filter(ee.Filter.eq('PresAbs', 1)).size()
            pseudo_absence_size = fc.filter(ee.Filter.eq('PresAbs', 0)).size()
            return ee.List([presence_size, pseudo_absence_size])
        sizes_info = (
            ee.List.sequence(0, ee.Number(numiter).subtract(1), 1)
            .map(get_pres_abs_size)
            .getInfo()
        )
        for i, sizes in enumerate(sizes_info):
            presence_size = sizes[0]
            pseudo_absence_size = sizes[1]
            print(
                f'Iteration {i + 1}: Presence Size = {presence_size}, Pseudo-absence Size = {pseudo_absence_size}'
            )
    # Extracting the Testing Datasets
    testing_datasets = ee.List.sequence(
        3, ee.Number(numiter).multiply(4).subtract(1), 4
    ).map(lambda x: results.get(x))
    
    print_pres_abs_sizes(testing_datasets, numiter)
    #
    def get_acc(hsm, t_data, pixel_size):
        pr_prob_vals = hsm.sampleRegions(
            collection=t_data, properties=['PresAbs'], scale=pixel_size
        )
        seq = ee.List.sequence(start=0, end=1, count=25)  # Divide 0 to 1 into 25 intervals
        def calculate_metrics(cutoff):
            # Each element of the seq list is passed as cutoff(threshold value)
            # Observed present = TP + FN
            pres = pr_prob_vals.filterMetadata('PresAbs', 'equals', 1)
            # TP (True Positive)
            tp = ee.Number(
                pres.filterMetadata('classification', 'greater_than', cutoff).size()
            )
            # TPR (True Positive Rate) = Recall = Sensitivity = TP / (TP + FN) = TP / Observed present
            tpr = tp.divide(pres.size())
            # Observed absent = FP + TN
            abs = pr_prob_vals.filterMetadata('PresAbs', 'equals', 0)
            # FN (False Negative)
            fn = ee.Number(
                pres.filterMetadata('classification', 'less_than', cutoff).size()
            )
            # TNR (True Negative Rate) = Specificity = TN  / (FP + TN) = TN / Observed absent
            tn = ee.Number(abs.filterMetadata('classification', 'less_than', cutoff).size())
            tnr = tn.divide(abs.size())
            # FP (False Positive)
            fp = ee.Number(
                abs.filterMetadata('classification', 'greater_than', cutoff).size()
            )
            # FPR (False Positive Rate) = FP / (FP + TN) = FP / Observed absent
            fpr = fp.divide(abs.size())
            # Precision = TP / (TP + FP) = TP / Predicted present
            precision = tp.divide(tp.add(fp))
            # SUMSS = SUM of Sensitivity and Specificity
            sumss = tpr.add(tnr)
            return ee.Feature(
                None,
                {
                    'cutoff': cutoff,
                    'TP': tp,
                    'TN': tn,
                    'FP': fp,
                    'FN': fn,
                    'TPR': tpr,
                    'TNR': tnr,
                    'FPR': fpr,
                    'Precision': precision,
                    'SUMSS': sumss,
                },
            )
        return ee.FeatureCollection(seq.map(calculate_metrics))
        #
    def calculate_and_print_auc_metrics(images, testing_datasets, pixel_size, numiter, path):
        # Calculate AUC-ROC and AUC-PR
        def calculate_auc_metrics(x):
            hsm = ee.Image(images.get(x))
            t_data = ee.FeatureCollection(testing_datasets.get(x))
            acc = get_acc(hsm, t_data, pixel_size)
            # Calculate AUC-ROC
            x = ee.Array(acc.aggregate_array('FPR'))
            y = ee.Array(acc.aggregate_array('TPR'))
            x1 = x.slice(0, 1).subtract(x.slice(0, 0, -1))
            y1 = y.slice(0, 1).add(y.slice(0, 0, -1))
            auc_roc = x1.multiply(y1).multiply(0.5).reduce('sum', [0]).abs().toList().get(0)
            # Calculate AUC-PR
            x = ee.Array(acc.aggregate_array('TPR'))
            y = ee.Array(acc.aggregate_array('Precision'))
            x1 = x.slice(0, 1).subtract(x.slice(0, 0, -1))
            y1 = y.slice(0, 1).add(y.slice(0, 0, -1))
            auc_pr = x1.multiply(y1).multiply(0.5).reduce('sum', [0]).abs().toList().get(0)
            return (auc_roc, auc_pr)
        auc_metrics = (
            ee.List.sequence(0, ee.Number(numiter).subtract(1), 1)
            .map(calculate_auc_metrics)
            .getInfo()
        )
        # Print AUC-ROC and AUC-PR for each iteration
        df = pd.DataFrame(auc_metrics, columns=['AUC-ROC', 'AUC-PR'])
        df.index = [f'Iteration {i + 1}' for i in range(len(df))]
        df.to_csv(path + 'auc_metrics.csv', index_label='Iteration')
        print(df)
        # Calculate mean and standard deviation of AUC-ROC and AUC-PR
        mean_auc_roc, std_auc_roc = df['AUC-ROC'].mean(), df['AUC-ROC'].std()
        mean_auc_pr, std_auc_pr = df['AUC-PR'].mean(), df['AUC-PR'].std()
        print(f'Mean AUC-ROC = {mean_auc_roc:.4f} ± {std_auc_roc:.4f}')
        print(f'Mean AUC-PR = {mean_auc_pr:.4f} ± {std_auc_pr:.4f}')
    
    #%%time
    # Calculate AUC-ROC and AUC-PR
    metr = calculate_and_print_auc_metrics(images, testing_datasets, pixel_size, numiter, path = path)
    metr
# End function
#%%
if __name__ == "__main__":
    main()


# =============================================================================
# This tutorial has provided a practical example of using Google Earth Engine (G
# EE) for Species Distribution Modeling (SDM). An important takeaway is the vers
# atility and flexibility of GEE in the field of SDM. Leveraging Earth Engine's 
# powerful geospatial data processing capabilities opens up endless possibilitie
# s for researchers and conservationists to understand and preserve biodiversity
#  on our planet. By applying the knowledge and skills gained from this tutorial
# , individuals can explore and contribute to this fascinating field of ecologic
# al research.
# 
# Was this helpful?
# 
# Send feedback
# Except as otherwise noted, the content of this page is licensed under the Crea
# tive Commons Attribution 4.0 License, and code samples are licensed under the 
# Apache 2.0 License. For details, see the Google Developers Site Policies. Java
#  is a registered trademark of Oracle and/or its affiliates.
# Last updated 2026-04-19 UTC.
# =============================================================================
