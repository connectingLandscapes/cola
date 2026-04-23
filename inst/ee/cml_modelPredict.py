# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 21:00:16 2023
Cinnects EE from command line
@author: ig299
"""

#%%
# IMPORTS
import ee
import sys
import csv
import os 

# os.chdir('N:/My Drive/git/cola/inst/ee')
#/home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_covariates_extraction.py /home/shiny/sat_ts_fusion/params_ext_cov.csv
# C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_modelPredict.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv"
def main() -> None:
    #%%

    # INPUTS
    paramcsv = sys.argv[1] # project name
    # paramcsv = "N:/My Drive/git/cola/inst/ee/params.csv"
    with open(paramcsv) as f:
        dic = dict(filter(None, csv.reader(f)))
    dic
    #
    eeproject = dic['eeProject']
    #
    try:
        ee.Authenticate( )
        ee.Initialize(project = eeproject)
        print('EE initialized')

    except:
        print(' ERROR: no EE initialized')
        exit(1)
    
    #########################################333
    from sat_ts_fusion.imagery import lsat_utils
    from sat_ts_fusion.ccdc import ccdc_utils
    from sat_ts_fusion.fusion import ancillary_covariates
    from sat_ts_fusion.fusion import covariates as covs
    # from sat_ts_fusion.fusion.models import models


    gedi_metrics = dic.get('gediMetrics')#['agbd_a0', 'rh_98_a0']
    gedi_metric_int16_mult = [int(x) for x in dic.get('gediMetricsInt16Mult', '0,0').strip().split(',')] #[10, 100]
    projection_years = [int(x) for x in dic.get('projectionYears', '0,0').strip().split(',')]#[2019,2020,2021,2022]
    num_folds = int(dic.get('numFolds', '0')) # 10
    num_subforests=int(dic.get('numSubForests', '0')) #10
    annual_date_frac = float(dic.get('annualDateFrac', '0')) #.5
    export_scale = int(dic.get('exportScale', '0')) # 30
    export_crs = dic.get('exportCrs', 'error')# 'EPSG:6350'
    nodata_value = int(dic.get('noDataValue', -9999)) # -9999
    
    region_name = dic.get('regionOfInterest')
    gcs_mod_path  = dic.get('gcsModPath') #     gcs_mod_path = 'gs://rf_mods/'+region_name+'/'
    
    landcover_ic = dic.get('landcoverIc')
    pts_path = dic.get('ptsWithCovs')
    
    eeFolderPredic = dic.get('eeFolderPredic')
    saveIndividualResults = dic.get('saveIndividualResults')
    seed = dic.get('seed', 'noseed')
    
    if seed == 'noseed' or int(seed) is not int:
        import random
        seed = random.randint(1, 100000)
        print ( 'Seed used is ', seed)
    #
    
    landcoverKeepClasses = [int(x) for x in dic.get('landcoverKeepClasses').strip().split(',')]
    landcoverPhenoVariesClasses = [int(x) for x in dic.get('landcoverPhenoVariesClasses').strip().split(',')]    
    fromre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
    tore = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 10, 0, 0, 0, 0, 0, 11, 11, 11, 11, 12, 255]
    landcoverKeepClasses =  list(set(tore))
    landcoverPhenoVariesClasses = landcoverKeepClasses
    
    
    lcic = ee.Image( landcover_ic ) #   
    landcover_ic = lcic.remap(from_= fromre, to=tore, bandName='b1',  defaultValue = None ).rename('lc')
    #
    aoi = ee.FeatureCollection(region_name)
    pts = ee.FeatureCollection(pts_path)
    ffirst = pts.first(); 
    figi = ffirst.getInfo() 
    type(figi)
    
    anc_static_cov_img = ancillary_covariates.make_covar_topo_stack()
    predictors = anc_static_cov_img.addBands(landcover_ic).clip(aoi)
    
    task2 = ee.batch.Export.image.toAsset(image = sui2.clip(aoi),
                                         assetId = eeFolderPredic + '/seed' + str(seed) + 'bin',
                                         description = ' Export binary ' + str(seed) ,
                                         scale=export_scale, 
                                         region=aoi.geometry(),
                                         #crs=export_crs,
                                         maxPixels=1e13)
    
    bands = predictors.bandNames()
    # bands.getInfo()
    #
    #
    ######
    def sdm(x):
        #seed = ee.Number(x)
        # seed = 287
        
        fctosplit = pts.randomColumn( columnName = 'random', seed = seed, distribution = 'uniform')
        fcgi = fctosplit.first().getInfo() 
        
        training = fctosplit.filter(ee.Filter.lt('random', 0.7))
        testing = fctosplit.filter(ee.Filter.gt('random', 0.7))
        
        # pointssize = pts.size().getInfo()
        # training.size().getInfo() # 62
        # testing.size().getInfo() # 24
    
        # Random Forest classifier
        classifier = ee.Classifier.smileRandomForest(
            numberOfTrees=500,
            variablesPerSplit=None,
            minLeafPopulation=10,
            bagFraction=0.5,
            maxNodes=None,
            seed=seed,
        )
        
        # Presence probability: Habitat suitability map
        classifier_pr = classifier.setOutputMode("PROBABILITY").train(
            training, "preabs", bands
        )
        
        #
        classified_img_pr = predictors.select(bands).classify(classifier_pr)
    
        # Binary presence/absence map: Potential distribution map
        classifier_bin = classifier.setOutputMode("CLASSIFICATION").train(
            training, "preabs", bands)
        
        classified_img_bin = predictors.select(bands).classify(classifier_bin)
        
        return classified_img_pr, classified_img_bin, classifier_pr
            # , training, testing 
    #
    siu1, sui2, alg = sdm(seed)
    
    try:
        task1 = ee.batch.Export.image.toAsset(image = siu1.clip(aoi),
                                             assetId = eeFolderPredic + '/seed' + str(seed) + 'cont',
                                             description = ' Export continuos ' + str(seed) ,
                                             scale=export_scale, 
                                             region=aoi.geometry(),
                                             #crs=export_crs,
                                             maxPixels=1e13)
        #
        task2 = ee.batch.Export.image.toAsset(image = sui2.clip(aoi),
                                             assetId = eeFolderPredic + '/seed' + str(seed) + 'bin',
                                             description = ' Export binary ' + str(seed) ,
                                             scale=export_scale, 
                                             region=aoi.geometry(),
                                             #crs=export_crs,
                                             maxPixels=1e13)
                                             # skipEmptyTiles=True,fileFormat='GeoTIFF')
        task1.start()
        task2.start()
        #
        #task3 = ee.batch.Export.classifier.toAsset(alg, assetId = eeFolderPredic + '/seed' + str(seed) + 'class',
        #                                     description = ' Export classifier ' + str(seed))
        # task3.start()
        print('\n Cola2 EE Tasks submited \n') 
    except ValueError as e:
        print('Cola 2 ERROR submitting task: ', e)
        exit(1)
        
    
 
# =============================================================================
#     # TODO: predictors are renamed in R script to save space in text file; need a lookup table with covariate names and shortened IDs
#     cov_numbers_0 = list(range(64))
#     old_cov_names = []
#     for n in cov_numbers_0:
#         old_cov_names.append(f'gse_a{n:02d}')
# 
#     cov_numbers_1 = list(range(1,65))
#     new_cov_names = []
#     for n in cov_numbers_1:
#         new_cov_names.append(f'p{n}')
# 
#     cov_img_rn = cov_img.select(old_cov_names, new_cov_names)
# 
#     # clip and mask the covariate image
#     cov_img_cm = models.clip_mask_cov_img(cov_img=cov_img_rn,
#                                           year=year,
#                                           region_geom=region_geom,
#                                           landcover_ic=ee.ImageCollection(landcover_ic),
#                                           landcover_keep_classes=landcoverKeepClasses)
# 
#     # project the model
#     projected_mean_sd = models.project_rf_text_mean_sd(gcs_mod_path=gcs_mod_path,
#                                                        num_subforests=num_subforests,
#                                                        metric=metric,
#                                                        int16_metric_mult=gedi_metric_int16_mult[i],
#                                                        cov_img=cov_img_cm)
# 
# 
#     task = ee.batch.Export.image.toDrive(image=projected_mean_sd.unmask(nodata_value),
#                                          description=region_metric_year_name+' to Google Drive',
#                                          fileNamePrefix=region_metric_year_name,
#                                          scale=export_scale,
#                                          region=region_geom,
#                                          crs=export_crs,
#                                          maxPixels=1e13,
#                                          skipEmptyTiles=True,
#                                          fileFormat='GeoTIFF',
#                                          formatOptions={'cloudOptimized': True,
#                                                         'noData': nodata_value})
# 
#     task.start()
#     
#     i = i+1
# =============================================================================
    ############################################### ----
# End function




if __name__ == "__main__":
    main()

