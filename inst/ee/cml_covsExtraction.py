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
# C:\\Users\\gonza\\AppData\\Local\\r-miniconda\\envs\\cola\\python.exe C:\\cola\\cola2\\cml_covariates_extraction.py C:\\cola\\cola2\\params_ext_cov.csv
#/home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_covariates_extraction.py /home/shiny/sat_ts_fusion/params_ext_cov.csv
# C:/Users/ig299/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_covsExtraction.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv" 
# C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe "N:/My Drive/git/cola/inst/ee/cml_covsExtraction.py" "N:/My Drive/git/cola/inst/ee/params_extcovs.csv"

def main() -> None:
    #%% Params

    # INPUTS
    # paramcsv to file holding xy coordinates
    paramcsv = sys.argv[1] # project name
    # paramcsv = "N:/My Drive/git/cola/inst/ee/params.csv"
    try:
        eelogpath = str(sys.argv[2]) # 
        if not os.path.isdir(eelogpath):
            os.makedirs(eelogpath, exist_ok=True)
        if not os.path.isdir(eelogpath):
            eelogpath = ""
    except IndexError:
        eelogpath = ""
    # eelogpath = 'C:/Users/ig299/cola'
        
    
    with open(paramcsv) as f:
        print(f)
        dic = dict(filter(None, csv.reader(f)))
    #  
    print(dic)
    
    eeproject = dic['eeProject']

    try:
        ee.Authenticate()
        ee.Initialize(project = eeproject)
        print('\n Cola2: EE initialized')
    except ValueError as e:
        print('Cola 2 ERROR: ', e)
        exit(1)
    
    # param1 = 'gonzalezivan'
    # param2 = 'projects/gonzalezivan/assets/cola/anoamanual'
    # param3 = 10000
    # param4 = 'projects/gonzalezivan/cola/sim1'
    
    
    #########################################
    # os.chdir('N:/My Drive/git/cola/inst/ee')
    from sat_ts_fusion.fusion import ancillary_covariates
    from sat_ts_fusion.imagery import lsat_utils
    from sat_ts_fusion.ccdc import ccdc_utils
    from sat_ts_fusion.fusion import covariates as covs
    from sat_ts_fusion import eecolatools as ect


    # point feature collection to use for covariate extraction
    # needs to have a 'date' field with dates like '2020-06-25'
    points_path = dic['points']
    pt_fc_id = points
    
    # maximum number of points to extract covariates
    # ext_limit = 200000
    ext_limit = int(dic['maxNumbPoints'])
    
    # export asset ID path
    
    pt_fc_ext_id = dic['outExtCovs']
    #pt_fc_ext_id = 'projects/ee-gedibio/assets/gediv002_l2l4a_20190417to20230316_covext_'+region_name
    #pt_fc_ext_id = 'projects/gonzalezivan/assets/cola/a_covext'
    
    do_ccdc = False
    
    if do_ccdc:  
        # path to CCDC results
        # set to None to exclude CCDC covariates
        ccdc_result_path = dic['ccdc_result_path']
        #'projects/ee-gedibio/assets/ccdc/results/workshop_sonoma/'
        # CCDC start year
        ccdc_start_year = int(dic['ccdcStartYear']) #2000
        # CCDC bands to extract segments coefficients for
        ccdc_seg_bands = dic['ccdcSegBands'].strip().split(',') # ['NBR2', 'NDMI', 'NDVI']
        # CCDC bands to extract synthetic (i.e. fitted) values for
        ccdc_syn_bands = dic['ccdcSynBands'].strip().split(',') # ['NBR2', 'NDMI', 'NDVI']
        #ccdc_syn_bands = ['NIR', 'SWIR1', 'SWIR2']
        # CCDC bands to extract time since break values for
        ccdc_tsin_bands = dic['ccdcTsinBands'].strip().split(',') 
        # break threshold values associated with ccdc_tsin_bands
        ccdc_tsin_thresh = [float(x) for x in dic['ccdcTsinThres'].strip().split(',')]
        #ccdc_tsin_thresh = [0.1, 0.25]
        #
        if ccdc_result_path == 'None':
            print( 'ccdc_result_path  is None')
            ccdc_result_path = None
        #
        #
        ## make the CCDC mosaic image
        if ccdc_result_path is not None:
            ccdc_img = ccdc_utils.mosaic_ccdc_img_tiles(ccdc_result_path)
        else:
            ccdc_img = None
        #
        # ccdc_img = ccdc_utils.mosaic_ccdc_img_tiles(ccdc_result_path)
        print('  ++ ccdc_img:')
        print(ccdc_img)
    
    
    stack_type = dic['stackType'] 
    landcover_path  = dic['landcoverIc']
    landcoverKeepClasses = [int(x) for x in dic['landcoverKeepClasses'].strip().split(',')]
    landcoverPhenoVariesClasses = [int(x) for x in dic['landcoverPhenoVariesClasses'].strip().split(',')]
    include_gse_covs = dic['includeGseCovs'].replace('False', '')
    ext_start_year = int(dic['extStartYear'])
    ext_end_year = int(dic['extEndYear'])
    ext_pheno_start_doy = int(dic['extPhenoStartDoy'])
    ext_pheno_end_doy = int(dic['extPhenoEndDoy'])
    ext_scale = int(dic['extScale'])
    temporal_match = dic['temporalMatch']
    annual_date_frac = float(dic['annualDateFrac'])
    absences = int(dic['absencesToSimulate'])
    roi = dic['regionOfInterest']
    
    ###############
    #%% PROCESSING
    ###############
    
    # load the point FeatureCollection, and limit collection size
    pt_fc = ee.FeatureCollection(pt_fc_id).limit(maximum=ext_limit, prop='random_1')
    columns = pt_fc.first().propertyNames().getInfo()
    #
    if 'preabs' not in columns:
        print (  ' Creating the missing "preabs" column, assigning 1 to all features')        
        def add_preabs(feature):
            return feature.set('preabs', ee.Number(1))
        #
        pt_fc = pt_fc.map(add_preabs)
        #
    #
    if 'date' not in columns:
        print (  ' Creating the missing "date" column, assigning 2000-01-01 to all features')        
        def add_missingdate(feature):
            # feature = pt_fc.first()
            fea = feature.set('date', ee.String('2000-01-01')) # fea.getInfo()
            return fea 
        #
        pt_fc = pt_fc.map(add_preabs)
        #
    #
    # make the ancillary covariate image (only elev, slope, and aspect in this case)
    anc_static_cov_img = precooked_mosaic(stack = 1)
    #
    # Remap clc
    # Original values
    fromre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
    tore = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 10, 0, 0, 0, 0, 0, 11, 11, 11, 11, 12, 255]
    lcic = ee.Image( landcover_path ) # 
    
    landcoverKeepClasses =  list(set(tore))
    landcoverPhenoVariesClasses = landcoverKeepClasses
    
    landcover_ic = lcic.rename('lc')#.remap(from_= fromre, to=tore, bandName='b1',  defaultValue = None ).rename('lc')
    #landcover_ic.getInfo()
    ## add absences
    if absences != 0:
        # Generate 5,000 random points
        eeaoi = ee.FeatureCollection(roi)
        abs_pts = landcover_ic.sample(region=eeaoi, scale=30, numPixels=absences, geometries=True)
        # Extract predictor variable values
        def add_preabs0(feature):
            return feature.set('preabs', ee.Number(0))
        #
        abs_pts = pt_fc.map(add_preabs0)
        pt_fc = pt_fc.merge(abs_pts)
    #
    # giptsfc = pt_fc.getInfo()      # test
    #
    #
    # Run the covariate extraction
    # fc=pt_fc; ancillary_static_cov_img=anc_static_cov_img; landcover_keep_classes = landcoverKeepClasses; landcover_pheno_varies_classes=landcoverPhenoVariesClasses; ext_start_doy = 1; ext_end_doy = 366; 
    fc_ext = covs.extract_covariates_at_points2(fc=pt_fc,
                                                ancillary_static_cov_img=anc_static_cov_img,
                                                landcover_ic = landcover_ic,
                                                landcover_keep_classes = landcoverKeepClasses,
                                                landcover_pheno_varies_classes=landcoverPhenoVariesClasses,
                                                include_gse_covs=include_gse_covs,
                                                ext_start_year=ext_start_year,
                                                ext_end_year=ext_end_year,
                                                ext_pheno_start_doy=ext_pheno_start_doy,
                                                ext_pheno_end_doy=ext_pheno_end_doy,
                                                ext_scale=ext_scale,
                                                annual_date_frac=annual_date_frac)
                                                    
    # =============================================================================
    #                                               temporal_match=temporal_match,
    #                                               ccdc_img=ccdc_img,
    #                                               ccdc_start_year=ccdc_start_year,
    #                                               ccdc_seg_bands=ccdc_seg_bands,
    #                                               ccdc_syn_bands=ccdc_syn_bands,
    #                                               ccdc_tsin_bands=ccdc_tsin_bands,
    #                                               ccdc_tsin_thresh=ccdc_tsin_thresh,
    # =============================================================================
    
    # export the feature collection with covariates
    # testtask = ee.batch.Export.table.toAsset(collection=fc_ext, description='test', assetId='projects/gonzalezivan/assets/cola/test'); testtask.start()
    
    try:
        task = ee.batch.Export.table.toAsset(collection=fc_ext,
                                             description='Cola2 '+ os.path.basename(points_path)+' FC cov extract',
                                             assetId=pt_fc_ext_id)
        task.start()
        print('\n Cola2 EE Task submited \n', task)
        if not eelogpath == "" :
            ect.writeTask(task, eelogpath)
    except ValueError as e:
        print('Cola 2 ERROR submitting task: ', e)
        exit(1)
    

# End function

if __name__ == "__main__":
    main()
