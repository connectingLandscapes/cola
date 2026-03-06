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

#/home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/sat_ts_fusion/cml_covariates_extraction.py /home/shiny/sat_ts_fusion/params_ext_cov.csv
def main() -> None:
    #%%

    # INPUTS
    # paramcsv to file holding xy coordinates
    paramcsv = sys.argv[1] # project name
    #paramcsv = '/home/shiny/sat_ts_fusion/params_ext_cov.csv'
    
    with open(paramcsv) as f:
      dic = dict(filter(None, csv.reader(f)))
      
    print(dic)
    
    eeproject = dic['eeProject']
    
    try:
        ee.Authenticate()
        ee.Initialize(project = eeproject)
        print('EE initialized')

    except:
        print(' ERROR: no EE initialized')
        exit(1)
    
    # param1 = 'gonzalezivan'
    # param2 = 'projects/gonzalezivan/assets/cola/anoamanual'
    # param3 = 10000
    # param4 = 'projects/gonzalezivan/cola/sim1'
    
    
    #########################################333
    from sat_ts_fusion.imagery import lsat_utils
    from sat_ts_fusion.ccdc import ccdc_utils
    from sat_ts_fusion.fusion import ancillary_covariates
    from sat_ts_fusion.fusion import covariates as covs
    from sat_ts_fusion.fusion.models import models

    # export asset ID path
    #pt_fc_ext_id = 'projects/ee-gedibio/assets/gediv002_l2l4a_20190417to20230316_covext_'+region_name
    

    #//////  dic['']
    region_name = dic['regionName']
    gedi_metrics = dic['gediMetrics']#['agbd_a0', 'rh_98_a0']
    gedi_metric_int16_mult = [int(x) for x in dic['gediMetricsInt16Mult'].strip().split(',')] #[10, 100]
    projection_years = [int(x) for x in dic['projectionYears'].strip().split(',')]#[2019,2020,2021,2022]
    num_folds=int(dic['numFolds']) # 10
    num_subforests=int(dic['numSubForests']) #10
    annual_date_frac = float(dic['annualDateFrac']) #.5
    export_scale = int(dic['exportScale']) # 30
    export_crs = dic['exportCrs']# 'EPSG:6350'
    nodata_value = = int(dic['noDataValue']) # -9999
    
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
    
    region_name = dic['regionName']
    
    landcover_ic = dic['landcoverIc']
    landcoverKeepClasses = [int(x) for x in dic['landcoverKeepClasses'].strip().split(',')]
    landcoverPhenoVariesClasses = [int(x) for x in dic['landcoverPhenoVariesClasses'].strip().split(',')]
    gcs_mod_path  = dic['gcsModPath'] #     gcs_mod_path = 'gs://rf_mods/'+region_name+'/'
    
        
    ##########
    # PROCESSING
    ##########
    try:
        ee_reg = ee.FeatureCollection(region_name)
        print(collection.getMapId().get('mapid'))
        print('Layer' + region_name + 'found')
    except Exception as e:
      # Handle errors (e.g., asset doesn't exist)
      print(' ERROR, ee region cant loaded: ' + region_name)
      print(e)
      return False
    
    
    
    # region path
    region_geom = (ee.Feature(ee_reg.union(maxError=30).first()).geometry())
    
    # make the ancillary covariate image (only elev, slope, and aspect in this case)
    anc_static_cov_img = anc_cov.make_covar_topo_stack()
    
    # make the CCDC mosaic image
    if ccdc_result_path is not None:
        ccdc_img = ccdc_utils.mosaic_ccdc_img_tiles(ccdc_result_path)
    else:
        ccdc_img = None
    
    fromre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
    tore = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 10, 0, 0, 0, 0, 0, 11, 11, 11, 11, 12, 255]
    lcic = ee.Image( landcover_path ) # 
    
    landcoverKeepClasses =  list(set(tore))
    landcoverPhenoVariesClasses = landcoverKeepClasses
    
    landcover_ic = lcic.remap(from_= fromre, to=tore, bandName='b1',  defaultValue = None ).rename('lc')
        


    ## Original code here
    i = 0
    for metric in gedi_metrics:
        for year in projection_years:
    metric = 0
    year = 0
    region_metric_year_name = (f'rf_projection_{region_name}_{metric}_{year}')
    year_frac = ee.Number(year).add(ee.Number(annual_date_frac))
    

    # covariate image for this year
    cov_img = covs.make_covariate_stack_for_year(ancillary_static_cov_img=anc_static_cov_img,
                                                 year_frac=year_frac,
                                                 ccdc_img=ccdc_img,
                                                 ccdc_start_year=ccdc_start_year,
                                                 ccdc_seg_bands=ccdc_seg_bands,
                                                 ccdc_syn_bands=ccdc_syn_bands,
                                                 ccdc_tsin_bands=ccdc_tsin_bands,
                                                 ccdc_tsin_thresh=ccdc_tsin_thresh,
                                                 include_gse_covs=True)

    # TODO: predictors are renamed in R script to save space in text file; need a lookup table with covariate names and shortened IDs
    cov_numbers_0 = list(range(64))
    old_cov_names = []
    for n in cov_numbers_0:
        old_cov_names.append(f'gse_a{n:02d}')

    cov_numbers_1 = list(range(1,65))
    new_cov_names = []
    for n in cov_numbers_1:
        new_cov_names.append(f'p{n}')

    cov_img_rn = cov_img.select(old_cov_names, new_cov_names)

    # clip and mask the covariate image
    cov_img_cm = models.clip_mask_cov_img(cov_img=cov_img_rn,
                                          year=year,
                                          region_geom=region_geom,
                                          landcover_ic=ee.ImageCollection(landcover_ic),
                                          landcover_keep_classes=landcoverKeepClasses)

    # project the model
    projected_mean_sd = models.project_rf_text_mean_sd(gcs_mod_path=gcs_mod_path,
                                                       num_subforests=num_subforests,
                                                       metric=metric,
                                                       int16_metric_mult=gedi_metric_int16_mult[i],
                                                       cov_img=cov_img_cm)


    task = ee.batch.Export.image.toDrive(image=projected_mean_sd.unmask(nodata_value),
                                         description=region_metric_year_name+' to Google Drive',
                                         fileNamePrefix=region_metric_year_name,
                                         scale=export_scale,
                                         region=region_geom,
                                         crs=export_crs,
                                         maxPixels=1e13,
                                         skipEmptyTiles=True,
                                         fileFormat='GeoTIFF',
                                         formatOptions={'cloudOptimized': True,
                                                        'noData': nodata_value})

    task.start()
    
    i = i+1
    ############################################### ----
# End function




if __name__ == "__main__":
    main()

