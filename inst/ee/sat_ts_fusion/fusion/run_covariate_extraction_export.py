import ee
# initialize with Earth Engine project id
gee_project_name = 'ee-gedibio'
ee.Initialize(project=gee_project_name)

from sat_ts_fusion.ccdc import ccdc_utils
import ancillary_covariates as anc_cov
import covariates as covs



##########
# INPUTS
##########

region_name = 'reshape_wdm'

# point feature collection to use for covariate extraction
# needs to have a 'date' field with dates like '2020-06-25'
pt_fc_id = 'projects/ee-gedibio/assets/gediv002_l2l4a_20190417to20230316_'+region_name

# maximum number of points to extract covariates
ext_limit = 200000

# export asset ID path
pt_fc_ext_id = 'projects/ee-gedibio/assets/gediv002_l2l4a_20190417to20230316_covext_'+region_name

# path to CCDC results
# set to None to exclude CCDC covariates
ccdc_result_path= None
#'projects/ee-gedibio/assets/ccdc/results/workshop_sonoma/'

# CCDC start year
ccdc_start_year = 2000
# CCDC bands to extract segments coefficients for
ccdc_seg_bands = ['NBR2', 'NDMI', 'NDVI']
# CCDC bands to extract synthetic (i.e. fitted) values for
ccdc_syn_bands = ['NIR', 'SWIR1', 'SWIR2']
# CCDC bands to extract time since break values for
ccdc_tsin_bands = ['NBR2', 'NDVI']
# break threshold values associated with ccdc_tsin_bands
ccdc_tsin_thresh = [0.1, 0.25]



##########
# PROCESSING
##########

# load the point FeatureCollection, and limit collection size
pt_fc = ee.FeatureCollection(pt_fc_id).limit(maximum=ext_limit, prop='random_1')

# make the ancillary covariate image (only elev, slope, and aspect in this case)
anc_static_cov_img = anc_cov.make_covar_topo_stack()

# make the CCDC mosaic image
if ccdc_result_path is not None:
    ccdc_img = ccdc_utils.mosaic_ccdc_img_tiles(ccdc_result_path)
else:
    ccdc_img = None


# run the covariate extraction
fc_ext = covs.extract_covariates_at_points(fc=pt_fc,
                                          ancillary_static_cov_img=anc_static_cov_img,
                                          ccdc_img=ccdc_img,
                                          ccdc_start_year=ccdc_start_year,
                                          ccdc_seg_bands=ccdc_seg_bands,
                                          ccdc_syn_bands=ccdc_syn_bands,
                                          ccdc_tsin_bands=ccdc_tsin_bands,
                                          ccdc_tsin_thresh=ccdc_tsin_thresh,
                                          landcover_ic=ee.ImageCollection("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER"),
                                          landcover_keep_classes=[21,41,42,43,52,71,90,95],
                                          landcover_pheno_varies_classes=[41,43,90],
                                          include_gse_covs=True,
                                          ext_start_year=2019,
                                          ext_end_year=2022,
                                          ext_pheno_start_doy=152,
                                          ext_pheno_end_doy=243,
                                          ext_scale=30,
                                          temporal_match='annual',
                                          annual_date_frac=0.5)

# export the feature collection with covariates
task = ee.batch.Export.table.toAsset(collection=fc_ext,
                                     description=region_name+' FC covariate extract to GEE Asset',
                                     assetId=pt_fc_ext_id)
task.start()