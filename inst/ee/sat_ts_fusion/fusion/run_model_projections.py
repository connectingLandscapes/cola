import ee
# initialize with Earth Engine project id
gee_project_name = 'ee-gedibio'
ee.Initialize(project=gee_project_name)

from sat_ts_fusion.ccdc import ccdc_utils
import ancillary_covariates as anc_cov
import covariates as covs
import sat_ts_fusion.fusion.models as models

# initialize with Earth Engine project id
gee_project_name = 'ee-gedibio'
ee.Initialize(project=gee_project_name)


##########
# INPUTS
##########

region_name = 'reshape_mwcf'

gedi_metrics = ['agbd_a0', 'rh_98_a0']
gedi_metric_int16_mult = [10, 100]
projection_years = [2019,2020,2021,2022]
num_folds=10
num_subforests=10
annual_date_frac = 0.5
export_scale = 30
export_crs = 'EPSG:6350'
nodata_value = -9999

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

# region path
region_geom = (ee.Feature(ee.FeatureCollection('projects/'+gee_project_name+'/assets/region_'+region_name)
                            .union(maxError=30).first())
                 .geometry())

# make the ancillary covariate image (only elev, slope, and aspect in this case)
anc_static_cov_img = anc_cov.make_covar_topo_stack()

# make the CCDC mosaic image
if ccdc_result_path is not None:
    ccdc_img = ccdc_utils.mosaic_ccdc_img_tiles(ccdc_result_path)
else:
    ccdc_img = None

gcs_mod_path = 'gs://rf_mods/'+region_name+'/'

i = 0
for metric in gedi_metrics:
    for year in projection_years:
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
                                              landcover_ic=ee.ImageCollection("projects/sat-io/open-datasets/USGS/ANNUAL_NLCD/LANDCOVER"),
                                              landcover_keep_classes=[21, 41, 42, 43, 52, 71, 90, 95])

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