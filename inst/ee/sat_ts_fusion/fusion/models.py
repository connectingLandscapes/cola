import ee
# initialize with Earth Engine project id
# gee_project_name = 'ee-gedibio'
# ee.Initialize(project=gee_project_name)


def project_rf_text_mean_sd(gcs_mod_path: str = None,
                            num_subforests: int = 10,
                            metric: str = None,
                            int16_metric_mult: int = None,
                            cov_img: ee.Image = None):
    """Project the mean and SD of all RF model subforests that are stored in Google Cloud Storage

    :param gcs_mod_path: the Google Cloud Storage path where model text files are stored
    :param num_subforests: number of subforest model text files
    :param metric: the metric to be modeled
    :param int16_metric_mult: an integer to multiply the projected values by to convert to int16 (to reduce space)
    :param cov_img: covariate image to use for projection
    :return: ee.Image of the projected RF model mean and SD across n folds
    """

    gcs_path_template = (f'{gcs_mod_path}{metric}_rf1_mod_subforest')
    gcs_paths_fold = [gcs_path_template + str(i) + '.txt' for i in list(range(1,num_subforests+1))]

    def subforest_path_to_str(subforest_path):
        """Convert subforest path to string

        :param subforest_path:
        :return:
        """
        return ee.Blob(subforest_path).string()

    model_strings_list = list(map(subforest_path_to_str, gcs_paths_fold))
    classifier_regr_raw = ee.Classifier.decisionTreeEnsemble(model_strings_list).setOutputMode('RAW_REGRESSION')

    projected_mn_img = cov_img.classify(classifier_regr_raw).arrayReduce(ee.Reducer.mean(), ee.List([0])).arrayGet(0)
    projected_sd_img = cov_img.classify(classifier_regr_raw).arrayReduce(ee.Reducer.stdDev(), ee.List([0])).arrayGet(0)

    projected_mean_sd = (ee.Image.cat([projected_mn_img, projected_sd_img])
                           .multiply(int16_metric_mult)
                           .int16()
                           .rename([f'{metric}_mn', f'{metric}_sd']))

    return projected_mean_sd


def clip_mask_cov_img(cov_img: ee.Image = None,
                      year: int = None,
                      region_geom: ee.Geometry = None,
                      landcover_ic: ee.ImageCollection = None,
                      landcover_keep_classes: list = None):
    """Clip and mask a model projection using a region geometry and landcover dataset

    :param cov_img: covariate image
    :param year: projection year for aligning with landcover
    :param region_geom: ee.Geometry for the region of interest
    :param landcover_ic: ee.ImageCollection landcover dataset
    :param landcover_keep_classes: landcover classes to keep for the mask
    :return: ee.Image model projection that has been clipped and masked
    """
    img_c = cov_img.clip(region_geom)

    if landcover_ic is not None:
        lc_y_mask = (ee.Image(landcover_ic.filter(ee.Filter.eq('year', year)).first())
                       .remap(landcover_keep_classes, [1] * len(landcover_keep_classes)))
        img_cm = img_c.updateMask(lc_y_mask)
    else:
        img_cm = img_c

    return img_cm



