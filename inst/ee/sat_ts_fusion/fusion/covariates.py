import ee
# initialize with Earth Engine project id
# gee_project_name = 'ee-gedibio'
# ee.Initialize(project=gee_project_name)

from sat_ts_fusion.ccdc import ccdc_utils
from sat_ts_fusion.fusion import ancillary_covariates
#as anc_cov



def make_covariate_stack_for_year(ancillary_static_cov_img: ee.Image = None,
                                  year_frac: float = None,
                                  ccdc_img = None,
                                  ccdc_start_year: int = None,
                                  ccdc_seg_bands: list = None,
                                  ccdc_syn_bands: list = None,
                                  ccdc_tsin_bands: list = None,
                                  ccdc_tsin_thresh: list = None,
                                  include_gse_covs: bool = False):
    """Make a multiband image of covariates
    :param ancillary_static_cov_img: an ee.Image with ancillary covariates that are temporally static
    :param year_frac: the fractional year to extract covariate values at
    :param ccdc_img: a CCDC result image
    :param ccdc_start_year: the year CCDC started
    :param ccdc_seg_bands: CCDC segment bands list
    :param ccdc_syn_bands: CCDC synthetic bands list
    :param ccdc_tsin_bands: CCDC time since break bands list
    :param ccdc_tsin_thresh: CCDC time since break thresholods list
    :param include_gse_covs: whether to include Google Satellite Embeddings covariates
    :return: multiband ee.Image of covariates
    """
    cov_stack = ancillary_static_cov_img
    if ccdc_img is not None:
        # CCDC segment coefficients
        ccdc_seg_coefs = ccdc_utils.get_segment_coefs(ccdc_img=ccdc_img, year_frac=year_frac,bands=ccdc_seg_bands, normalize=True, use_next=False)
        # CCDC segment phase and amplitude
        ccdc_seg_phamp = ccdc_utils.get_segment_phase_amp(coefs_img=ccdc_seg_coefs)
        # CCDC synthetic values
        ccdc_syn_vals = ccdc_utils.get_synthetic(ccdc_img=ccdc_img, year_frac=year_frac,
                                                 bands=ccdc_syn_bands)
        # CCDC segment breaks
        # number of breaks before the site observation
        ccdc_n_breaks_before = ccdc_utils.get_n_breaks_img(ccdc_img=ccdc_img,
                                                           start_year=ccdc_start_year,
                                                           year_frac=year_frac).rename('ccbreak_nbrk_before')
        ccdc_time_since = ccdc_utils.get_time_since_largest_mag_break_bands(ccdc_img=ccdc_img, start_year=ccdc_start_year,year_frac=year_frac, bands_list=ccdc_tsin_bands,thresh_list=ccdc_tsin_thresh)
        ccdc_cov_img = ee.Image.cat([ccdc_seg_coefs, ccdc_seg_phamp, ccdc_syn_vals,
                                     ccdc_n_breaks_before,
                                     ccdc_time_since])
        # combine all images into one, including ancillary covariates and embeddings
        cov_stack = cov_stack.addBands(ccdc_cov_img)
        # close if
        
    
    if include_gse_covs:
        sat_embed = ancillary_covariates.get_google_sat_embed(ee.Number(year_frac).floor())
        cov_stack = cov_stack.addBands(sat_embed)
    return cov_stack


def extract_covariates_at_points(fc: ee.FeatureCollection = None,
                                 ancillary_static_cov_img: ee.Image = None,
                                 ccdc_img: ee.Image = None,
                                 ccdc_start_year: int = None,
                                 ccdc_seg_bands: list = None,
                                 ccdc_syn_bands: list = None,
                                 ccdc_tsin_bands: list = None,
                                 ccdc_tsin_thresh: list = None,
                                 landcover_ic: ee.ImageCollection = None,
                                 landcover_keep_classes: list = None,
                                 landcover_pheno_varies_classes: list = None,
                                 include_gse_covs: bool = False,
                                 ext_start_year:int = None,
                                 ext_end_year:int = None,
                                 ext_start_doy:int = 1,
                                 ext_end_doy:int = 366,
                                 ext_pheno_start_doy: int = None,
                                 ext_pheno_end_doy: int = None,
                                 ext_scale: float = 30.0,
                                 temporal_match:str = 'annual',
                                 annual_date_frac:float = 0.5):
    """Extract CCDC and ancillary covariates at point geometries
    :param fc: ee.FeatureCollection with point geometries and date field
    :param ccdc_img: CCDC array image
    :param ccdc_start_year: starting year for CCDC run
    :param ccdc_seg_bands: the bands/indices to use for getting CCDC segment coefficients
    :param ccdc_syn_bands: the bands/indices to use for getting CCDC synthethic values
    :param ccdc_tsin_bands: the bands/indices to use for getting time since largest absolute change between CCDC segments
    :param ccdc_tsin_thresh: the thresholds to use for each band/index in ccdc_tsin_bands
    :param landcover_ic
    :param landcover_keep_classes
    :param include_gse_covs
    :param ext_start_year: first year to use for extracting covariates
    :param ext_end_year: last year to use for extracting covariates
    :param ext_start_doy: first day of year to use for extracting covariates
    :param ext_end_doy: last day of year to use for extracting covariates
    :param temporal_match: the scenario to use for temporally matching dates with CCDC. Options are 'annual' or 'day' (the slower option)
    :param annual_date_frac: if temporal_match is 'annual', this is the fraction of the year to add to each year
    :return: ee.FeatureCollection with covariate fields added
    """
    
    # add relevant date fields
    def add_date_fields(feat: ee.Feature = None):
        """Add additional date fields to a feature
        :param feat: ee.Feature
        :return: ee.Feature with additional date properties
        """
        date_fm = ee.Date.parse(format='yyyy-MM-dd', date=feat.get('date'))
        date_millis = date_fm.millis()
        year_frac = date_fm.get('year').add(date_fm.getFraction('year'))
        return feat.set('system:time_start', date_millis, 'year_frac', year_frac)
    # filter the point feature collection temporally
    fc_f = (fc.map(add_date_fields)
              .filter(ee.Filter.calendarRange(ext_start_year, ext_end_year, 'year'))
              .filter(ee.Filter.calendarRange(ext_start_doy, ext_end_doy, 'day_of_year')))
    # extract covariates at points
    # there are different extraction scenarios based on temporal matching strategy
    if temporal_match=='annual':
        years_list = ee.List.sequence(ext_start_year, ext_end_year)
        def extract_covariates_annually(year: int = None):
            """Extract values of covariates at points. The points will span the filtered date range specified above,
            but the CCDC covariates will correspond to one fractional year value.
            :param year: year
            :return: ee.FeatureCollection with covariate fields added
            """
            # get features for a year
            fc_y = fc_f.filter(ee.Filter.calendarRange(year, year, 'year'))
            # specify the fractional year value to use for getting CCDC coefficients and synthetics
            ext_year_frac = ee.Number(year).add(ee.Number(annual_date_frac))
            # grab landcover for the year
            landcover_img_y = ee.Image(landcover_ic.filter(ee.Filter.calendarRange(year, year, 'year')).first()).rename('lc')
            def extract_lc_feat(feat: ee.Feature = None):
                """Extract landcover at a point
                :param feat: ee.Feature
                :return: ee.Feature with landcover field added
                """
                lc_rr = landcover_img_y.reduceRegion(reducer=ee.Reducer.mode(),
                                                     geometry = ee.Feature(feat).geometry(),
                                                     scale=ext_scale,
                                                     tileScale=8)
                return feat.set(lc_rr)
            fc_y_lc = (ee.FeatureCollection(fc_y.map(extract_lc_feat))
                         .filter(ee.Filter.inList('lc',ee.List(landcover_keep_classes))))
            fc_y_pheno_stable = ee.FeatureCollection(fc_y_lc.filter(ee.Filter.inList('lc', landcover_pheno_varies_classes).Not()))
            fc_y_pheno_varies = ee.FeatureCollection(fc_y_lc.filter(ee.Filter.And(ee.Filter.inList('lc', landcover_pheno_varies_classes),ee.Filter.gte('doy', ext_pheno_start_doy),ee.Filter.lte('doy', ext_pheno_end_doy))))
            fc_y_lc_ph = fc_y_pheno_stable.merge(fc_y_pheno_varies)
            # make the covariate stack for the year
            cov_stack = make_covariate_stack_for_year(year_frac=ext_year_frac,
                                                      ancillary_static_cov_img=ancillary_static_cov_img,
                                                      ccdc_img=ccdc_img,
                                                      ccdc_start_year=ccdc_start_year,
                                                      ccdc_seg_bands=ccdc_seg_bands,
                                                      ccdc_syn_bands=ccdc_syn_bands,
                                                      ccdc_tsin_bands=ccdc_tsin_bands,
                                                      ccdc_tsin_thresh=ccdc_tsin_thresh,
                                                      include_gse_covs=include_gse_covs)
            def extract_covariates_feat(feat: ee.Feature = None):
                """Extract values of covariate at a point
                :param feat: ee.Feature
                :return: ee.Feature with covariate fields added
                """
                cov_rr = cov_stack.reduceRegion(reducer=ee.Reducer.mean(),
                                                geometry=ee.Feature(feat).geometry(),
                                                scale=ext_scale,
                                                tileScale=8)
                return feat.set(cov_rr)
            return fc_y_lc_ph.map(extract_covariates_feat)
        fc_cov_ext = ee.FeatureCollection(years_list.map(extract_covariates_annually)).flatten()
    elif temporal_match=='day':
        fyears_list = fc_f.aggregate_array('year_frac').distinct().sort()
        def extract_covariates_day(year_frac: float = None):
            """Extract values of covariates at points. The CCDC covariates will correspond to the date associated with the feature.
            :param year_frac: fractional year
            :return: ee.FeatureCollection with covariate fields added
            """
            # get features for associated with a day
            fc_d = fc_f.filter(ee.Filter.eq('year_frac', year_frac))
            # grab landcover for the year
            year = ee.Number(year_frac).floor()
            landcover_img_y = ee.Image(landcover_ic.filter(ee.Filter.calendarRange(year, year, 'year')).first()).rename('lc')
            def extract_lc_feat(feat: ee.Feature = None):
                """Extract landcover at a point
                :param feat: ee.Feature
                :return: ee.Feature with landcover field added
                """
                lc_rr = landcover_img_y.reduceRegion(reducer=ee.Reducer.mode(),
                                                     geometry=ee.Feature(feat).geometry(),
                                                     scale=ext_scale,
                                                     tileScale=8)
                return feat.set(lc_rr)
            
            fc_d_lc = (ee.FeatureCollection(fc_d.map(extract_lc_feat))
                         .filter(ee.Filter.inList('lc', ee.List(landcover_keep_classes))))
            
            fc_d_pheno_stable = ee.FeatureCollection(fc_d_lc.filter(ee.Filter.inList('lc', landcover_pheno_varies_classes).Not()))
            fc_d_pheno_varies = ee.FeatureCollection(fc_d_lc.filter(ee.Filter.And(ee.Filter.inList('lc', landcover_pheno_varies_classes), ee.Filter.gte('doy', ext_pheno_start_doy), ee.Filter.lte('doy', ext_pheno_end_doy))))
            fc_d_lc_ph = fc_d_pheno_stable.merge(fc_d_pheno_varies)
            
            # make the covariate stack
            cov_stack = make_covariate_stack_for_year(year_frac=year_frac,
                                                      ancillary_static_cov_img=ancillary_static_cov_img,
                                                      ccdc_img=ccdc_img,
                                                      ccdc_start_year=ccdc_start_year,
                                                      ccdc_seg_bands=ccdc_seg_bands,
                                                      ccdc_syn_bands=ccdc_syn_bands,
                                                      ccdc_tsin_bands=ccdc_tsin_bands,
                                                      ccdc_tsin_thresh=ccdc_tsin_thresh,
                                                      include_gse_covs=include_gse_covs)
            
            def extract_covariates_feat(feat: ee.Feature = None):
                """Extract values of covariates at a point
                :param feat: ee.Feature
                :return: ee.Feature with covariate fields added
                """
                cov_rr = cov_stack.reduceRegion(reducer=ee.Reducer.mean(),
                                                geometry=ee.Feature(feat).geometry(),
                                                scale=ext_scale, tileScale=8)
                return feat.set(cov_rr)
            return fc_d_lc_ph.map(extract_covariates_feat) 
        fc_cov_ext = ee.FeatureCollection(fyears_list.map(extract_covariates_day)).flatten()
    return fc_cov_ext



def extract_covariates_at_points2(fc: ee.FeatureCollection = None,
                                 ancillary_static_cov_img: ee.Image = None,
                                 ccdc_img: ee.Image = None,
                                 ccdc_start_year: int = None,
                                 ccdc_seg_bands: list = None,
                                 ccdc_syn_bands: list = None,
                                 ccdc_tsin_bands: list = None,
                                 ccdc_tsin_thresh: list = None,
                                 landcover_ic: ee.ImageCollection = None,
                                 landcover_keep_classes: list = None,
                                 landcover_pheno_varies_classes: list = None,
                                 include_gse_covs: bool = False,
                                 ext_start_year:int = None,
                                 ext_end_year:int = None,
                                 ext_start_doy:int = 1,
                                 ext_end_doy:int = 366,
                                 ext_pheno_start_doy: int = None,
                                 ext_pheno_end_doy: int = None,
                                 ext_scale: float = 30.0,
                                 annual_date_frac:float = 0.5):
    """Extract CCDC and ancillary covariates at point geometries
    :param fc: ee.FeatureCollection with point geometries and date field
    :param ccdc_img: CCDC array image
    :param ccdc_start_year: starting year for CCDC run
    :param ccdc_seg_bands: the bands/indices to use for getting CCDC segment coefficients
    :param ccdc_syn_bands: the bands/indices to use for getting CCDC synthethic values
    :param ccdc_tsin_bands: the bands/indices to use for getting time since largest absolute change between CCDC segments
    :param ccdc_tsin_thresh: the thresholds to use for each band/index in ccdc_tsin_bands
    :param landcover_ic
    :param landcover_keep_classes
    :param include_gse_covs
    :param ext_start_year: first year to use for extracting covariates
    :param ext_end_year: last year to use for extracting covariates
    :param ext_start_doy: first day of year to use for extracting covariates
    :param ext_end_doy: last day of year to use for extracting covariates
    :param temporal_match: the scenario to use for temporally matching dates with CCDC. Options are 'annual' or 'day' (the slower option)
    :param annual_date_frac: if temporal_match is 'annual', this is the fraction of the year to add to each year
    :return: ee.FeatureCollection with covariate fields added
    """
    
    # add relevant date fields
    def add_date_fields(feat: ee.Feature = None):
        """Add additional date fields to a feature
        :param feat: ee.Feature
        :return: ee.Feature with additional date properties
        """
        # feat = fc.first(); feat.getInfo()
        
        date_fm = ee.Date.parse(format='yyyy-MM-dd', date=feat.get('date')) # date_fm.getInfo() 
        year_frac = date_fm.get('year').add(date_fm.getFraction('year'))
        date_millis = date_fm.millis()
        finalfeat = feat.set('system:time_start', date_millis, 'year_frac', year_frac) # finalfeat.getInfo()
        return finalfeat 
    
    # filter the point feature collection temporally
    fc_f = (fc.map(add_date_fields)
          .filter(ee.Filter.calendarRange(ext_start_year, ext_end_year, 'year'))
          .filter(ee.Filter.calendarRange(ext_start_doy, ext_end_doy, 'day_of_year')))
    
    #print( 'fc_f.size().getInfo()' )
    #print( fc_f.size().getInfo() )
    # fc_f.size().getInfo() 
    # fc_f_gi = fc_f.getInfo() 
    
    # extract covariates at points
    # there are different extraction scenarios based on temporal matching strategy
    # years_list = ee.List.sequence(ext_start_year, ext_end_year)
    #years_list = 0
    years_list = ee.List.sequence(ext_end_year, ext_end_year) # TEMPORAL DEBUG
    #print('years_list.getInfo()')
    #print(years_list.getInfo())
    
    # years_list.getInfo() 
    def extract_covariates_annually(year: int = None):
        """Extract values of covariates at points. The points will span the filtered date range specified above,
        but the CCDC covariates will correspond to one fractional year value.
        :param year: year
        :return: ee.FeatureCollection with covariate fields added
        """
        # get features for a year
        # year = 2019
        # fc_y = fc_f.filter(ee.Filter.calendarRange(year, year, 'date'))
        fc_y = fc_f ## TEMP SOLUTION
        # specify the fractional year value to use for getting CCDC coefficients and synthetics
        ext_year_frac = ee.Number(year).add(ee.Number(annual_date_frac))
        # grab landcover for the year
        # landcover_img_y = ee.Image(landcover_ic.filter(ee.Filter.calendarRange(year, year, 'year')).first()).rename('lc')
        landcover_img_y = landcover_ic
        def extract_lc_feat(feat: ee.Feature = None):
            """Extract landcover at a point
            :param feat: ee.Feature
            :return: ee.Feature with landcover field added
            """
            # feat = fc.first() # feat.getInfo() 
            lc_rr = landcover_img_y.reduceRegion(reducer=ee.Reducer.mode(),
                                                 geometry = ee.Feature(feat).geometry(),
                                                 scale=ext_scale,
                                                 tileScale=8)
            # lc_rr.getInfo()
            return feat.set(lc_rr)
        #
        fromre = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
        tore = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 10, 0, 0, 0, 0, 0, 11, 11, 11, 11, 12, 255]
        
        #fc2 = ee.FeatureCollection(fc).map(extract_lc_feat)
        #fc2.first().getInfo()
        #fc3 = fc2.remap(fromre, tore, 'lc');
        #fc3.first().getInfo()
        fc_y_lc_ph = ee.FeatureCollection(fc).map(extract_lc_feat)
        fc_y_lc_ph = fc_y_lc_ph.remap(fromre, tore, 'lc') 
        fc_y_lc_ph = fc_y_lc_ph.filter(ee.Filter.inList('lc',ee.List(landcover_keep_classes)))
        # fc_y_lc.size().getInfo() 
        # fc_y_pheno_stable = ee.FeatureCollection(fc_y_lc.filter(ee.Filter.inList('lc', landcover_pheno_varies_classes).Not()))
        # fc_y_pheno_varies = ee.FeatureCollection(fc_y_lc.filter(ee.Filter.And(ee.Filter.inList('lc', landcover_pheno_varies_classes),ee.Filter.gte('doy', ext_pheno_start_doy),ee.Filter.lte('doy', ext_pheno_end_doy))))
        # fc_y_lc_ph = fc_y_pheno_stable.merge(fc_y_pheno_varies)
        # fc_y_lc_ph.size().getInfo()
        # print( 'fc_y_lc_ph.size().getInfo()')
        # print( fc_y_lc_ph.size().getInfo() )
        #print( 'year:' + str(year))

        # fc_y_lc_gi = fc_y_lc_ph.getInfo() 
        # make the covariate stack for the year
        cov_stack = make_covariate_stack_for_year(year_frac=ext_year_frac,
                                                  ancillary_static_cov_img=ancillary_static_cov_img,
                                                  include_gse_covs=include_gse_covs)
# =============================================================================
#                                                   ccdc_img=ccdc_img,
#                                                   ccdc_start_year=ccdc_start_year,
#                                                   ccdc_seg_bands=ccdc_seg_bands,
#                                                   ccdc_syn_bands=ccdc_syn_bands,
#                                                   ccdc_tsin_bands=ccdc_tsin_bands,
#                                                   ccdc_tsin_thresh=ccdc_tsin_thresh,
# =============================================================================
        
        def extract_covariates_feat(feat: ee.Feature = None):
            """Extract values of covariate at a point
            :param feat: ee.Feature
            :return: ee.Feature with covariate fields added
            """
            cov_rr = cov_stack.reduceRegion(reducer=ee.Reducer.mean(),
                                            geometry=ee.Feature(feat).geometry(),
                                            scale=ext_scale,
                                            tileScale=8)
            return feat.set(cov_rr)
        # extract_covariates_annually
        return fc_y_lc_ph.map(extract_covariates_feat)
        
        
    fc_cov_ext = ee.FeatureCollection(years_list.map(extract_covariates_annually)).flatten()
    # fc_cov_ext.size().getInfo() 
    # fc_cov_ext_gi =  fc_cov_ext.getInfo() 
    return fc_cov_ext
