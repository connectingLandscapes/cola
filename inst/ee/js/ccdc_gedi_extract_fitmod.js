/*
By: Patrick Burns [pb463@nau.edu], Northern Arizona University
    Some code provided by USFS (Zhiqiang Yang and GTAC) and Paulo Arevalo

About: explore various dimensions of CCDC result images

*/



// ----------------------
// ----- IMPORTS -----
// ----------------------
//var palettes = require('users/gena/packages:palettes');
var pal_viridis = ['#440154', '#433982', '#30678D', '#218F8B', '#36B677', '#8ED542', '#FDE725']
var ccdc_gg = require('users/pb463/NAU_GoetzGroup:01_utilities/modules/ccdc')
Map.setOptions('SATELLITE')



// ----------------------
// ----- INPUTS -----
// ----------------------
// region of interest
var region_geom = ee.Feature(ee.FeatureCollection('users/pb463/S2L/Sonoma_cty_v2_PBcleaned').first())
                    .geometry()

// path to CCDC assets (a directory of CCDC array images)
var ccdc_assets_path = 'projects/ee-gedibio/assets/ccdc/results/workshop_sonoma'

// names of bands that were included in the CCDC run
var bands = ['green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
var indices = ['NDVI', 'NDMI', 'NBR2']

// CCDC start and end years
var start_year = 1985
var end_year = 2024

// false color visualiztion params for synthetic Landsat images
var lsat_fcolor_viz = {bands: ['SWIR1', 'NIR', 'red'], min:0, max:0.45}

// GEDI asset
var gedi_fc = ee.FeatureCollection('projects/ee-gedibio/assets/gedi/vector_uploads/L2AB_002_2019-04-01_2021-09-30_DOY91_273_US-CA-Sonoma3_shifted')
                .filterBounds(region_geom)
                .filter(ee.Filter.gte('doy',150))
                .filter(ee.Filter.lt('doy',270))

var gedi_y_lim = 1000

var gedi_predext_asset = 'projects/ee-gedibio/assets/gedi_ccdc_predext_sonoma_tswork_15k'

var als_chm = ee.Image('users/pb463/als_ca-sonoma2013_chm_mos_3m')
                .rename('ALS_CHM')



// ----------------------
// ----- PROCESSING -----
// ----------------------

// load the ccdc result assets as a mosaiced image
var ccdc_assets = ee.data.listAssets(ccdc_assets_path).assets
                         .map(function(i){
                           return i.id
                         })
var ccdc_img = ee.ImageCollection(ccdc_assets)
                 .mosaic()

//print('CCDC mosaic image structure:', ccdc_img)
//Map.addLayer(ccdc_img, null, 'CCDC Output', 0)

var ccdc_syn = ccdc_gg.getSynthetic(ccdc_img, 2020.6, bands)
Map.addLayer(ccdc_syn, lsat_fcolor_viz, 'LS Synthetic, Summer 2020')
Map.addLayer(gedi_fc, {'color': 'red'}, 'GEDI Shots', 0)

var coefs_img = ccdc_gg.getSegmentParams(ccdc_img, 2020.6, indices, true, true)
var phamp_ind = ccdc_gg.getPhaseAmplitude(coefs_img, '.*SIN.*', '.*COS.*')
var tsince_brk_nbr0p1 = ccdc_gg.get_largestmag_tbreak_thresh_img(ccdc_img, start_year, 2020.6, 'NBR2', 0.1)
                              .rename('ccbrk_tsin_nbr0p1')

function rename_coef_bands(img){
  var new_names = img.bandNames().map(function(b){return ee.String('ccsegco_').cat(ee.String(b))})
  return img.rename(new_names)
}

function rename_phamp_bands(img){
  var new_names = img.bandNames().map(function(b){return ee.String('ccsegpa_').cat(ee.String(b))})
  return img.rename(new_names)
}

function rename_syn_bands(img){
  var new_names = img.bandNames().map(function(b){return ee.String('ccsynth_').cat(ee.String(b))})
  return img.rename(new_names)
}

if (gedi_predext_asset == null){
var gedi_fc_lim = ee.FeatureCollection(gedi_fc.filter(ee.Filter.eq('year', 2019)).limit(gedi_y_lim))
                    .merge(ee.FeatureCollection(gedi_fc.filter(ee.Filter.eq('year', 2020)).limit(gedi_y_lim)))
                    .merge(ee.FeatureCollection(gedi_fc.filter(ee.Filter.eq('year', 2021)).limit(gedi_y_lim)))

// prepare the predictor stack and extract for different years
var gedi_years_ls = ee.List([2019.5, 2020.5, 2021.5])
var gedi_predext_fc = ee.FeatureCollection(gedi_years_ls.map(function(y){
  var year_frac = ee.Number(y)
  // coefs and synthetics for the summer
  var coefs_all_img = rename_coef_bands(ccdc_gg.getSegmentParams(ccdc_img, year_frac, bands, true, true))
  var coefs_ind_img = ccdc_gg.getSegmentParams(ccdc_img, year_frac, indices, true, true)
  var phamp_ind = rename_phamp_bands(ccdc_gg.getPhaseAmplitude(coefs_ind_img, '.*SIN.*', '.*COS.*'))
  var syn = rename_syn_bands(ccdc_gg.getSynthetic(ccdc_img, year_frac, bands))

  // break info
  var n_breaks_before = ccdc_gg.get_n_breaks_img(ccdc_img, start_year, year_frac).rename('ccbreak_nbrk_before')
  var year_brk_nbr0p1 = ccdc_gg.get_largestmag_tbreak_thresh_img(ccdc_img, start_year, year_frac, 'NBR2', 0.1)
  var tsin_brk_nbr0p1 = ee.Image(year_frac).subtract(year_brk_nbr0p1.unmask(start_year))
                          .rename('ccbreak_tsin_nbr0p1')
  var ccdc_y_stack = ee.Image.cat([coefs_all_img, phamp_ind, syn, n_breaks_before, tsin_brk_nbr0p1])

  // extract at GEDI shots
  var gedi_fc_y = gedi_fc_lim.filter(ee.Filter.eq('year', year_frac.floor()))
  var gedi_fc_y_ext = ee.FeatureCollection(gedi_fc_y.map(function(f){
    var gedi_f_rr = ccdc_y_stack.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: ee.Feature(f).buffer(20,5).geometry(),
      scale: 30,
      tileScale: 8
    })
    return(f.set(gedi_f_rr))
  }))

  return gedi_fc_y_ext
})).flatten()

/*
Export.table.toAsset({
  collection: gedi_predext_fc,
  description: 'gedi_predext_sonoma_tswork_'+gedi_y_lim,
  assetId: 'projects/ee-gedibio/assets/gedi_ccdc_predext_sonoma_tswork_'+gedi_y_lim})
*/
} else {
  var gedi_predext_fc = ee.FeatureCollection(gedi_predext_asset)
}

// CCDC predictor stack for the prediction year
var year_frac_pred = 2013.5
// coefs and synthetics for the summer
var coefs_all_img = rename_coef_bands(ccdc_gg.getSegmentParams(ccdc_img, year_frac_pred, bands, true, true))
var coefs_ind_img = ccdc_gg.getSegmentParams(ccdc_img, year_frac_pred, indices, true, true)
var phamp_ind = rename_phamp_bands(ccdc_gg.getPhaseAmplitude(coefs_ind_img, '.*SIN.*', '.*COS.*'))
var syn = rename_syn_bands(ccdc_gg.getSynthetic(ccdc_img, year_frac_pred, bands))
// break info
var n_breaks_before = ccdc_gg.get_n_breaks_img(ccdc_img, start_year, year_frac_pred).rename('ccbreak_nbrk_before')
var year_brk_nbr0p1 = ccdc_gg.get_largestmag_tbreak_thresh_img(ccdc_img, start_year, year_frac_pred, 'NBR2', 0.1)
var tsin_brk_nbr0p1 = ee.Image(year_frac_pred).subtract(year_brk_nbr0p1.unmask(start_year))
                        .rename('ccbreak_tsin_nbr0p1')
var pred_stack = ee.Image.cat([coefs_all_img, phamp_ind, syn, n_breaks_before, tsin_brk_nbr0p1])
var pred_names = pred_stack.bandNames()

/* all pred names
[
  "ccsegco_green_INTP",  "ccsegco_green_SLP",  "ccsegco_green_COS",  "ccsegco_green_SIN",  "ccsegco_green_COS2",  "ccsegco_green_SIN2",  "ccsegco_green_COS3",  "ccsegco_green_SIN3",
  "ccsegco_red_INTP",  "ccsegco_red_SLP",  "ccsegco_red_COS",  "ccsegco_red_SIN",  "ccsegco_red_COS2",  "ccsegco_red_SIN2",  "ccsegco_red_COS3",  "ccsegco_red_SIN3",
  "ccsegco_NIR_INTP",  "ccsegco_NIR_SLP",  "ccsegco_NIR_COS",  "ccsegco_NIR_SIN",  "ccsegco_NIR_COS2",  "ccsegco_NIR_SIN2",  "ccsegco_NIR_COS3",  "ccsegco_NIR_SIN3",
  "ccsegco_SWIR1_INTP",  "ccsegco_SWIR1_SLP",  "ccsegco_SWIR1_COS",  "ccsegco_SWIR1_SIN",  "ccsegco_SWIR1_COS2",  "ccsegco_SWIR1_SIN2",  "ccsegco_SWIR1_COS3",  "ccsegco_SWIR1_SIN3",
  "ccsegco_SWIR2_INTP",  "ccsegco_SWIR2_SLP",  "ccsegco_SWIR2_COS",  "ccsegco_SWIR2_SIN",  "ccsegco_SWIR2_COS2",  "ccsegco_SWIR2_SIN2",  "ccsegco_SWIR2_COS3",  "ccsegco_SWIR2_SIN3",
  "ccsegco_NDVI_INTP",  "ccsegco_NDVI_SLP",  "ccsegco_NDVI_COS",  "ccsegco_NDVI_SIN",  "ccsegco_NDVI_COS2",  "ccsegco_NDVI_SIN2",  "ccsegco_NDVI_COS3",  "ccsegco_NDVI_SIN3",
  "ccsegco_NDMI_INTP",  "ccsegco_NDMI_SLP",  "ccsegco_NDMI_COS",  "ccsegco_NDMI_SIN",  "ccsegco_NDMI_COS2",  "ccsegco_NDMI_SIN2",  "ccsegco_NDMI_COS3",  "ccsegco_NDMI_SIN3",
  "ccsegco_NBR2_INTP",  "ccsegco_NBR2_SLP",  "ccsegco_NBR2_COS",  "ccsegco_NBR2_SIN",  "ccsegco_NBR2_COS2",  "ccsegco_NBR2_SIN2",  "ccsegco_NBR2_COS3",  "ccsegco_NBR2_SIN3",
  "ccsegpa_NDVI_PHASE",  "ccsegpa_NDVI_PHASE2",  "ccsegpa_NDVI_PHASE3",  "ccsegpa_NDMI_PHASE",  "ccsegpa_NDMI_PHASE2",  "ccsegpa_NDMI_PHASE3",  "ccsegpa_NBR2_PHASE",  "ccsegpa_NBR2_PHASE2",  "ccsegpa_NBR2_PHASE3",
  "ccsegpa_NDVI_AMPLITUDE",  "ccsegpa_NDVI_AMPLITUDE2",  "ccsegpa_NDVI_AMPLITUDE3",  "ccsegpa_NDMI_AMPLITUDE",  "ccsegpa_NDMI_AMPLITUDE2",  "ccsegpa_NDMI_AMPLITUDE3",  "ccsegpa_NBR2_AMPLITUDE",  "ccsegpa_NBR2_AMPLITUDE2",  "ccsegpa_NBR2_AMPLITUDE3",
  "ccsynth_green",  "ccsynth_red",  "ccsynth_NIR",  "ccsynth_SWIR1",  "ccsynth_SWIR2",  "ccsynth_NDVI",  "ccsynth_NDMI",  "ccsynth_NBR2",
  "ccbreak_nbrk_before",  "ccbreak_tsin_nbr0p1"
]
*/

var sel_pred_names = ["ccsegco_green_INTP",  "ccsegco_green_SLP", "ccsegco_red_INTP",  "ccsegco_red_SLP", "ccsegco_SWIR1_INTP",  "ccsegco_SWIR1_SLP", "ccsegco_SWIR2_INTP",  "ccsegco_SWIR2_SLP",
                     "ccsegco_NDVI_INTP",  "ccsegco_NDVI_SLP", "ccsegco_NDMI_INTP",  "ccsegco_NDMI_SLP", "ccsegco_NBR2_INTP",  "ccsegco_NBR2_SLP",
                     "ccsegpa_NDVI_PHASE", "ccsegpa_NDVI_AMPLITUDE", "ccsegpa_NDMI_PHASE", "ccsegpa_NDMI_AMPLITUDE", "ccsegpa_NBR2_PHASE", "ccsegpa_NBR2_AMPLITUDE",
                     "ccsynth_green",  "ccsynth_red",  "ccsynth_NIR",  "ccsynth_SWIR1",  "ccsynth_SWIR2",  "ccsynth_NDVI",  "ccsynth_NDMI",  "ccsynth_NBR2",
                     "ccbreak_nbrk_before",  "ccbreak_tsin_nbr0p1"]

// RF model
var gedi_predext_fc = gedi_predext_fc.randomColumn('rid',123)

// Train RF
var rf_train = ee.Classifier.smileRandomForest({numberOfTrees: 200, seed: 123})
                 .setOutputMode('REGRESSION')
                 .train({features: gedi_predext_fc.filter(ee.Filter.notNull(sel_pred_names))
                                                  .limit(14000, 'rid', true),
                         classProperty: 'rh_98_a0',
                         inputProperties: sel_pred_names})
//print(rf_train.explain())


var rf_pred_fc = gedi_predext_fc.filter(ee.Filter.notNull(sel_pred_names))
                                .limit(1000, 'rid', false)
                                .classify(rf_train)
                                .filter(ee.Filter.notNull(['classification']))


// Make a chart comparing model predictions with observed data
    var pred_fc = rf_pred_fc.aggregate_array('classification')
    var obs = rf_pred_fc.aggregate_array('rh_98_a0')
    var scatter_chart = ui.Chart.array.values(pred_fc, 0, obs)
      .setSeriesNames(['RH98'])
      .setOptions({
        title: "RF Prediction of GEDI RH98",
        hAxis: {title: "Observed",
        viewWindow: {min: 0, max: 50}
        },
        vAxis: {title: "Predicted",
        viewWindow: {min: 0, max: 50}
        },
        pointSize: 6,
        dataOpacity: 0.4,
        trendlines: {
          0: {
                type: 'linear',
                showR2: true,
                color: 'black',
                visibleInLegend: true
          }
        }
    });
    print(scatter_chart)


var rf_pred_img = pred_stack.classify(rf_train)
                            .clip(region_geom)
                            .rename('GEDI_RH98')
Map.addLayer(rf_pred_img, {min:0, max: 50, palette:pal_viridis}, 'RF Predicted RH98')


Export.image.toAsset({
  image: rf_pred_img,
  description: 'gedi_ccdc_rfpred_sonoma_workshop_test',
  assetId: 'projects/ee-gedibio/assets/gedi_ccdc_rfpred_sonoma_workshop_test',
  region: region_geom,
  scale: 30,
  })

Map.addLayer(als_chm,{min:0,max:50, palette:pal_viridis},'ALS CHM', 0)

    var RF_varimp = ee.Feature(null, ee.Dictionary(rf_train.explain()).get('importance'));
    var RF_imp_chart =
      ui.Chart.feature.byProperty(RF_varimp)
        .setChartType('ColumnChart')
        .setOptions({
          title: 'Random Forest Variable Importance',
          legend: {position: 'none'},
          hAxis: {title: 'Predictors'},
          vAxis: {title: 'Importance'}
      });
    print(RF_imp_chart)

print(
  ui.Chart.image.histogram({
    image: ee.Image.cat([rf_pred_img, als_chm]),
    region: region_geom,
    scale: 90,
    maxBuckets: 64,
    minBucketWidth: 2,
    //maxRaw: 1e5,
    maxPixels: 1e10})
)