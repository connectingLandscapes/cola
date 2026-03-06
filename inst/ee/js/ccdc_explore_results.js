/*
By: Patrick Burns [pb463@nau.edu], Northern Arizona University
    Kathleen Orndahl
    Some code provided by USFS (Zhiqiang Yang and GTAC) and Paulo Arevalo

About: explore various dimensions of CCDC result images

*/



// ----------------------
// ----- IMPORTS -----
// ----------------------
var palettes = require('users/gena/packages:palettes');
Map.setOptions('SATELLITE')



// ----------------------
// ----- INPUTS -----
// ----------------------
// path to CCDC assets (a directory of CCDC array images)
var ccdc_assets_path = 'projects/ee-gedibio/assets/ccdc/results/workshop_sonoma'

// names of bands that were included in the CCDC run
var bands = ['green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
// CCDC start and end years
var start_year = 1985
var end_year = 2024

// a fractional year for visualizing CCDC segments and synthetic imagery
var year_frac = 2000.5

// false color visualiztion params for synthetic Landsat images
var lsat_fcolor_viz = {bands: ['SWIR1', 'NIR', 'red'], min:0, max:0.45}



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

print('CCDC mosaic image structure:', ccdc_img)
Map.addLayer(ccdc_img, null, 'CCDC Output', 0)


////
// code from Zhiqiang Yang (USFS)
////

// normalize the intercept coefficient
var normalizeIntercept = function(ccdc_img, bands) {
  var result = []
  bands.forEach(function(band) {
    result.push(normalizeInterceptByBand(ccdc_img, band))
  })

  return ccdc_img.addBands(ee.Image(result), null, true)
}

var normalizeInterceptByBand = function(ccdc_img, band) {
  var mid = ccdc_img.select('tStart').add(ccdc_img.select('tEnd')).divide(2).toArray(1)

  var coefs = ccdc_img.select(band + '_coefs')
  var adjusted = coefs.arraySlice(1, 0, 1).add(coefs.arraySlice(1, 1, 2).multiply(mid))

  return adjusted.arrayCat(coefs.arraySlice(1, 1), 1).rename(band + '_coefs')
}

// get a multiband coefficient image associated with a time/segment
// TODO: possible to get closest segment in addition to previous or next?
var getSegmentParamsByBandImpl = function(ccdc_img, year_frac, band, useNext) {
  var tStart = ccdc_img.select('tStart')
  var tEnd = ccdc_img.select('tEnd')
  var coefs = ccdc_img.select(band + '_coefs')

  //use previous segments if requested year_frac is not belong to a segment.
  var coefNames = ["INTP", "SLP", "COS", "SIN", "COS2", "SIN2", "COS3", "SIN3"]
  coefNames = coefNames.map(function(v) {
    return band + '_' + v.toLowerCase()
  })

  //use next segment
  var vicNext = tEnd.gt(year_frac)
  //is year_frac pass the end of time series
  var fixNext = vicNext.arrayReduce(ee.Reducer.sum(), [0]).arrayFlatten([['count']])
  var selectedNext = coefs.arrayMask(vicNext.toArray(1)).arraySlice(0, 0, 1).arrayProject([1]).arrayFlatten([coefNames])

  //use previous segment
  var vicPrev = tStart.lte(year_frac)
  //is year_frac before the beginning of time series
  var fixPrev = vicPrev.arrayReduce(ee.Reducer.sum(), [0]).arrayFlatten([['count']])
  var selectedPrev = coefs.arrayMask(vicPrev.toArray(1)).arraySlice(0, -1).arrayProject([1]).arrayFlatten([coefNames])

  //fix beinning and end of time series issue
  var fixedSelectedNext = selectedNext.where(fixNext.eq(0), selectedPrev)
  var fixedSelectedPrev = selectedPrev.where(fixPrev.eq(0), selectedNext)

  return useNext ? fixedSelectedNext : fixedSelectedPrev;
}

var getSegmentParams = function(ccdc_img, year_frac, bands, normalize, useNext) {
  if (normalize) {
    ccdc_img = normalizeIntercept(ccdc_img, bands)
  }

  var result = []
  bands.forEach(function(band) {
    result.push(getSegmentParamsByBandImpl(ccdc_img, year_frac, band, useNext))
  })

  return ee.Image(result).float()
}

var coefs_img_zy = getSegmentParams(ccdc_img, 2015.5, bands, true, true)
print('CCDC multiband coefficient image:', coefs_img_zy)
Map.addLayer(coefs_img_zy, null, 'CCDC coefficients (ZY)', 0)

// get phase and amplitude associated with a time/segment
function getPhaseAmplitude(coefs_img, sinExpr, cosExpr){
    var sin = coefs_img.select(sinExpr)
    var cos = coefs_img.select(cosExpr)

    var phase = sin.atan2(cos)
      // Scale to [0, 1] from radians.
      .unitScale(-3.14159265359, 3.14159265359)
      .multiply(365) // To get phase in days!

    var amplitude = sin.hypot(cos)
    var ap_ratio = amplitude.divide(phase)

    var phaseNames = phase.bandNames().map(function(x){return ee.String(x).replace('_sin', '_phase')})
    var amplitudeNames = amplitude.bandNames().map(function(x){return ee.String(x).replace('_sin', '_amplitude')})
    var apNames = ap_ratio.bandNames().map(function(x){return ee.String(x).replace('_sin', '_amphratio')})
    return phase.rename(phaseNames).addBands(amplitude.rename(amplitudeNames)).addBands(ap_ratio.rename(apNames))
}
var phase_amp_img_zy = getPhaseAmplitude(coefs_img_zy, '.*sin.*', '.*cos.*')
print('CCDC multiband phase+amplitude image:', phase_amp_img_zy)
Map.addLayer(phase_amp_img_zy, null, 'CCDC phase+amp (ZY)', 0)


// get synthetic image (i.e. predictions of spectra/indices at a given time)
// TODO: possible to extrapolate a segment? KO looking into Wiell
var getSyntheticByBandImpl = function(ccdc_img, year_frac, band, useNext) {
  var omega = 2.0 * Math.PI
  var year_frac = ee.Number(year_frac)
  var imageT = ee.Image.constant([1, year_frac,
                                year_frac.multiply(omega).cos(),
                                year_frac.multiply(omega).sin(),
                                year_frac.multiply(omega * 2).cos(),
                                year_frac.multiply(omega * 2).sin(),
                                year_frac.multiply(omega * 3).cos(),
                                year_frac.multiply(omega * 3).sin()])

  var params = getSegmentParamsByBandImpl(ccdc_img, year_frac, band, useNext)
  return imageT.multiply(params).reduce('sum').rename(band)
}
var getSynthetic = function(ccdc_img, year_frac, bands, useNext) {
  var result = []
  bands.forEach(function(b) {
    result.push(getSyntheticByBandImpl(ccdc_img, year_frac, b, useNext))
  })
  return ee.Image(result)
}

var syn_zy = getSynthetic(ccdc_img, year_frac, bands, true)
print('CCDC synthetic image:', syn_zy)
Map.addLayer(syn_zy, lsat_fcolor_viz, 'Synthetic (ZY)')


// breaks (change) info
// number of breaks
function get_n_breaks_img(ccdc_img, start_year, end_year){
  var tBreak = ccdc_img.select('tBreak')
  var vic = tBreak.lt(ee.Number(end_year)).and(tBreak.gte(ee.Number(start_year))).toArray(0)
  return tBreak.arrayMask(vic).arrayLength(0)
}

// number of break greater than an absolute magnitude threshold
function get_n_breaks_thresh_img(ccdc_img, start_year, end_year, band_name, mag_thresh_abs){
  var band_mag_name = ee.String(band_name).cat(ee.String('_magnitude'))
  var band_mag = ccdc_img.select(band_mag_name)
  var tBreak = ccdc_img.select('tBreak')
  var vic = tBreak.lt(ee.Number(end_year))
                  .and(tBreak.gte(ee.Number(start_year))
                  .and(band_mag.abs().gte(mag_thresh_abs))).toArray(0)
  return tBreak.arrayMask(vic).arrayLength(0)
}
var n_breaks_thresh_img = get_n_breaks_thresh_img(ccdc_img, start_year, end_year, 'NDVI', 0.2)
Map.addLayer(n_breaks_thresh_img.selfMask(), {palette: palettes.matplotlib.viridis[7], min:1, max: 10}, 'N breaks above threshold')

// year of first break
function get_first_tbreak_img(ccdc_img, start_year, end_year){
  var tBreak = ccdc_img.select('tBreak')
  var vic = tBreak.lt(ee.Number(end_year)).and(tBreak.gte(ee.Number(start_year))).toArray(0)
  return tBreak.arrayMask(vic).arraySlice(0, 0, 1).arrayProject([0]).arrayFlatten([['tBreak']])
}

// year of largest break greater than an absolute magnitude threshold
function get_largestmag_tbreak_thresh_img(ccdc_img, start_year, end_year, band_name, mag_thresh_abs){
  var band_mag_name = ee.String(band_name).cat(ee.String('_magnitude'))
  var band_mag = ccdc_img.select(band_mag_name)
  var tBreak = ccdc_img.select('tBreak')
  var mask = tBreak.lt(ee.Number(end_year))
                   .and(tBreak.gte(ee.Number(start_year))
                   .and(band_mag.abs().gte(mag_thresh_abs))).toArray(0)
  var dates = tBreak.arrayMask(mask).arrayPad([1]);
  var magnitudes = ccdc_img.select(band_mag_name)
						               .arrayMask(mask)
						               .arrayPad([1]);

	var maxIndex = magnitudes.abs()
						               .arrayArgmax()
						               .arrayFlatten([['index']]);
	var sel_mag = magnitudes.arrayGet(maxIndex);
	var sel_tBreak = dates.arrayGet(maxIndex).selfMask();

	return sel_tBreak
}

var largest_mag_thresh = get_largestmag_tbreak_thresh_img(ccdc_img, start_year, end_year, 'NDMI', 0.2)
Map.addLayer(largest_mag_thresh, {palette: palettes.matplotlib.inferno[7], min: start_year, max: end_year}, 'Year of Largest Change')

var mtbs_fires = ee.FeatureCollection("USFS/GTAC/MTBS/burned_area_boundaries/v1")
Map.addLayer(mtbs_fires.style({color: 'red', fillColor:'ffffff00'}), null, 'Fire Polys')

throw('e')
////
// USFS GTAC methods
// see: https://geeviz.org/
////
// get a multiband coefficient image associated with a time/segment
function getCCDCSegCoeffs(coefs_img, year_frac, fillGaps) {
  var coeffKeys = [".*_coefs"];
  var tStartKeys = ["tStart"];
  var tEndKeys = ["tEnd"];
  var tBreakKeys = ["tBreak"];

  //Get coeffs and find how many bands have coeffs
  var coeffs = coefs_img.select(coeffKeys);
  var bns = coeffs.bandNames();
  var nBns = bns.length();
  var harmonicTag = ee.List(["INTP", "SLP", "COS1", "SIN1", "COS2", "SIN2", "COS3", "SIN3"]);

  //Get coeffs, start and end times
  coeffs = coeffs.toArray(2);
  var tStarts = coefs_img.select(tStartKeys);
  var tEnds = coefs_img.select(tEndKeys);
  var tBreaks = coefs_img.select(tBreakKeys);

  //If filling to the tBreak, use this
  tStarts = ee.Image(ee.Algorithms.If(fillGaps, tStarts.arraySlice(0, 0, 1).arrayCat(tBreaks.arraySlice(0, 0, -1), 0), tStarts));
  tEnds = ee.Image(ee.Algorithms.If(fillGaps, tBreaks.arraySlice(0, 0, -1).arrayCat(tEnds.arraySlice(0, -1, null), 0), tEnds));

  var timeImg = ee.Image(ee.Number(year_frac))
  //Set up a mask for segments that the time band intersects
  var tMask = tStarts.lt(timeImg).and(tEnds.gte(timeImg)).arrayRepeat(1, 1).arrayRepeat(2, 1);
  coeffs = coeffs.arrayMask(tMask).arrayProject([2, 1]).arrayTranspose(1, 0).arrayFlatten([bns, harmonicTag]);

  //If time band doesn't intersect any segments, set it to null
  coeffs = coeffs.updateMask(coeffs.reduce(ee.Reducer.max()).neq(0));

  return coeffs;
}
var coefs_img_gtac = getCCDCSegCoeffs(ccdc_img, year_frac, true)

// get synthetic image (i.e. predictions of spectra/indices at a given time)
function get_syn_img(coefs_img, year_frac, whichHarmonics, whichBands) {
  //Unit of each harmonic (1 cycle)
  var omega = ee.Number(2.0).multiply(Math.PI);

  //Pull out the time band in the yyyy.ff format
  var tBand = ee.Image(ee.Number(year_frac))

  //Pull out the intercepts and slopes
  var intercepts = coefs_img.select([".*_INTP"]);
  var slopes = coefs_img.select([".*_SLP"]).multiply(tBand);

  //Set up the omega for each harmonic for the given time band
  var tOmega = ee.Image(whichHarmonics).multiply(omega).multiply(tBand);
  var cosHarm = tOmega.cos();
  var sinHarm = tOmega.sin();

  //Set up which harmonics to select
  var harmSelect = whichHarmonics.map(function (n) {
    return ee.String(".*").cat(ee.Number(n).format());
  });

  //Select the harmonics specified
  var sins = coefs_img.select([".*_SIN.*"]);
  sins = sins.select(harmSelect);
  var coss = coefs_img.select([".*_COS.*"]);
  coss = coss.select(harmSelect);

  //Set up final output band names
  var outBns = whichBands.map(function (bn) {
    return ee.String(bn).cat("_CCDC_fitted");
  });

  //Iterate across each band and predict value
  var predicted = ee
    .ImageCollection(
      whichBands.map(function (bn) {
        bn = ee.String(bn);
        return ee.Image([intercepts.select(bn.cat("_.*")), slopes.select(bn.cat("_.*")), sins.select(bn.cat("_.*")).multiply(sinHarm), coss.select(bn.cat("_.*")).multiply(cosHarm)])
                 .reduce(ee.Reducer.sum());
      })
    )
    .toBands()
    .rename(bands);
  return predicted;
}
var syn_gtac = get_syn_img(coefs_img_gtac, year_frac, [1,2,3], bands)
Map.addLayer(syn_gtac, lsat_fcolor_viz, 'Synthetic (GTAC)', 0)