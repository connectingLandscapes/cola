/*
By: Patrick Burns [pb463@nau.edu], Northern Arizona University

About: export Landsat Col. 2 CCDC results corresponding to a grid

*/



// ----------------------
// ----- IMPORTS -----
// ----------------------

// Landsat processing module developed by the GEODE lab
var lsat_proc = require('users/pb463/NAU_GoetzGroup:01_utilities/modules/landsat_proc_mod')



// ----------------------
// ----- INPUTS -----
// ----------------------

// a geometry for the region of interest
var region_geom = ee.Feature(ee.FeatureCollection('users/pb463/S2L/Sonoma_cty_v2_PBcleaned').first())
                    .geometry()
Map.addLayer(region_geom)

// a grid covering the AOI
var proj_epsg = 'EPSG:3310'
var grid = region_geom.coveringGrid({proj: proj_epsg, scale: 10000})

var ccdc_start = '1984-01-01'
var ccdc_end = '2024-12-31'
var max_cloud_cover_land = 50
var min_sun_elev = 30
var bands = ['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'NDVI', 'NDMI', 'NBR2']
var sensors_meta = 'landsat_45789'

// export asset folder
var asset_folder = 'projects/ee-gedibio/assets/ccdc/results/workshop_sonoma/'



// ----------------------
// ----- PROCESSING -----
// ----------------------

// Landsat Collection 2 image collections
var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")

// Merge landsat collections
var l4to9 = l5.merge(l7).merge(l8).merge(l9)

// Filter the merged collection spatially and temporally.
// Mask low/high SR, exclude clouds, and calculate spectral indices
var l4to9_p = l4to9.filterBounds(region_geom)
                   .filterDate(ccdc_start, ccdc_end)
                   .filter(ee.Filter.gte('SUN_ELEVATION', min_sun_elev))
                   .filter(ee.Filter.lte('CLOUD_COVER_LAND', max_cloud_cover_land))
                   .map(lsat_proc.l4to9_c2_rename_bands)
                   .map(lsat_proc.l4to9_c2_scaleoff)
                   .map(lsat_proc.l4to9_c2_maskSR)
                   .map(lsat_proc.l4to9_c2_qa_maskClouds)
                   .map(lsat_proc.l4to9_c2_indices)
                   .select(bands)

//print('Landsat merged, filtered, and masked collection (first 100):', l4to9_p.limit(100))

// ccdc input params
var ccdcParams = {
  collection: l4to9_p,
  breakpointBands: ['green', 'red', 'NIR', 'SWIR1', 'SWIR2'],
  tmaskBands: ['green', 'SWIR1'],
  minObservations: 6,
  chiSquareProbability: 0.99,
  minNumOfYearsScaler: 1.33,
  dateFormat: 1,
  lambda: 0.002,
  maxIterations: 10000}

// Run CCDC
var ccdc =  ee.Algorithms.TemporalSegmentation.Ccdc(ccdcParams)

// Show the structure of the ccdc image
//print('CCDC structure:', ccdc)



// ----------------------
// ----- EXPORT -----
// ----------------------

// Create a metadata dictionary with the parameters and arguments used.
  var metadata=ccdcParams;
   metadata['breakpointBands']=metadata['breakpointBands'].toString();
   metadata['tmaskBands']=metadata['tmaskBands'].toString();
   metadata['startDate']=ccdc_start;
   metadata['endDate']=ccdc_end;
   metadata['min_sun_elev']=min_sun_elev;
   metadata['max_cloud_cover_land']=max_cloud_cover_land;
   metadata['bands']=bands.toString();
   metadata['sensors']=sensors_meta;

  // Export multiple tiled images using a pre-computed grid
  Map.addLayer(grid, {'color': 'green'}, 'grids for CCDC export')

  // Convert the grid to a list and loop through tiles for export
  var grid_list = grid.toList(1000)

  grid_list.size().evaluate(function(s){
    print('Number of square grids: ', s)
    for (var i = 0; i<s; i++){
      var outGeo = ee.Feature(grid_list.get(i)).geometry()
      Export.image.toAsset({
        image: ccdc.set(metadata),
        scale: 30,
        crs: proj_epsg,
        description: 'ccdc-lsc2_' + ccdc_start + '_' + ccdc_end + '_sq-grid-' + i,
        maxPixels: 1e13,
        region: outGeo,
        assetId: asset_folder + 'ccdc-lsc2_' + ccdc_start + '_' + ccdc_end + '_grid-' + i ,
        pyramidingPolicy: { '.default': 'sample' } })
    }
  })
