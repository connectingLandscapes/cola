## Description

Identify and rank existing corridors, assigning a higher rank to the shorter links that connect the bigger nodes. You can use it by prioritize_core_conn.py in the command line or pri_py( ) in R.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :--: |:--------------:  | :----: | :-----------------: | :--: |:---------- |
|Surface resistance| |tif| |String| |File path to the input raster. Requires a GeoTIFF file with square pixels and cost units as distance units|
|Kernels| |incrk| |String| |File path to the input kernel GeoTiff raster file.| 
|Corridors| |inlcc| |String| |File path to the input corridors GeoTiff raster file.| 
|Masked output raster| |maskedcsname| |String| |The output raster file path of the resulting nodes, derived from the kernels raster and the quantile ([10] threshold)| 
|Output point layer| |outshppoint| |String| |The output point layer file path of the resulting corridor centroids in ESRI Shapefile format.|
|Output polygon layer| |outshppol| |String| |The output polygon layer file path of the resulting corridors in ESRI Shapefile format.|
|Output patch polygon layer| |outshppatch| |String| |The output polygon layer file path of the resulting patches in ESRI Shapefile format.|
|Output patch raster layer| |outtifpatch| |String. The output raster layer file path of the resulting patches in GeoTIFF format.|
|Output corridor raster layer| |outtif| |String| |The output raster layer file path of the resulting corridors in GeoTIFF format.|
|Threshold quantile| |threshold| |Decimal| |Numeric value between zero and one (0 - 1) to convert continuous kernels into discrete patches.|
|Maximum distance| |maxdist| |Numeric| |This is the maximum distance to consider when calculating kernels and should correspond to the maximum dispersal distance of the focal species. Values greater than this will be converted to 0 before summing kernels. For example, if the maximum dispersal distance of the focal species is 10 km, set this value to 10000. This should be the same value used for the parameter [4], Max. dispersal distance in cost units in the least cost path function.|
