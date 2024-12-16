## Description

Creates a dispersal resistance matrix among a set of points. create_cdmat.py in the command line or cdmat_py( ) in R.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :--: |:--------------:  | :----: | :-----------------: | :--: |:---------- |
| | | | | | | | |

|Source points |      | inshp|       | String|      |File path to the point layer. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file|
| Surface resistance |   |intif|     | String|      | File path to the input raster. Requires a projected file with square pixels. Not LonLat projection allowed|
| | | | | | | |
|Output csv file name |   | outcsv|      | String|      | File path of the outputcsv matrix|
| | | | | | | |
|Max. dispersal distance in cost units|   |maxdist|   |Numeric|   |This is the maximum resistance value after transformation from suitability. The default value is 100. The minimum value is set to 1.|
|Number of cores|   | ncores|      | Numeric|      | Number of CPU cores to run the analysis|
|Projection parameter|  |prj|  |String|  |Projection information in the case the input raster  has no spatial projection. For GeoTiffs, this is automatically determined. For text-based files like ASCII or RSG rasters, the user must input them. Provide it as EPSG or ESRI string, e.g. "ESRI:102028". The default value is ‘None’|


