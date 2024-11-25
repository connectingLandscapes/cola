## Description

Creates a dispersal resistance matrix among a set of points. create_cdmat.py in the command line or cdmat_py( ) in R.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :--: |:--------------:  | :----: | :-----------------: | :--: |:---------- |
| | | | | | | | |
|Source points |      | inshp|       | String|      |File path to the point layer. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file|
| Output point layer|   |outtif|  | String| |File path to the output point layer. Written in ESRI Shapefile format.|
| Surface resistance |   |intif|     | String|      | File path to the input raster. Requires a projected file with square pixels. Not LonLat projection allowed|
| | | | | | | | |
|Minimum value| |minval| | Numeric| | The lower value of the pixels in the raster to consider to simulate the points.
|Maximum value| |maxval| | Numeric| | The upper value of the pixels in the raster to consider to simulate the points. 
|Number of points | |npoints|  |Numeric|  | Number of points to simulate.
|Is it suitable?| | | issuit|  |String|  |‘Yes’ (default) or ‘No’. Indicates if the provided raster [1]  is suitability. If so, the script will likely sample higher value pixels. If ‘No’, will assume it is resistance and will sample more likely lower values|
|Update CRS | | upcrs| |String| |Projection information in the case the input raster [1] has no spatial projection. For GeoTiffs, this is automatically determined. For text-based files like ASCII or RSG rasters, this must be input by the user. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.