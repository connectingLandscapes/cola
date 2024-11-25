## Description

Creates a continuous surface of cumulative paths that connects a set of points. crk.py in the command line, or crk_py( ) in R.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :----: |:--------------:  | :----: | :-----------------: | :----:  |:---------- |
| Source points |      | inshp|       | String|      |File path to the point layer. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file|
| Surface resistance |   |intif|     | String|      | File path to the input raster. Requires a projected file with square pixels. Not LonLat projection allowed|
| Output raster filename|  |outtif|  | string| |File path of the output raster layer written in GeoTiff format.|
|Max. dispersal distance in cost units| |maxdist| | Numeric| | This is the maximum distance to consider when calculating kernels and should correspond to the maximum dispersal distance of the focal species. Values greater than this will be converted to 0 before summing kernels. For example, if the maximum dispersal distance of the focal species is 10 km, set this value to 10000.|
|Kernel shape| |shape| | Numeric| | String. This determines how the probability of dispersal declines with distance from the focal point. 'linear' implements the function 1 - (1/dThreshold) * d * where dThreshold is the specified distance threshold and d is the distance from the focal point. 'gaussian' implements the function exp(-1*((d^2)/(2*(dispScale^2)))) where d is the distance from the focal point and dispScale is equal to dThreshold/4. In future versions, users will be able to specify dispScale. |
|Use transformation?  | |transf|  |String|  | ’yes’ or ‘no’. If it is set to ‘no’, then no kernel volume transformation is applied and the next argument [8] is ignored. If it is set to ‘yes’, then the value in argument volume [8] is used to transform the kernel volume.|
|Kernel volume |  | volume|  |Numeric|  | If 1, the default, the resistant kernel value at the origin is 1 and no kernel volume transformation is applied. If > 1, the parameter value is used to scale distance values by a constant that is determined by the equation kVol * 3/(pi*dThreshold^2) where kVol is the kernel volume parameter, dThreshold is the specified distance threshold, and pi is the mathematical constant pi. The constant is then multiplied by the distances to the focal point resulting in a scaled kernel volume.|
|Number of cores | | ncores| |Numeric| |Number or cores of your computer to run the analysis | 
|Projection string  | | crs| |String| |Projection information in the case the input raster [2] has no spatial projection. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.|