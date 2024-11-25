## Description

Creates a continuous surface of cumulative paths that connects a set of points. lcc.py in the command line or lcc_py( ) in R. Also, you can use two extra versions of this function to generate equivalent results. The first is a memory-safe version for big files: lccHeavy_py( ) and lcc_heavy.py. The second one is for parallel computing: lccJoblib_py( ) and lcc_joblib.py. lccjoblib.py and lccJoblib_py will work best if you read and write results to your primary hard drive, usually your C: drive. Your primary hard drive should be SSD (most are) and should have at least 50GB of free space (more if you are processing very large landscapes with many points). You can make a rough estimate of free space needed according to the number of cells in your landscape x the number of points x 64 bit x the compression ratio x conversion factor to GB. For example, a 10 million cell landscape with 5000 points will require ~80GB of temporary storage: ((10000000*5000*64)*0.2)/8e9 = 80.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :--: |:--------------:  | :----: | :-----------------: | :--: |:---------- |
| | | | | | | | |
| Source points |      | inshp|       | String|      |File path to the point layer. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file|
| Surface resistance |   |intif|     | String|      | File path to the input raster. Requires a projected file with square pixels. Not LonLat projection allowed|
| Output raster filename|  |outtif|  | string| |File path of the output raster layer written in GeoTiff format.|
|Max. dispersal distance in cost units| |maxdist| | Numeric| | This is the maximum distance to consider when calculating corridors and should correspond to the maximum dispersal distance of the focal species. For example, if the maximum dispersal distance of the focal species is 10 km, set this value to 10000. Values greater than this will be converted to 0 before summing corridors.|
| | | | | | | | |
|Corridor smoothing factor| |smooth| | Numeric| |The width of the window, in the number of cells, is used to smooth the output corridor surface. If no smoothing is desired, set it to 0. This parameter allows backward compatibility with the original UNICOR functionality, which runs a smoothing window over the least-cost path surface. |
| Corridor tolerance in cost units  | |tolerance|  |String|  | This is the distance beyond the least-cost path that an animal might traverse when moving between source points. Larger values result in wider corridors.|
|Number of cores | | ncores| |Numeric| |PNumber or cores of your computer to run the analysis | 
|Projection string  | | crs| |String| |Projection information in the case the input raster [2] has no spatial projection. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.|
| First temporary h5 file (not used in R)| |NA| |String| |Not used in R, but required in Python. Created as a temporal file internally in R but required in the Python command line. This file is created to store temporal structures in disk and avoid allocate the information on the RAM. This argument is used only in the R and Python version of lcc_heavy( ) and lcc_joblib( ).|
| Second temporary h5 file (not used in R)| | NA | String| |Not used in R, but required in Python. Created as a temporal file internally in R but required in the Python command line. This file is created to store temporal structures in disk and avoid allocate the information on the RAM. This argument is used only in the R and Python version of lcc_heavy( ) and lcc_joblib( ).| 
|Max GB ram allowed| |maxram| |Numeric| |Maximum gibabytes to use in your machine. The default value is 6.|