## Description

Compares kernels (crk) or corridors (lcc) baseline layers with other scenarios. You can use them by crk_compare.py and lcc_compare.py in the command line or crk_compare_py( ) and lcc_compare_py( ) in R.

###  Parameters

| Variable          | .  _  . | Name            | .  _  .  |   Type | .  _  . | Description |
| :---------------: | :--: |:--------------:  | :----: | :-----------------: | :--: |:---------- |
| Baseline raster| |intif| |String| |File path to the input baseline kernel or corridor raster. Requires a GeoTIFF file with square pixels and cost units as distance units|
|Layers to compare| |intifs| |String| |String of raster layer paths to compare. Requires to be passed as a comma separated single string as “layer1.tif, layer2.tif,layer3.tif,...”. Accepts relative or full paths. Pass the baseline raster again [1] as the first raster.|
|Absolute difference out table| |outcsvabs| |String| |File path of the table to be written that contains the absolute values of the layers. Returns one value for argument [1] plus one of each of the layers passed on the argument [2]|
|Relative difference out table| |outcsvrel| |String| |File path of the table to be written that contains the relative values of the layers. Returns one value for each of the layers passed on the argument [2] relative to the argument [1]|
|Absolute difference out png| |outpngabs| |String| |File path of the barplot image to be written that contains the absolute values of the layers. Returns one value for argument [1] plus one of each of the layers passed on the argument [2]|
|Relative difference out png| |outpngrel| |String| |File path of the barplot image to be written that contains the relative values of the layers. Returns one value for each of the layers passed on the argument [2] relative to the argument [1]|
|Out folder path| |outfolder| |String| |Data folder path where the comparison rasters will be saved. Writes as many layers as passed in the argument [2]|
|Cutline vectorial layer path| |inshp| |String| |File path to a vectorial layer used to cutline the comparison. Default value is 'None' for ignoring this parameter.| 
|Cutline vectorial layer field| |shpfield| |String| |Attribute or column name of the vectorial layer provided if you want to regionalize the comparison. By providing a valid column name the csv table will have a second dimension of the length of the unique values of the column. Default value is 'None'.|
