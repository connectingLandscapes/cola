###   `cola` functions

`cola` contains the following functions:
  * A. Habitat suitability (HS) to surface resistance (SR)
  * B. Create source points
  * C. Cost distance matrix
  * D. Least corridors paths
  * E. Cumulative resistance dispersal kernels
  * F. Prioritization
  * H. CDPOP -- soon
  * G. Landscape genetics -- soon

-----------

#####  **A. Habitat suitability (HS) to surface resistance (SR)**

Converts a a given raster into a new given transformation function and outputs a numeric range.
  
  **Python file:** s2res.py
  
  **Parameters:**
  
[A] py: Python executable path or instruction R runs in console to call cola python.

[B] src: Python script path

[1] initif: Habitat Suitability (HS) georreferenced TIF raster

[2] outif: Surface resistance (SR) output file name TIF raster

[3] param3: Suitability grid min value. Numeric value to cutoff the original layer. 

[4] param4: Suitability grid max value. Numeric value to cutoff the original layer. 

[5] param5: Maximum resistance value to have the new raster. Min value is set to 1.

[6] param6: Shape Parameter. Between 1 and X

[7] param7: No data value of suitability raster. Default NULL. To be guessed from metadata if not provided

[8] param8: CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"


  **Usage in R:**
```
library(cola)
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_tif <- 'outputsurfaceresistance.tif'

surface_resistance <- s2res_py(intif = input_tif, outtif = out_tif,
   param3 = 0
   param4 =  100
   param5 = 100
   param6 = 1)
  
if(file.exists(out_tif)){
  terra::plot(terra::rast(out_tif))  
}
```

  **Usage in console:**
```
# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe s2res.py myhabitatsuitability.tif outputsurfaceresistance.tif 0 100 100 1

# Linux:
/home/user/anaconda3/envs/cola/bin/python s2res.py myhabitatsuitability.tif outputsurfaceresistance.tif 0 100 100 1
```

**Creates:**  *outputsurfaceresistance.tif*, a raster layer with a minimum value of 1 and maximum value given the parameter [5]

**Common errors:**
Add errors here

-----------


#####  **B. Simulate source points**

Creates a vectorial layer with a given number of points. The points are generated , Converts a a given raster into a new given transformation function and outputs a numeric range.

**Python file:** create_source_points.py


  **Parameters:**
[A] py: Python executable path or instruction R runs in console to call cola python.
[B] src: Python script path
[1] initif: Reference georreferenced TIF raster layer
[2] outshp:  Vectorial output layer
[3] param3: Minimum value to mask the raster layer
[4] param4: Maximun value to mask the raster layer
[5] param5: Number of source points to create
[6] param6: CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"

  **Usage in R:**
```
library(cola)
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_shp <- 'outputPoints.shp'

simulated_points <- points_py(intif = input_tif, outshp = out_shp,
   param3 = 0 param4 =  100 param5 = 50)
  
if(file.exists(out_tif)){
  plot(sf::read_sf(simulated_points))  
}
```

  **Usage in console:**
```

# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe create_source_points.py sampleTif.tif outputPoints.shp 0 100 50

# Linux:
/home/user/anaconda3/envs/cola/bin/python create_source_points.py sampleTif.tif outputPoints.shp 0 100 50
```

**Creates:**  *outputPoints.shp*, a vectorial layer with 100 points located over the pixels between 1 and 100

**Common errors:**

Add errors here


-----------

#####  **C. Cost distance matrix**

Creates a two way interaction matrix path from a set of points and a surface resistance raster layer


  **Python file:** create_cdmat.py
  
  **Parameters:**
[A] py: Python executable path or instruction R runs in console to call cola python.
[B] src: Python script path
[1] inshp: Vectorial layer with points file path
[2] initif: Georreferenced TIF raster file path
[3] outcsv: Resulting two-way matrix output file path
[4] param4: Distance threshold (in cost distance units)
[5] param5: Number of processors. Default is one (1)
[6] param6: User provided CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028". Default 'None'

  **Usage in R:**
```
library(cola)
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
input_shp <- system.file(package = 'cola', 'sampledata/SourcePoints_50.shp')
out_csv <- 'outputCSVmatrix.csv'

out_matrix <- cdmat_py(inshp = input_shp, intif = input_tif, 
 outcsv = out_csv, param3 = 25000)
  
if(file.exists(out_matrix)){
  head(read.csv(out_matrix))  
}
```

  **Usage in console:**
```
# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe create_cdmat.py /path/to/sampleTif.tif /path/to/SourcePoints_50.shp /path/to/outputCSVmatrix.csv 500 1 'None'

# Linux:
/home/user/anaconda3/envs/cola/bin/python create_cdmat.py /path/to/sampleTif.tif /path/to/SourcePoints_50.shp /path/to/outputCSVmatrix.csv 500 1 'None'
```

**Creates:**  *outputCSVmatrix.csv*, a CSV two-way matrix file with the cumulative resistance among each pair of points that occurs below param4

**Common errors:**
Add errors here

-----------


-----------

#####  **D. Least cost corridor paths (LCC)**

Creates a continuous raster with the cumulative corridor (distance * value) bewteen the points


  **Python file:** lcc.py and lcc_heavy.py for bigger rasters
  
  **Parameters:**
[A] py: Python executable path or instruction R runs in console to call cola python.
[B] src: Python script path
[1] inshp: Spatial source point layer (any ORG driver), CSV (X, Y files), or *.xy file
[2] initif: Habitat Suitability (HS) georreferenced TIF raster
[3] outif: Least cost path corridors (LCC) output file name TIF raster
[4] param4: Distance threshold to connect points. It should be in map distane units (such as meters)
[5] param5: Corridor radius smoothing factor in number of cells. 
[6] param6: Corridor tolerance in cost distance units. This is the amount to add to the least cost path to create least cost corridors (should be in meters)
[7] param7: Number of cores. Default is one (1)
[8] param8: CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"

  **Usage in R:**
```
library(cola)
input_shp <- system.file(package = 'cola', 'sampledata/SourcePoints_50.shp')
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_tif <- 'outputcorridors.tif'


# Use lcc_py( ) or lccHeavy_py( ) in the same way
corridors <- lcc_py( inshp = input_shp, intif = input_tif, outtif = out_tif,
   param4 =  25000
   param5 = 100
   param6 = 5000)
  
if(file.exists(corridors$file)){
  terra::plot(terra::rast(corridors$file))  
}
```

  **Usage in console:**
```
# use lcc.py or lcc_heavy.py

# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe lcc.py SourcePoints_50.shp  myhabitatsuitability.tif outputcorridors.tif 25000 5 5

# Linux:
/home/user/anaconda3/envs/cola/bin/python lcc.py SourcePoints_50.shp  myhabitatsuitability.tif outputcorridors.tif 25000 5 5

```

**Creates:**  *outputcorridors.tif*, a raster layer with continuous values of cumulative resistance between points that fail below distance threshold 

**Common errors:**
Raster has no define NODATA values
Raster has no or uncommon projection
Raster has no squared pixels
Raster has multiple bands


  

-----------

#####  **E. Create connectivity dispersal resistance kernels (CRK)**

Creates a continuous raster with the cumulative resistance kernels arround points

  **Python file:** crk.py
  
  **Parameters:**
[A] py: Python executable path or instruction R runs in console to call cola python.
[B] src: Python script path
[1] inshp: Spatial source point layer (any ORG driver), CSV (X, Y files), or *.xy file
[2] initif: Habitat Suitability (HS) georreferenced TIF raster
[3] outif: Least cost path corridors (LCC) output file name TIF raster
[4] param4: Distance threshold to connect points. It should be in map distane units (such as meters)
[5] param5: Kernel shape. Accepted: linear, gaussian
[6] param6: Kernel volume. Set to 1 if no kernel volume multiplier is desired.
[7] param7: Number of cores. Default is one (1)
[8] param8: CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"

  **Usage in R:**
```
library(cola)
input_shp <- system.file(package = 'cola', 'sampledata/SourcePoints_50.shp')
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_tif <- 'outputkernels.tif'


kernels <- crk_py( inshp = input_shp, intif = input_tif, outtif = out_tif,
   param4 =  25000, 
   param5 = 'linear'
   param6 = 1)
  
if(file.exists(kernels$file)){
  terra::plot(terra::rast(kernels$file))  
}
```

  **Usage in console:**
```
# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe crk.py SourcePoints_50.shp  myhabitatsuitability.tif outputkernels.tif 25000 linear 1

# Linux:
/home/user/anaconda3/envs/cola/bin/python crk.py SourcePoints_50.shp  myhabitatsuitability.tif outputkernels.tif 25000 linear 1

```

**Creates:**  *outputkernels.tif*, a raster layer with continuous values of cumulative resistance kernels around the points

**Common errors:**
Raster has no define NODATA values
Raster has no or uncommon projection
Raster has no squared pixels
Raster has multiple bands

-----------




-----------

#####  **F. Prioritization**

Prioritize the resulting corridors (lcc) based on the distance and costs of the cumulative kernels (crk)

  **Python file:** prioritize_core_conn.py
  
  **Parameters:**
[A] py: Python executable path or instruction R runs in console to call cola python.
[B] src: Python script path
[1] initif: Habitat Suitability (HS) georreferenced TIF raster. Mandatory.
[2] incrk: Cumulative resistance kernels input layer. Mandatory.
[3] inlcc: Cumulative resistance kernels input layer. Mandatory.
[4] maskedcsname: Output masked cost surface (could be a temp file)
[5] outshp: Output vectorial centroids of the prioritized corridors. Mandatory.
[6] outtif: Output raster of the prioritized corridors. Mandatory.
[7] param7: Quantile threshold for identifying high value patches from resistant kernel surface. Should be between 0 and 1. Mandatory.
[8] param8: Corridor tolerance. Should be the same tolerance used in the lcc script/function

  **Usage in R:**
```
library(cola)
input_shp <- system.file(package = 'cola', 'sampledata/SourcePoints_50.shp')
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_tif <- 'outputkernels.tif'


kernels <- crk( inshp = input_shp, intif = input_tif, outtif = out_tif,
   param4 =  25000, 
   param5 = 'linear'
   param6 = 1)
  
if(file.exists(kernels$file)){
  terra::plot(terra::rast(kernels$file))  
}
```

  **Usage in console:** 
```
# Windows:
C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe crk.py SourcePoints_50.shp  myhabitatsuitability.tif outputkernels.tif 25000 linear 1

# Linux:
/home/user/anaconda3/envs/cola/bin/python crk.py SourcePoints_50.shp  myhabitatsuitability.tif outputkernels.tif 25000 linear 1

```

**Creates:**  *outputkernels.tif*, a raster layer with continuous values of cumulative resistance kernels around the points

**Common errors:**
Raster has no define NODATA values
Raster has no or uncommon projection
Raster has no squared pixels
Raster has multiple bands

-----------
