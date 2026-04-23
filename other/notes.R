"

############################ installing #################################
detach('package:cola', unload=TRUE)
.rs.restartR()
remove.packages('cola')
.rs.restartR()
if(!require(devtools)){install.packages('devtools')}
library(devtools)
.rs.restartR()
devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') ## option 3: None
devtools::install('N:/My Drive/git/cola')
"

" # command
su shiny;
cd /home/shiny/cola; git pull
R -e "devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never')"
sudo rm /srv/shiny-server/connecting-landscapes -R; sudo cp /home/shiny/cola/inst/app/* /home/shiny/colashiny/connecting-landscapes/. -R
R -e "shinyParallel::installShinyParallel('/home/shiny/colashiny/connecting-landscapes', max.sessions = 25)"
"

"
library(cola)
cola::diagnose_cola()
cola::setup_cola()
cola::setup_cola(ask = FALSE, zarr = TRUE, dss = TRUE, cola2 = TRUE)
file.edit(file.path(Sys.getenv('HOME'), '.Renviron'))
.rs.restartR()
Sys.getenv(c('COLA_PYTHON_PATH', 'COLA_SCRIPTS_PATH'))
Sys.getenv(c('COLA_DATA_PATH', 'COLA_SCRIPTS_PATH', 'COLA_DSS_UPL_MB', 'COLA_VIZ_THRES_PIX', 'COLA_VIZ_RES_NCOL', 'COLA_VIZ_RES_NROW', 'COLA_NCORES' ))
cola::setup_cola_dss()
.rs.restartR()

library(cola)
cola::cola_dss()
# install.packages('C:/Users/gonza/Downloads/cola-main.zip', repos = NULL, type = 'win.binary')
# devtools::install_local('C:/Users/gonza/Downloads/cola-main.zip')


shinyParallel::installShinyParallel('/home/shiny/colashiny/connecting-landscapes', max.sessions = 25)
install.packages('C:/temp/cola-main.zip')
cola::setup_cola()

library(reticulatete)
reticulate::miniconda_path()
miniconda_uninstall()

reticulate::conda_remove('cola')



if(!require(devtools)){
  install.packages('devtools')
}

library(devtools)
devtools::install_github('connectingLandscapes/cola', dependencies=NA, upgrade = 'never')
devtools::install('N:/My Drive/git/cola')
install.packages('N:/My Drive/git/cola', repos = NULL, type = 'source')

# install.packages('C:/Users/gonza/Downloads/cola-main.zip', repos = NULL, type = 'win.binary')
# 'C:/Users/gonza/Downloads/cola-main.zip'
# reticulate::conda_remove('cola')

library(cola)
# Option A: 1 instruction
setup_cola(force = TRUE,# force any problems
           ask = FALSE, # Allow install conda
           dss = TRUE) # Install DSS libs

# Option B. Getting ask and split steps
setup_cola() # Wait for 2 questions
setup_cola_dss()
cola_dss()


## Solve the user name space problem

library(reticulate)
install_miniconda(path = 'C:/temp/R', update = TRUE, force = FALSE)
Sys.setenv(RETICULATE_MINICONDA_PATH = 'C:/temp/R')


## Edit your parameters
file.edit(file.path(Sys.getenv('HOME'), '.Renviron'))

cola::diagnose_cola()


# Launch cola
library(cola)
cola::cola_dss()



inshp <- 'C:/cola/colaCMT2025080721174905//out_simpts_WRE2025080721182905.shp'
intif <- 'C:/cola/colaCMT2025080721174905//in_points_RSH2025080721182905.tif'
outtif <- 'C:/cola/colaCMT2025080721174905//2.tif'
maxdist <- 10000
shape <- 'linear'
transform <- 'no'
volume <- 1
ncores = as.numeric(Sys.getenv('COLA_NCORES'))
crs = 'None'
py = Sys.getenv('COLA_PYTHON_PATH')
pyscript = system.file(package = 'cola', 'python/crk.py')
cml = TRUE
show.result = TRUE

isProjected(intif)


###
library(cola)


zarr <- lccZarr_py(
  inshp = inshp, intif = intif,
  outtif = outtif,
  tempFolder  = tempFolder,
  maxdist = maxdist,
  smooth = smooth, tolerance = tolerance,
  ncores = ncores, maxram = maxram,
  crs = crs, sci = sci, eci = eci
)

system.file(package = 'cola', 'sampledata/points_sabah_50.shp')

/c/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe
import os
os.chdir('N:/My Drive/git/cola/inst/python')

inshp = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/points_sabah_50.shp'
rg = intif = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/sampleSR.tif'
ofile = outtif = 'C:/cola/out_lfcc2.tif'
tempFolder  = 'C:/cola',
maxdist = 1000000
gRad = smooth = 0
corrTolerance = tolerance = 0
nThreads = ncores = 8
maxram = 6
upCRS = crs = 'None'
sci = 'None'
eci = 'None'
dazarr = 'C:/cola/filea7a845ad590b_dazarr.zarr'
reOrderFile = 'C:/cola/filea7a845ad590b_reOrderFile.csv'
nodeidsFile = 'C:/cola/filea7a845ad590b_nodeidsFile.csv'

## CRK
# # Path to file holding xy coordinates
# xyf = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/points_sabah_50.shp'
# # rg = 'C:\\Users\\gonza\\AppData\\Local\\Temp\\RtmpSYYvAb/resistance.tif'
# rg = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/sampleSR.tif'
# ofile = 'C:\\Users\\gonza\\AppData\\Local\\Temp\\RtmpSYYvAb/kernels.tif'
# dThreshold = '100000'
# tForm = 'linear'
# tkv = 'no'
# kvol = '1'
# nThreads = '8'
# upCRS = 'None'
"


## Step7. Errors ----------------------------------------------

# https://github.com/rstudio/reticulate/issues/838


## Error Warning 3: Cannot find header.dxf (GDAL_DATA is not defined)
# https://stackoverflow.com/questions/45883445/how-to-fix-the-enviroment-variable-gdal-data-path-set

## >   tryB <- tryCatch(system( tryBcmd ), error = function (e) e)
# Traceback (most recent call last):
#   File "C:/Users/Admin/Documents/R/win-library/4.0/cola/python/welcome.py", line 1, in <module>
#   import osgeo
# File "C:\Users\Admin\AppData\Local\R-MINI~1\envs\cola\lib\site-packages\osgeo\__init__.py", line 21, in <module>
#   _gdal = swig_import_helper()
# File "C:\Users\Admin\AppData\Local\R-MINI~1\envs\cola\lib\site-packages\osgeo\__init__.py", line 17, in swig_import_helper
# _mod = imp.load_module('_gdal', fp, pathname, description)
# File "C:\Users\Admin\AppData\Local\R-MINI~1\envs\cola\lib\imp.py", line 242, in load_module
# return load_dynamic(name, filename, file)
# File "C:\Users\Admin\AppData\Local\R-MINI~1\envs\cola\lib\imp.py", line 342, in load_dynamic
# return _load(spec)
# ImportError: DLL load failed while importing _gdal: No se encontrÃ³ el proceso especificado.


# Error.
# > reticulate::py_run_file( system.file("python/welcome.py", package = envName))
# Error in py_module_import(module, convert = convert) :
#   ModuleNotFoundError: No module named 'rpytools'


### Error. Required for IG desktop. Error before installing some libs
# https://stackoverflow.com/questions/64261546/how-to-solve-error-microsoft-visual-c-14-0-or-greater-is-required-when-inst


## Error -- geopandas takes so long
# [1] " --- Installing geopandas"
# Collecting package metadata (current_repodata.json): ...working... done
# Solving environment: ...working... failed with initial frozen solve. Retrying with flexible solve.
# Solving environment: ...working... failed with repodata from current_repodata.json, will retry with next repodata source.

## Error -- problem with rasterio
# https://gis.stackexchange.com/questions/417733/unable-to-import-python-rasterio-package-even-though-it-is-installed

## Error -- no name of conda under "conda info --envs"
# (base) C:\Users\Admin>conda activate C:\Users\Admin\AppData\Local\r-miniconda\envs\cola # activate unnamed env
# conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs ## add unamed envs
#https://stackoverflow.com/questions/57527131/conda-environment-has-no-name-visible-in-conda-env-list-how-do-i-activate-it-a

# C:\Users\Admin\DOCUME~1\VIRTUA~1\colaR3\Scripts\python.exe' -c 'import io, os, sys, setuptools, tokenize; sys.argv[0] = '"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"'; __file__='"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"';f



# Server ON

# <table align="center" border="0" >
#   <tr>
#   <td align="center" width="60%">
#     <a href="http://34.57.191.163:3838/connecting-landscapes" target="_blank">
#       <img src="https://github.com/connectingLandscapes/cola/blob/main/other/servericon_small_wh.png?raw=true" alt="DON'T FORGET THIS">
#         </a>
#         </td>
#         <td align="center" width="40%"> Wait some seconds<br>while the server<br>fully loads</td>
#           </tr>
#           </table>
#
# <table align="center" border="0" >
#   <tr>
#   <td align="center" width="60%">
#     <a href="34.121.114.48:3838/connecting-landscapes" target="_blank">
#       <img src="https://github.com/connectingLandscapes/cola/blob/main/other/servericon_small_wh.png?raw=true" alt="DON'T FORGET THIS">
#         </a>
#         </td>
#         <td align="center" width="40%"> Wait some seconds<br>while the server<br>fully loads</td>
#           </tr>
#           </table>
#
### Server off
#
# <table align="center" border="0" >
#   <tr>
#   <td align="center" width="60%">
#     <a href="" target="_blank">
#       <img src="" alt="DON'T FORGET THIS">
#         </a>
#         </td>
#         <td align="center" width="40%"> No public server at this time.<br>Contact us for specific server time requests<br> at ig299@nau.edu </td>
#           </tr>
#           </table>
#
