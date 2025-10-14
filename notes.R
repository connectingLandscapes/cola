"

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


## Example of back end
python script.py input output parameter1 parameter2 …



C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe 
C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/python/s2res.py
C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/sampleTif.tif 
C:/cola/colaMVX2025070801495305//out_surface_FYG2025070801542005.tif
0.068 0.999 150 1 -9999 None 


# CML for corridors
C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe
C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/python/lcc.py

C:/cola/colaMVX2025070801495305//out_simpts_ENP2025070801540405.shp 
C:/cola/colaMVX2025070801495305//out_surface_FYG2025070801542005.tif
C:/cola/colaMVX2025070801495305//out_lcc_NJX2025070801575205.tif
1000000 0 5 4 None 2>&1

## Run lcc
corridors <- lcc_py(inshp = 'C:/cola/colaMVX2025070801495305//out_simpts_ENP2025070801540405.shp',
                    intif = 'C:/cola/colaMVX2025070801495305//out_surface_FYG2025070801542005.tif', 
                    outtif = ,'C:/cola/colaMVX2025070801495305//out_lcc_fromR.tif',
                    maxdist = 1000000, smooth = 20, tolerance = 100,
                    ncores = 4, cml = TRUE, show.result = TRUE)
 
## Run crk
library(cola)
kernels <- crkJoblib_py(
  inshp = 'C:/cola/singye/sampled_points_tiger_211_proj.shp',
                    intif = 'C:/cola/singye/crk_tiger_impact_c10_proj.tif', 
                    outtif = ,'C:/cola/singye/out_crk.tif',
                    maxdist = 440837.4, shape = 'linear', transform = 'yes',
                    volume = 100000000, ncores = 4, cml = TRUE, show.result = TRUE)

#################################3
inshp = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/points_sabah_50.shp'
intif = 'C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sampledata/sampleSR.tif'


inshp = 'C:/cola/colaSWV2025100802211305//out_simpts_EHB2025100802213205.shp'
intif = 'C:/cola/colaSWV2025100802211305//out_surface_VKG2025100802213505.tif'

inshp = 'C:/cola/colaSWV2025100802211305/out_simpts_EHB2025100802213205.shp' 
intif = 'C:/cola/colaSWV2025100802211305/out_surface_VKG2025100802213505.tif'


out_crk <<- tryCatch(
  crk_py(
    inshp = inshp, intif = intif,
    outtif = 'C:/cola//crk.tif',
    maxdist = 1000,
    transform = 'linear',
    shape = 'no',
    volume = '1'),
  error = function(e) list(log = e$message, file = ''))

out_cola <<- tryCatch(
  cola::crk_py(
    inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp'), 
    intif = system.file(package = 'cola', 'sampledata/sampleSR.tif'),
    outtif = tempfile(fileext = '.tif'),
    maxdist = 1000,
    transform = 'linear',
    shape = 'no',
    volume = '1'),
  error = function(e) list(log = e$message, file = ''))

C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/python/crk.py C:/cola/colaSWV2025100802211305//out_simpts_EHB2025100802213205.shp C:/cola/colaSWV2025100802211305//out_surface_VKG2025100802213505.tif C:/cola/colaSWV2025100802211305//out_crk_IUV2025100802221305.tif 100000 linear no 1 4 None 2>&1 

C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/python/crk.py C:/cola/colaSWV2025100802211305/out_simpts_EHB2025100802213205.shp C:/cola/colaSWV2025100802211305/out_surface_VKG2025100802213505.tif C:/cola/colaSWV2025100802211305//out_crk_IUV.tif 10000 no linear 1 4 None 2>&1 
C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/python/crk.py C:/cola/colaCXZ2025100801562805/out_simpts_KAS2025100801574405.shp C:/cola/colaCXZ2025100801562805/out_surface_PRY2025100801574705.tif C:/cola/colaSWV2025100802211305//crk.tif 1000 no linear 1 4 None 2>&1 


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