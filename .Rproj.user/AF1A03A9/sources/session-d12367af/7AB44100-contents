
## Step 0 --- Install reticulate
if (!require(reticulate)){
  install.packages('reticulate')
} else {
  library(reticulate)
}


## Step 1 --- Install miniconda
(sys <- reticulate::import("sys", convert = TRUE))
(instMiniConda <- tryCatch(reticulate::install_miniconda(), error = function (e) e))
(updMiniConda <- tryCatch(reticulate::miniconda_update(path = miniconda_path()), 
                          error = function (e) e))

# reticulate::install_miniconda(force = TRUE)
(pyConf <- reticulate::py_config())
(condaLists <- reticulate::conda_list())
py_discover_config() ##  Python version
(pyDiscover <- py_discover_config(use_environment = 'base'))

envname <- 'cola'


# Step2. Install your environment
if (! envname %in% condaLists$name){
  system.time(conda_create("cola"))
  
  # reticulate::py_list_packages(envname = 'cola') # dont use
  # reticulate::conda_create vs  py_install vs conda_install
  
  # This step does not work -- creating cola with more arguments
  #reticulate::conda_create( envname = 'cola',forge = TRUE,   channel = "conda-forge",
  #  packages = c( 'pandas',  'cython', 'geopandas', 'numba' ) )
  
}

(condaLists <- reticulate::conda_list())
# conda_remove(envname = 'cola')

## Python version
(pyCola <- subset(condaLists, name == 'cola')$python)

## list packages
reticulate::py_list_packages(envname = 'cola')
{
#            package      version                requirement     channel
# 1            bzip2        1.0.8                bzip2=1.0.8 conda-forge
# 2  ca-certificates     2024.2.2   ca-certificates=2024.2.2 conda-forge
# 3           libffi        3.4.2               libffi=3.4.2 conda-forge
# 4        libsqlite       3.45.2           libsqlite=3.45.2 conda-forge
# 5          libzlib       1.2.13             libzlib=1.2.13 conda-forge
# 6          openssl        3.2.1              openssl=3.2.1 conda-forge
# 7              pip         24.0                   pip=24.0 conda-forge
# 8           python       3.8.19              python=3.8.19 conda-forge
# 9       setuptools       69.2.0          setuptools=69.2.0 conda-forge
# 10              tk       8.6.13                  tk=8.6.13 conda-forge
# 11            ucrt 10.0.22621.0          ucrt=10.0.22621.0 conda-forge
# 12              vc         14.3                    vc=14.3 conda-forge
# 13    vc14_runtime  14.38.33130   vc14_runtime=14.38.33130 conda-forge
# 14  vs2015_runtime  14.38.33130 vs2015_runtime=14.38.33130 conda-forge
# 15           wheel       0.43.0               wheel=0.43.0 conda-forge
# 16              xz        5.2.6                   xz=5.2.6 conda-forge
}


### Required for IG desktop:
# https://stackoverflow.com/questions/64261546/how-to-solve-error-microsoft-visual-c-14-0-or-greater-is-required-when-inst

# Step3. Install packages
## Install more packages
libs2Install <- c(  'gdal', 'h5py',
                    # 'osgeo', 
                    'rasterio', 'pytables',
                    'pandas',  'cython', 'numba' , 
                    'networkit', 'fiona', 'shapely',
                    # 'geopandas-base==0.14.0 ', 
                    #'geopandas-base',
                    #'geopandas==0.14.0',
                    'geopandas',
                    'KDEpy', # problem
                    'scikit-image')

avLibs <- reticulate::py_list_packages(envname = envname)
for( l in 1:length(libs2Install)){ # l = 8
  (lib2inst <- libs2Install[l])
  if( ! lib2inst %in% avLibs$package ){
    print(paste0(' --- Installing ',  libs2Install[l]))
    logPkg <- tryCatch(reticulate::py_install( envname = envname, 
                                     # python_version = pyCola,
                                     channel = "conda-forge",
                                     packages = lib2inst), 
             error = function (e) e)
    avLibs <- reticulate::py_list_packages(envname = envname)
  }
}

avLibs <- reticulate::py_list_packages(envname = envname)
## Installed
libs2Install[libs2Install %in% avLibs$package]
## No installed
libs2Install[!libs2Install %in% avLibs$package]

# https://github.com/rstudio/reticulate/issues/838


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


# Step4. Get Git package
# curl::curl_download()
# download.packages()
tempPy <- paste0(tempfile(), '.py')
tempPyFun <- paste0(tempfile(), '.py')
download.file('https://raw.githubusercontent.com/gonzalezivan90/connectscape/main/welcome.py', destfile = tempPy)
download.file('https://raw.githubusercontent.com/gonzalezivan90/connectscape/main/cola_functions.py', destfile = 'cola_functions.py')
readLines(tempPy)

(test_cmd <- paste( pyCola, tempPy))
system( test_cmd )


(test_cmd2 <- paste( pyCola, ' -c "import ', 
                     paste0(
                       gsub(replacement = 'skimage', pattern = 'scikit-image', libs2Install), 
                       collapse = ', '), '"'))

(test_cmd2 <- paste( pyCola, ' -c "import os; print(os.getcwd())'))
system( test_cmd2 )
(test_cmd3 <- paste( pyCola, ' -c "import cola_functions as cf; print(1)'))
system( test_cmd3 )

if(!require("devtools")){
  install.packages("devtools")
}

# devtools:::install_github("gearslaboratory/gdalUtils")
#devtools:::install_github("gonzalezivan90/connectscape")

# C:\Users\Admin\DOCUME~1\VIRTUA~1\colaR3\Scripts\python.exe' -c 'import io, os, sys, setuptools, tokenize; sys.argv[0] = '"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"'; __file__='"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"';f 