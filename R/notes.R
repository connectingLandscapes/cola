############################# installing #################################
# .rs.restartR()
# remove.packages('cola')
# .rs.restartR()
# devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') ## option 3: None
# library(cola)
# cola::diagnose_cola()
# cola::setup_cola()
# Sys.getenv(c('COLA_PYTHON_PATH', 'COLA_SCRIPTS_PATH'))
# Sys.getenv(c('COLA_DATA_PATH', 'COLA_SCRIPTS_PATH', 'COLA_DSS_UPL_MB', 'COLA_VIZ_THRES_PIX', 'COLA_VIZ_RES_NCOL', 'COLA_VIZ_RES_NROW', 'COLA_NCORES' ))
# cola::setup_cola_dss()
# .rs.restartR()
# library(cola)
# cola::cola_dss()

# origLibs <- installed.packages()
# cola::setup_cola_dss()




############################# Shiny end #################################
# """""""" ========== -------------------
#
# sudo cp /home/shiny/connecting-landscapes/R/cola_tools.R /srv/shiny-server/cola/cola_tools.R
# sudo cp /home/shiny/connecting-landscapes/R /srv/shiny-server/cola -R
# sudo cp /home/shiny/connecting-landscapes/R/*.R /srv/shiny-server/cola -R
# sudo cp /home/shiny/connecting-landscapes/R/* /home/shiny/cola/connecting-landscapes/.; sudo rm /srv/shiny-server/connecting-landscapes -R
# cp /home/shiny/cola/inst/app/* -R /srv/shiny-server/cola/.
# http://18.190.126.82:3838/cola

## orig
# cp /home/shiny/connecting-landscapes/R/* /home/shiny/cola/connecting-landscapes/. -R; sudo rm /srv/shiny-server/connecting-landscapes -R
# shinyParallel::installShinyParallel('/home/shiny/cola/connecting-landscapes/', max.sessions = 20, users.per.session = 10)

# sudo rm /srv/shiny-server/connecting-landscapes -R; sudo cp /home/shiny/cola/inst/app/* /home/shiny/colashiny/connecting-landscapes/. -R
# shinyParallel::installShinyParallel('/home/shiny/colashiny/connecting-landscapes', max.sessions = 20, users.per.session = 10)
# http://18.190.126.82:3838/connecting-landscapes
# http://18.190.126.82:3838/connecting-landscapes/?admin
# http://18.190.126.82:3838/connecting-landscapes/?admin1

# sudo su shiny; cd /home/shiny/cola; git add . ; git commit -m " "; git push
# sudo su shiny; cd /home/shiny/connecting-landscapes; git add . ; git commit -m "Change coordiantes()"; git push
# git pull main --rebase --autostash
# sudo chown -R shiny:shiny .
# git stash
# remove before commit, split or lost it
# git pull connectscape |||  git pull --rebase --autostash || git pull origin HEAD


# https://github.com/settings/tokens/1354187156/regenerate

# git pull connectscape
# cd /home/shiny/connecting-landscapes/; git pull .

# R -e "shinyParallel::installShinyParallel('/home/shiny/cola/connecting-landscapes', max.sessions = 25)"
# ##sudo su - -c "R -e \"shinyParallel::installShinyParallel('/home/shiny/cola/connecting-landscapes', max.sessions = 25)\"" #
# sudo chown -R shiny:shiny .

# sudo cat /var/log/shiny-server/cola
# sudo rm /home/shiny/tmpR/leafSim.RDatasudo cp /home/vmuser/gedivis /srv/shiny-server/gedivis -R### COLA web app.


############################# Shiny end #################################


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

##
# reticulate::py_run_file( system.file("python/welcome.py", package = envName))
## Step6. For debugging ----------------------------------------------



if (FALSE){

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

}



# setwd(Sys.getenv('R_USER'))
# options(ZZZ = 'test')

'C:/Program Files/R/R-4.0.2/etc/Rprofile.site'
'C:/Program Files/R/R-4.0.2/library/base/R/Rprofile'
'C:/Users/Admin/Documents/R/win-library/4.0/'
'C:/Users/Admin/Documents/R/win-library/4.0/packrat/resources/init-rprofile.R' #init.R
'C:/Users/Admin/Documents/R/win-library/4.0/usethis/html/ini'

# Rprofile.site


