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