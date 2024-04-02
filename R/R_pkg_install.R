# Parameters for installation

cola_params <- list(
  ## Environment name
  envName = 'cola',
  ## Libraries to install
  libs2Install = c('gdal', 'h5py', # 'osgeo',
                   'numexpr', 
                   'rasterio', 'pytables',
                   'pandas',  'cython', 'numba' ,
                   'networkit', 'fiona', 'shapely',
                   'geopandas',
                   'kdepy', # 'KDEpy',
                   'scikit-image'),
  ## Number steps
  nSteps = 5
)
# attach(cola_params)


## Install miniconda
.onLoad <- function(libname = 'reticulate', pkgname = 'cola') {
  ## ask for miniconda
  user_permission <- utils::askYesNo("Install miniconda? downloads 50MB and takes time")
  
  if (isTRUE(user_permission)) {
    reticulate::install_miniconda()
  } else {
    message("You should run `reticulate::install_miniconda()` before using this package")
  }
}

## Errors
commonErrors <- function(envName = 'cola', libs2Install){
  
  cat(sep = '', 
      ' We found some errors. Running `', envName, '::setup_cola()` should help you to configure the package.\n',
      ' Please refer to https://github.com/gonzalezivan90/cola/blob/main/known-issues.md for more details') 
  
  if (!require(reticulate)){
    cat("  1. `reticulate` package is not installed. Try `install.packages('reticulate')` \n")
    
  } else{
    
    ## Reticulate available. Next check
    miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL)
    
    if (is.null(miniPath)){
      cat("  2. `miniconda` software is not installed. Try `reticulate::install_miniconda()`\n")  
    
      } else {
      
      ## Miniconda available. Next check
      avEnv <- tryCatch(reticulate::conda_list(), error = function (e) NULL)
      if (is.null(avEnv)){
        cat(sep = '', "  3. `", envName, "` conda environment not installed. Try `reticulate::conda_create('", 
            envName, "')`\n")
      } else {
        if (class(avEnv) == 'data.frame'){
          if( ! envName %in% avEnv$name){
            cat(sep = '', "  3. `", envName, "` conda environment not installed. Try `reticulate::conda_create('", 
                envName, "')`\n")
          } else {
            
            ## cola available. Next check
            avLibs <- reticulate::py_list_packages(envname = envName)
            
            ## No installed
            noInsLibs <- libs2Install[!libs2Install %in% avLibs$package]
            # noInsLibs <- libs2Install[3:6]
            if(length(noInsLibs) != 0){
              cat(sep = '', "  4. Some python libraries aren't installed. Try `reticulate::py_install(envname = '", envName, 
                  "', channel = 'conda-forge', packages = c('", 
                  paste0(noInsLibs, collapse= "', '"), "'))`") 
            } else {
              cat(sep = '', "  5. More problems (none yet, 4 tests passed)... to be completed '") 
            }
          }
        }
      }
    }
  }
  
  if( FALSE) {
    cat(' We found some errors. Please chek the following steps:\n',
        '  1. `reticulate` package installed. Try `require(reticulate)`. Must be TRUE\n',
        '  2. `miniconda` software installed. Try `miniconda_path()`. Must a valid path\n',
        '    2a. If error occurred, try: `miniconda_update()`\n',
        '    2b. If error occurred, try: `reticulate::miniconda_uninstall()` and then `reticulate::install_miniconda()`\n',
        '  3. `cola` conda environment installed. Try `reticulate::conda_list()` to see installed environments\n')
  }
}

setup_cola <- function(force = FALSE){
  
  ## Step 1 --- Install reticulate ----------------------------------------------
  cat(sep = '', '  Step 1/',nSteps, ': Installing & checking miniconda\n')
  
  if (!require(reticulate)){
    
    user_permission <- utils::askYesNo("Install `reticulate` package?")
    
    if (isTRUE(user_permission)) {
      cat(sep = '', '    Installing `reticulate`\n')
      install.packages('reticulate')
    } else {
      message("You should run `install.packages('reticulate')` before using this package")
      stop()
    }
    
  } else {
    cat(sep = '', '    `reticulate` installed already!\n')
    library(reticulate)
  }
  
  
  ## Step 2 - Install miniconda ----------------------------------------------

  cat (sep = '', '  Step 2/',nSteps, ' Installing & checking miniconda\n')
  # (sys <- reticulate::import("sys", convert = TRUE))
  # (instMiniConda <- tryCatch(reticulate::install_miniconda(force = TRUE), error = function (e) e))
  # (updMiniConda <- tryCatch(reticulate::miniconda_update(path = miniconda_path()), error = function (e) e))
  # miniconda_uninstall(miniconda_path())
  (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
  ## reticulate::miniconda_path() migth return value even if was uninstalled: "C:/Users/Admin/AppData/Local/r-miniconda"
  if ( is.null(miniPath) & dir.exists(miniPath) ){
    .onLoad(libname = 'reticulate', pkgname = envName)
  } else {
    cat(sep = '', '    miniconda found at ', miniPath, '!\n')
  }
  
  
  # (pyConf <- reticulate::py_config())
  # py_discover_config() ##  Python version
  # (pyDiscover <- py_discover_config(use_environment = 'base'))
  
  
  # Step3. Install your environment ----------------------------------------------
  cat (sep = '', '  Step 3/',nSteps, ' Installing & checking conda environment\n')
  ## Check again
  condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL)
  
  if (is.null(condaLists)){
    
    user_permission <- utils::askYesNo(paste0("Install ´", envName, "´ environment"))
    if (isTRUE(user_permission)) {
      system.time(conda_create(envName))
    } else {
      message("You should run `conda_create('", envName, "')` before using this package")
      stop()
    }
    
  } else {
    
    if (class(condaLists) == 'data.frame'){
      if( ! envName %in% condaLists$name){
        user_permission <- utils::askYesNo("Install ´cola´ environment")
        if (isTRUE(user_permission)) {
          instCondEnv <- tryCatch(conda_create(envName), error = function(e) e)
          if ( any(class(instCondEnv) == 'error') & isTRUE(force)){
            envDir <- file.path(miniconda_path(), 'envs', envName)
            if ( dir.exists(envDir) ){
              unlink( envDir, recursive = TRUE, force = TRUE )
              instCondEnv <- tryCatch(conda_create(envName), error = function(e) e)
            }
          }
          # + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "create" "--yes" "--name" "cola" "python=3.8" "--quiet" "-c" "conda-forge"
        } else {
          message('You should run `conda_create("conda")` before using this package')
          stop()
        }
      } else {
        colaexe <- condaLists$python[condaLists$name %in% envName]
        cat (sep = '', '    `', envName, '` conda environment installed in ', colaexe,'\n')
      }
    }
  }
  
  condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL)
  
  ## Confirm env name
  (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
  if (is.null(pyCola)){
    message('You should run `conda_create("conda")` before using this package')
    stop()
  } else {
    if(is.character(pyCola) & file.exists(pyCola)){
      cat (sep = '', '    `', envName, '` conda environment named correctly!\n')
      
    } else {
      ## Error -- no name of conda under "conda info --envs"
      
      (possiblePy <- grep(paste0(envName, '/python'), condaLists$python, value = TRUE))
      if( any(possiblePy) ) {
        
        (minibat <- file.path(miniconda_path(), 'condabin/conda.bat'))
        if(file.exists(minibat)){
          cat (sep = '', '    `', envName ,'` was found but in a different path: ', envPath2Add, '\n')
          
          ## Syntax example from reticulate: 
          # + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "install" "--yes" "--name" "cola" "-c" "conda-forge" "pytables"
          ## Actual conda cmd solution
          # conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs ## add unamed envs
          # system('"C:/Users/ig299/AppData/Local/r-miniconda/condabin/conda.bat" "info" "--envs"')
          cmd2add <- (paste0('"',minibat, '" "conda" "config" "--append" "envs_dirs" "', possiblePy, '"'))
          # cat(cmd2add)
          system(cmd2add)
          
          condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL)
          
          ## Confirm env name
          (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
          if (is.null(pyCola)){
            envPath2Add <- dirname(dirname(possiblePy))
            cat (sep = '', '    The `', envName ,'` was found but in a different path.\n',
               '    Try to run in your miniconda shell:\n         conda config --append envs_dirs ', envPath2Add, ' \n')
          }
        }
      }
    }
  }
  
  (test_pyCola <- paste0( pyCola, ' -V'))
  (pyColaVersion <- tryCatch(system( test_pyCola , intern = TRUE), error = function(e) NULL))
  if( is.null(pyColaVersion)){
    message('We can´t check python version.')
    stop()
  } else {
    cat(sep = '', '    The python version is ', pyColaVersion, '\n')
  }
  ## NULL if error +++++++++++++
  
  # conda_remove(envname = 'cola2')
  
  # reticulate::py_list_packages(envname = 'cola') # dont use
  # reticulate::conda_create vs  py_install vs conda_install
  
  # This step does not work -- creating cola with more arguments
  # reticulate::conda_create( envname = 'cola',forge = TRUE, channel = "conda-forge",
  #  packages = c( 'pandas',  'cython', 'geopandas', 'numba' ) )
    

  ## Step4. Install packages ----------------------------------------------
  cat (sep = '', '  Step 4/',nSteps, ' Installing conda modules\n')
  
  
  ## list packages
  avLibs <- reticulate::py_list_packages(envname = envName)
  for( l in 1:length(libs2Install)){ # l = 1
    (lib2inst <- libs2Install[l])
    if( ! lib2inst %in% avLibs$package ){
      cat(paste0(' --- Installing `',  libs2Install[l], '` module\n'))
      
      logPkg <- tryCatch(
        reticulate::py_install( 
          envname = envName, # python_version = pyCola,
          channel = "conda-forge", packages = lib2inst), 
        error = function (e) e)
      
      # reticulate::py_install(envname = envName, channel = "conda-forge", packages = 'numexpr')
      # -- Installing `gdal` module + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "install" "--yes" "--name" "cola" "-c" "conda-forge" "gdal"
      
      ## If there's a problem installing trought conda-forge, use PIP
      if( any(!is.null(logPkg)) ){
        print(lib2inst)
        logPkg2 <- tryCatch(reticulate::py_install( envname = envName, 
                                                    pip = TRUE,
                                                    packages = lib2inst), 
                            error = function (e) e)
      }
      avLibs <- reticulate::py_list_packages(envname = envName)
    }
  }
  
  avLibs <- reticulate::py_list_packages(envname = envName)
  
  ## Installed
  insLibs <- libs2Install[libs2Install %in% avLibs$package]
  ## No installed
  noInsLibs <- libs2Install[!libs2Install %in% avLibs$package]
  # noInsLibs <- c('a', 'f')
  
  if( length(noInsLibs) > 0 ){
    cat (sep = '', '  `', envName, '` conda environment installed!\n')
    cat(sep = '', "   Some python libraries aren't installed. Try `reticulate::py_install(envname = '", envName, 
        "', channel = 'conda-forge', packages = c('", 
        paste0(noInsLibs, collapse= "', '"), "')`") 
  } else {
    cat (sep = '', '  All required modules installed!\n')
  }
  
  
  cat (sep = '', '  Step 5/',nSteps, ' Setting up local variables\n')
  
  # reticulate::py_available()
  
  ## Using cola python as default
  reticulate::use_python(pyCola)
  
  ## Setting cola python as environmental variable
  #Sys.getenv()
  Sys.setenv("COLA_MINICONDA_PATH" = pyCola)
  
  
  ## Connecting lib paths
  (welcomepy <- system.file("python/welcome.py", package = "cola"))
  #(welcomepy <- file.path('N:/My Drive/git/cola/inst/python/welcome.py'))
  #read.delim(welcomepy)
  if (any(grep(' ', welcomepy))){
    welcomepy <- paste0('"', welcomepy, '"')
  }
  
  #tryA <- tryCatch(reticulate::py_exe(system.file("python/welcome.py", package = "cola")), error = function (e) e)
  (cmd2test <- paste0( pyCola, ' ', welcomepy)); #cat(tryBcmd)
  (cmdans <- tryCatch( system( cmd2test , intern = TRUE), error = function (e) e))
  
  if (any(grep('WELCOME ', cmdans))){
    #libP <- .libPaths()
    #cola_scripts_path <- file.path(libP, 'cola/python')
    cola_scripts_path <- dirname(welcomepy)
    if (dir.exists(cola_scripts_path)){
      Sys.setenv("COLA_SCRIPTS_PATH" = pyCola)
      cat (sep = '', '    Ready to connect landscapes!\n')
      
    }
  }
  ## 
  # reticulate::py_run_file( system.file("python/welcome.py", package = envName))
  # reticulate::use_python(pyCola)
  
  
  # (test_cmd1 <- paste( pyCola, ' -c "import os; print(os.getcwd())')); system( test_cmd2 )
  # (test_cmd2 <- paste0( pyCola, ' -V')); system( test_cmd2 )
  
  
  
  ## Step6. For debugging ----------------------------------------------
  

  if (FALSE){
    
    # tempPy <- paste0(tempfile(), '.py')
    # tempPyFun <- paste0(tempfile(), '.py')
    # download.file('https://raw.githubusercontent.com/gonzalezivan90/connectscape/main/welcome.py', destfile = tempPy)
    # download.file('https://raw.githubusercontent.com/gonzalezivan90/connectscape/main/cola_functions.py', destfile = 'cola_functions.py')
    # readLines(tempPy)
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
    
    # StepX. Get Git package
    # devtools:::install_github("gonzalezivan90/cola")## Select option 3: NONE
    # remove.packages( "cola" )
    # list.files("C:/Users/Admin/Documents/R/win-library/4.0/", recursive = TRUE)
    # C:\Users\Admin\Documents\R\win-library\4.0\cola no py files here
    # .libPaths()
  }
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
  
}

# C:\Users\Admin\DOCUME~1\VIRTUA~1\colaR3\Scripts\python.exe' -c 'import io, os, sys, setuptools, tokenize; sys.argv[0] = '"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"'; __file__='"'"'C:\\Users\\Admin\\AppData\\Local\\Temp\\pip-install-q2renvli\\networkit_bbb1ed5652414ced8de6bd7c807b6c54\\setup.py'"'"';f 
# devtools::install_github('gonzalezivan90/cola')
# library(cola)
# setup_cola(force = FALSE)