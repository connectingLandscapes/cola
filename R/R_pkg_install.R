# Parameters for installation

cola_params <<- list(
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
# .onLoad <- function(libname = 'reticulate', pkgname = 'cola') {
#   ## ask for miniconda
#   user_permission <- utils::askYesNo("Install miniconda? downloads 50MB and takes time")
#
#   if (isTRUE(user_permission)) {
#     reticulate::install_miniconda()
#   } else {
#     message("You should run `reticulate::install_miniconda()` before using this package")
#   }
# }

## Errors
commonErrors <- function(envName = 'cola',
                         libs2Install = cola::cola_params$libs2Install){

  cat(sep = '',
      ' \n We found some errors. Running `', envName, '::setup_cola()` should help you to configure the package.\n',
      ' Please refer to https://github.com/connectingLandscapes/cola/blob/main/known-issues.md for more details.\nHere the diagnostic: (please wait a moment)')

  if (!require(reticulate)){
    cat("  1. `reticulate` package is not installed. Try `install.packages('reticulate')` \n")

  } else{
    cat("  1. `reticulate` package installed. Evaluating next step \n")
    ## Reticulate available. Next check
    miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL)
    (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))

    if (is.null(miniPath) | is.null(miniBase)){
      cat("  2. `miniconda` software is not installed. Try `reticulate::install_miniconda()`\n")

    } else {
      cat("  2. `miniconda` software installed. Evaluating next step \n")
      ## Miniconda available. Next check
      avEnv <- tryCatch(reticulate::conda_list(), error = function (e) NULL)
      if (is.null(avEnv)){
        cat(sep = '', "  3. `", envName, "` conda environment not installed. Try `reticulate::conda_create('",
            envName, "')`\n")
      } else {
        if (class(avEnv) == 'data.frame'){
          if( ! envName %in% avEnv$name){
            cat(sep = '', "  3. `", envName, "` conda environment not installed. Try in R: \n\t    `reticulate::conda_create('",
                envName, "')`\n\t or `conda_create(", envName, ", f = '", system.file('python/python_conda_config.yml', package = "cola"),"')`")
          } else {

            cat("  3. `cola` conda environment installed. Evaluating next step \n")

            ## cola available. Next check
            avLibs <- reticulate::py_list_packages(envname = envName)

            ## No installed
            noInsLibs <- libs2Install[!libs2Install %in% avLibs$package]

            if(length(noInsLibs) != 0){
              cat(sep = '', "  4. Some python libraries aren't installed. Try:\nR: `reticulate::py_install(envname = '", envName,
                  "', channel = 'conda-forge', packages = c('",
                  paste0(noInsLibs, collapse= "', '"), "'))`\n conda CMD: conda install -n ", envName, " ",
                  paste0(noInsLibs, collapse= " ")
                  )
            } else {

              cat("  4. All `cola` conda environment packages installed. Evaluating next step \n")


              ## All libs installed

              (pyCola <- Sys.getenv('COLA_MINICONDA_PATH'))
              (pathCola <- Sys.getenv('COLA_SCRIPTS_PATH'))

              if( dir.exists(pyCola) & file.exists(pathCola) ){

                cat(sep = '', "   === All dependencies and requirements installed. Look for futher details in the repository documentation ===")

              } else {
                cat(sep = '', "  5. Can't connect to python scripts'. The scripts seems to exists, but are not saved ",
                    "by `cola::setup_cola( )` likely because it wasn't able to run the test we designed.\n",
                    "  `Sys.getenv('COLA_MINICONDA_PATH')` and `Sys.getenv('COLA_SCRIPTS_PATH')` should have the paths to `cola` python and folder path.\n",
                    "This error usually arises when some libraries can't be used by python. ",
                    "We solved this by getting the last version of R\nPlease go to:\n\t",
                    "https://github.com/connectingLandscapes/cola/blob/main/known-issues.md")

              }
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




setup_cola <- function(envName = 'cola', nSteps = 5, force = FALSE, yml = TRUE,
                       libs2Install =  c('gdal', 'h5py', 'numexpr', 'rasterio',
                                         'pytables', 'pandas',  'cython', 'numba' ,
                                         'networkit', 'fiona', 'shapely', 'geopandas',
                                         'kdepy', 'scikit-image', 'kdepy')
                       ){

  #envName = cola_params$envName, nSteps = cola_params$nSteps, force = FALSE, libs2Install =  cola_params$libs2Install


  ## Step 1 --- Install reticulate ----------------------------------------------
  cat(sep = '', '  +Step 1/',nSteps, ': Installing & checking reticulate R package\n')

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
    loadLib <- tryCatch(library(reticulate), error = function(e) NULL)
    if( is.null(loadLib) ) {
      commonErrors()
    }
    cat(sep = '', '    `reticulate` installed already!\n')
  }

  library(reticulate)


  ## Step 2 - Install miniconda ----------------------------------------------

  cat (sep = '', '  +Step 2/',nSteps, ' Installing & checking miniconda\n')

  # (sys <- reticulate::import("sys", convert = TRUE))
  # (instMiniConda <- tryCatch(reticulate::install_miniconda(force = TRUE), error = function (e) e))
  # (updMiniConda <- tryCatch(reticulate::miniconda_update(path = miniconda_path()), error = function (e) e))
  # miniconda_uninstall(miniconda_path())
  (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
  (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))


  ## reticulate::miniconda_path() migth return value even if was uninstalled: "C:/Users/Admin/AppData/Local/r-miniconda"
  if ( is.null(miniBase) ){

    #.onLoad(libname = 'reticulate', pkgname = envName)
    user_permission <- utils::askYesNo("Install miniconda? downloads 80MB and takes some minutes")
    if (isTRUE(user_permission)) {
      reticulate::install_miniconda()

      (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
      (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
      if ( is.null(miniBase) ){
        commonErrors()
        }
    } else {
      message("You should run `reticulate::install_miniconda()` before using this package")
      stop()
    }

  } else {
    cat(sep = '', '    miniconda found at ', miniPath, '!\n')
  }


  ## Check if miniconda exists
  (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
  if ( is.null(miniPath) & dir.exists(miniPath) ){
    commonErrors()
  }



  # (pyConf <- reticulate::py_config())
  # py_discover_config() ##  Python version
  # (pyDiscover <- py_discover_config(use_environment = 'base'))


  # Step3. Install your environment ----------------------------------------------
  cat (sep = '', '  +Step 3/',nSteps, ' Installing & checking conda environment\n')
  ## Check again
  (condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL))

  (ymlFile <- system.file('python/python_conda_config.yml', package = "cola"))

  if (is.null(condaLists)){

    user_permission <- utils::askYesNo(paste0("Install ´", envName, "´ environment"))
    if (isTRUE(user_permission)) {

      #  conda env create -f environment.yml
      if(file.exists(ymlFile)){
        system.time(conda_create(envName))
      } else {
        system.time(conda_create(envName, f = ymlFile))
      }

    } else {
      message("You should run `conda_create('", envName, "')` before using this package")
      stop()
    }
  } else {
    if (class(condaLists) == 'data.frame'){
      if( ! envName %in% condaLists$name ){
        user_permission <- utils::askYesNo(paste0("Install '", envName, "' conda environment? Migth take some minutes"))
        if (isTRUE(user_permission)) {
          if(file.exists(ymlFile) & yml){
            instCondEnv <- tryCatch(conda_create(envName, f = ymlFile), error = function(e) e)
          } else {
            instCondEnv <- tryCatch(conda_create(envName), error = function(e) e)
          }
          if ( any(class(instCondEnv) == 'error') & isTRUE(force)){
            envDir <- file.path(miniconda_path(), 'envs', envName)
            if ( dir.exists(envDir) ){
              unlink( envDir, recursive = TRUE, force = TRUE )
              if(file.exists(ymlFile) & yml){
                instCondEnv <- tryCatch(conda_create(envName, f = ymlFile), error = function(e) e)
              } else {
                instCondEnv <- tryCatch(conda_create(envName), error = function(e) e)
              }
            }
          }
          # + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "create" "--yes" "--name" "cola" "python=3.8" "--quiet" "-c" "conda-forge"
        } else {
          message(paste0('You should run `conda_create("', envName,'")` before using this package'))
          stop()
        }
      } else {
        colaexe <- condaLists$python[condaLists$name %in% envName]
        cat (sep = '', '    `', envName, '` conda environment installed in ', colaexe, '\n')
      }
    }
  }

  ## List conda after instaling
  (condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
  if (is.null(condaLists)){
    commonErrors()
  }

  ## Confirm env name
  (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
  if (is.null(pyCola)){
    message(paste0('You should run `conda_create(",', envName ,'")` before using this package'))
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
    ## List conda after instaling
    commonErrors()
    stop()
  } else {
    cat(sep = '', '    The python version is ', pyColaVersion, '\n')
  }




  ## Step4. Install packages ----------------------------------------------
  cat (sep = '', '  +Step 4/',nSteps, ' Installing & checking conda modules\n')

  # Try 3 times to install all the packages with yml file
  for(i in 1:3){

    ## list packages
    (avLibs <- reticulate::py_list_packages(envname = envName))

    ## Install conda packages
    (lib2inst <- libs2Install[! libs2Install %in% avLibs$package])

    if(length(lib2inst) != 0){
      logPkg <- tryCatch(
        reticulate::py_install(
          envname = envName, # python_version = pyCola,
          channel = "conda-forge", packages = lib2inst),
        error = function (e) e)
    } else {
      break()
    }
  }


  ## Try individually

  ## list packages
  (avLibs <- reticulate::py_list_packages(envname = envName))

  ## Install conda packages
  (lib2inst <- libs2Install[! libs2Install %in% avLibs$package])

  if(length(lib2inst) != 0){
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
    commonErrors()
    stop()
  } else {
    cat (sep = '', '    All required conda modules installed!\n')
  }


  cat (sep = '', '  +Step 5/',nSteps, ' Setting up local variables\n')
  # reticulate::py_available()

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
    (cola_scripts_path <- dirname(welcomepy))

    timeMark <- gsub('[[:punct:]]| ', '',
                     format(as.POSIXct(Sys.time(), tz="CET"),
                            tz="America/Bogota",usetz=TRUE))

    pyScript <- system.file("python/s2res.py", package = "cola")
    if ( any(grep(' ', pyScript)) ){
      pyScript <- paste0('"', pyScript, '"')
    }

    outTest <- paste0(tempfile(), timeMark, '.tif')
    (test_suit2res <- paste0(pyCola, ' ', # Python
                             pyScript, ' ', # script
                             system.file("tif/sampleTif.tif", package = "cola"), ' ', #in [1]
                             outTest, #out
                             ' 0 1 100', # min max scale-max
                             ' 1 -9999 None')) # shape nodata proj
    # cat(test_suit2res)
    intCMD <- tryCatch(system(test_suit2res, intern = TRUE, ignore.stdout = TRUE), error = function(e) e)

    if ( dir.exists(cola_scripts_path) & file.exists(outTest) ){

      ## Using cola python as default
      tryCatch(reticulate::use_python(pyCola), error = function(e) e)

      ## Setting cola python as environmental variable
      #Sys.getenv()
      Sys.setenv("COLA_PYTHON_PATH" = pyCola)
      Sys.setenv("COLA_SCRIPTS_PATH" = cola_scripts_path)

      cat (sep = '', '    === Ready to connect landscapes! ===\n')

    } else {
      Sys.unsetenv("COLA_SCRIPTS_PATH")
      Sys.unsetenv("COLA_PYTHON_PATH")
      cat (sep = '', "    -- Final test didn't run. System vars COLA_SCRIPTS_PATH and COLA_PYTHON_PATH removed.\n\t", intCMD)
      commonErrors()
      stop()
    }
  }  else {
    cat (sep = '', "    Error: Can't load conda modules.\n\t", cmdans)

  }
}

# devtools::install_github('connectingLandscapes/cola') ## option 3: None
# library(cola)
# cola::setup_cola()
# setup_cola(envName = 'cola2')
# # # # remove.packages('cola')
# Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))
# origLibs <- installed.packages()
# save(origLibs, file = 'origLibsBeforeDss.RData')
# cola::setup_cola_dss()
# libsafter1 <- installed.packages(); save(libsafter1, file = 'libsafter1.RData')
# (condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
