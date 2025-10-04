#' @title  Parameters for installation
#' @description Defines the default cola installation parameters
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
# cola_params <<- list(
#   ## Environment name
#   envName = 'cola',
#
#   ## Libraries to install in the python conda environment
#   libs2Install = c('geopandas', # this first to avoid problems
#                    'gdal', 'h5py',
#                    'numexpr', 'zarr',
#                    'rasterio', 'pytables',
#                    'pandas',  'cython', 'numba' ,
#                    'networkit', 'fiona', 'shapely',
#                    'kdepy', 'joblib',
#                    'scikit-image'),
#   yml = TRUE,
#   ## Number steps
#   nSteps = 5,
#
#   ## R packages for the DSS
#   libs2colaDSS = c(
#     'markdown', 'rmarkdown',
#     'knitr', 'units',
#     "reshape2", 'bit', 'digest', 'dplyr',
#     'tidyverse', 'DT', 'ggplot2',
#     # debug install order: htmltools >> shiny >> shinyWidgets
#     'htmlwidgets', 'htmltools', ## Before shiny
#     'magrittr', 'RColorBrewer', 'viridis',
#     ## Spatial
#     # "rgeos", "rgdal", 'raster', ## old
#     'sf', 'terra',
#     'rlang', "leaflet", "leaflet.extras",
#     'gdalUtilities',
#     ## Shiny
#     'shiny',  ## Before shiny plugins
#     "shinydashboard",  "shinycssloaders",
#     'shinydashboardPlus', 'shinyjs',
#     'shinyWidgets', 'dashboardthemes',
#     "highcharter", 'plotly')
# )
# # attach(cola_params)


#' @title Diagnose  \emph{COLA} installation
#' @description Runs some test assessing the correct cola installation
#' @param envName The environment name. 'cola' by default
#' @param libs2Install The name of python libraries
#' @return NULL. Prints in console five (5) logs regarding different steps
#' @examples
#' diagnose_cola()
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#'
#' @export

diagnose_cola <- function(envName = 'cola',
                          libs2Install = c('gdal', 'h5py', 'zarr',
                                           'numexpr',
                                           'rasterio', 'pytables',
                                           'pandas',  'cython', 'numba' ,
                                           'networkit', 'fiona', 'shapely',
                                           'geopandas',
                                           'kdepy',
                                           'scikit-image')){

  cat(sep = '',
      ' \n\n\tDiagnosing CoLa\nWe will look for errors. Running `cola::setup_cola()` should help you to configure the package.\n',
      ' Please refer to https://github.com/connectingLandscapes/cola/blob/main/known-issues.md for more details.\nHere the diagnostic: (please wait a moment)\n')

  if (!require(reticulate)){
    cat("  1. `reticulate` package is not installed. Try `install.packages('reticulate')` \n")

  } else{
    cat("  1. `reticulate` package installed. Evaluating next step \n")
    library(reticulate)
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
        if ( class(avEnv) == 'data.frame' ){
          if( !envName %in% avEnv$name ){
            cat(sep = '', "  3. `", envName, "` conda environment not installed. Try in R: \n\t    `reticulate::conda_create('",
                envName, "')`\n\t or in conda CMD `conda env create -f ", system.file('python/python_conda_config.yml', package = "cola"),"`\n")
          } else {

            cat("  3. `cola` conda environment installed. Evaluating next step \n")

            ## cola available. Next check
            avLibs <- reticulate::py_list_packages(envname = envName)

            ## No installed
            noInsLibs <- libs2Install[!libs2Install %in% avLibs$package]

            (noInsLibs <- libs2Install[
              !((libs2Install %in% avLibs$package) |
                  (gsub('==', '=', libs2Install) %in% avLibs$requirement))])

            if(length(noInsLibs) != 0){
              cat(sep = '', "  4. Some python libraries aren't installed. Try:\nR: `reticulate::py_install(envname = '", envName,
                  "', channel = 'conda-forge', packages = c('",
                  paste0(noInsLibs, collapse= "', '"), "'))`\n conda CMD: conda install -n ", envName, " ",
                  paste0(noInsLibs, collapse= " "), '\n'
              )
            } else {

              cat("  4. All `cola` conda environment packages installed. Evaluating next step \n")


              ## All libs installed

              (pyCola <- Sys.getenv('COLA_PYTHON_PATH'))
              (pathCola <- Sys.getenv('COLA_SCRIPTS_PATH'))

              (pyCola2check <- gsub(' .+', '', pyCola))

              if( file.exists(pyCola2check) & dir.exists(pathCola) ){

                cat(sep = '', "   === All dependencies and requirements installed === \nLook for futher details in the repository documentation\n\n")

              } else {
                cat(sep = '', "  5. Can't connect to python scripts'. The scripts seems to exists, but are not saved ",
                    "by `cola::setup_cola( )` likely because it wasn't able to run the test we designed.\n",
                    "  `Sys.getenv('COLA_PYTHON_PATH')` and `Sys.getenv('COLA_SCRIPTS_PATH')` should have the paths to `cola` python and folder path.\n",
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

  warning(paste0('We are using a new Python versiom. CoLs requires 3.12.11\n',
                 'To update your CoLa version, please run the following commands and try setup_cola() again:\n',
                 '\n\treticulate::conda_remove(envname = "cola") \n',
                 '\tfile.remove("', colaDir, '")\n',
                 '\treticulate::conda_update()\n',
                 '\tsetup_cola(pyVer = "3.12.11")\n'
  ))
  if( FALSE) {
    cat(' We found some errors. Please chek the following steps:\n',
        '  1. `reticulate` package installed. Try `require(reticulate)`. Must be TRUE\n',
        '  2. `miniconda` software installed. Try `miniconda_path()`. Must a valid path\n',
        '    2a. If error occurred, try: `miniconda_update()`\n',
        '    2b. If error occurred, try: `reticulate::miniconda_uninstall()` and then `reticulate::install_miniconda()`\n',
        '  3. `cola` conda environment installed. Try `reticulate::conda_list()` to see installed environments\n')
  }
}


#' @title Install \emph{cola} conda environment
#' @description Installs `cola` conda environment
#' @param envName The environment name. 'cola' by default
#' @param ymlFile Path to YML to use. Default is NULL
#' @return NULL. Prints in console logs regarding different steps
#' @examples
#' install_cond_env(envName = 'cola',
#'     ymlFile = system.file('python/python_conda_config.yml', package = "cola"))
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
install_cond_env <- function(envName, useYML = TRUE, ymlFile = NULL,
                             packages = NULL , pv = '3.12.11'){


  if(useYML & is.null(ymlFile) ){
    ymlFile = system.file('python/python_conda_config.yml', package = "cola")
  }
  newYmlFile <- ymlFile
  insCondLog <- ''

  if( useYML ){
    ## Create env using yml file
    (instCondEnv <- paste0(conda_binary(), ' "env" "create" "--file" "', newYmlFile, '"'))
    cat('   Creating conda using YML file:', instCondEnv, '\n')
    insCondLog <- tryCatch(system(instCondEnv, intern = TRUE), error = function(e) e$message) #

    if( any(grep('Could not solve for environment specs', insCondLog)) ){
      cat('   YML creation failed. Trying conda_create("', envName, '")\n')
      insCondLog <- tryCatch(conda_create(envname = envName, packages = packages, python_version = pv
                                          ), error = function(e) e$message)
    }

    if( any(grep(' prefix already exists', insCondLog)) ){
      cat( 'ERROR: ', insCondLog, '\n',
           'Try conda_remove(envname ="', envName, '"); conda_create("', envName, '")\n')
    }
  } else {
    ## Creating env with no yml
    insCondLog <- tryCatch(conda_create(envname =  envName,  packages = packages, python_version = pv), error = function(e) e$message)
    if( any(grep(' prefix already exists', insCondLog)) ){
      cat( 'ERROR: ', insCondLog, '\n',
           'Try conda_remove(envname ="', envName, '"); conda_create("', envName, '")\n')
    }
  }
  return(insCondLog)
}

#' @title Install  \emph{COLA}
#' @description Installs cola software. Includes `reticulate` R package, miniconda, cola conda environment and libraries,
#' and defines the ways R and python should interact
#' @param envName The environment name. 'cola' by default
#' @param libs2Install The name of python libraries
#' @param nSteps The number of steps for printing log in console
#' @param pyVer Python version to use. Default is '3.12.11'
#' @param force Force miniconda installation? Passed to `reticulate::install_miniconda()`
#' @param yml Use YML file to build the conda environment? Default TRUE
#' @param ask Ask before installing reticulate, minoconda, and cola conda environment?. Default TRUE. If FALSE, installations are made without asking
#' @param dss Install cola DSS as well? Default FALSE
#' @param onlyIndividual Try installing libraries one by one? Default FALSE
#' @return NULL. Prints in console logs regarding different steps
#' @examples
#' setup_cola()
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
setup_cola <- function( envName = 'cola', nSteps = 5, force = FALSE,
                        yml = FALSE, onlyIndividual = FALSE, ask = TRUE,
                        dss = FALSE, pyVer = '3.12.11',
                        libs2Install =
                          c('zarr', 'networkit',
                            # 'zarr==3.1.3', 'networkit==11.0',
                            'geopandas', 'gdal', 'h5py', 'numexpr', 'rasterio',
                            'pytables', 'pandas',  'cython', 'numba' ,
                            'fiona', 'shapely',
                            'kdepy', 'scikit-image', 'kdepy')
){

  # envName = 'cola'; nSteps = 5; force = FALSE; yml = FALSE; onlyIndividual = F; ask = FALSE; dss = TRUE
  if ( !ask ){
    user_permission <- TRUE
  }

  ## Step 1. Install reticulate ----------------------------------------------
  cat(sep = '', '  +Step 1/',nSteps, ': Installing & checking reticulate R package\n')

  if (!require(reticulate, quietly = TRUE)){

    if (ask){
      user_permission <- utils::askYesNo("Install `reticulate` package?")
    }

    if (isTRUE(user_permission) | !ask) {
      cat(sep = '', '    Installing `reticulate`\n')
      install.packages('reticulate')
    } else {
      message("You should run `install.packages('reticulate')` before using this package")
      stop()
    }

  } else {
    loadLib <- tryCatch(library(reticulate, quietly = TRUE), error = function(e) NULL)
    if( is.null(loadLib) ) {
      diagnose_cola()
    }
    cat(sep = '', '    `reticulate` package installed already!\n')
  }

  library(reticulate)


  ## Step 2. Install miniconda ----------------------------------------------

  cat (sep = '', '  +Step 2/', nSteps, ' Installing & checking miniconda\n')

  # (sys <- reticulate::import("sys", convert = TRUE))
  # (instMiniConda <- tryCatch(reticulate::install_miniconda(force = TRUE), error = function (e) e))
  # (updMiniConda <- tryCatch(reticulate::miniconda_update(path = miniconda_path()), error = function (e) e))
  # miniconda_uninstall(miniconda_path())
  (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
  (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))


  # Fix miniconda / conda bat not found.
  # https://stackoverflow.com/questions/60974507/error-unable-to-find-conda-binary-is-anaconda-installed
  if ( dir.exists(miniPath) & is.null(miniBase) ){
    # "C:/Users/Admin/AppData/Local/r-miniconda/" "C:/Users/Admin/miniconda3/condabin/conda.bat"
    (newCondaPath <- gsub('//', '/', gsub('\\', '/', fixed = TRUE, file.path(Sys.getenv('HOMEDRIVE'), Sys.getenv('HOMEPATH'), 'miniconda3/condabin/conda.bat'))))

    if( file.exists(newCondaPath) ){
      options(reticulate.conda_binary = newCondaPath)
      # options(reticulate.conda_binary = 'C:/Users/Admin/miniconda3/condabin/conda.bat')## Works
      (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
    }

    (newCondaPath2 <- file.path(miniPath, 'condabin/conda.bat'))
    if( file.exists(newCondaPath) & is.null(miniBase) ){
      options(reticulate.conda_binary = newCondaPath)
      (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
    }

    if( is.null(miniBase) ){
      message(
        paste0(
          "   Miniconda is likely to be installed but can't be found. Please try found 'conda' or 'conda.bat' files at:\n",
          '\t', newCondaPath2, '\n', '\t', newCondaPath , '\n',
          '   Use R the command `options(reticulate.conda_binary = "C:/path/to/conda.bat")`\n')
      )
      stop()
    }

  }


  ## reticulate::miniconda_path() migth return value even if was uninstalled: "C:/Users/Admin/AppData/Local/r-miniconda"
  if ( is.null(miniBase) ){

    #.onLoad(libname = 'reticulate', pkgname = envName)
    if(ask){
      user_permission <- utils::askYesNo("Install miniconda? downloads 80MB and takes some minutes")
    }

    if (isTRUE(user_permission) | !ask) {
      reticulate::install_miniconda()

      (miniPath <- tryCatch(reticulate::miniconda_path(), error = function (e) NULL))
      (miniBase <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
      if ( is.null(miniBase) ){
        diagnose_cola()
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
    diagnose_cola()
  }

  # (pyConf <- reticulate::py_config())
  # py_discover_config() ##  Python version
  # (pyDiscover <- py_discover_config(use_environment = 'base'))

  # Step 3. Install your environment ----------------------------------------------
  cat (sep = '', '  +Step 3/', nSteps, ' Installing & checking conda environment\n')

  ## List conda after instaling
  (condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL))
  if (is.null(condaLists)){
    message(paste0('Is miniconda installed? Please run `reticulate::conda_list()` and be sure it returns a list of conda environments'))
    diagnose_cola()
    stop()
  }

  ## Check version
  (pyBase <- tryCatch( subset(condaLists, name == 'base')$python, error = function (e) NULL) )
  (pyBaseVersion <- tryCatch(system( paste0( pyBase, ' -V') , intern = TRUE), error = function(e) NULL))

  if( is.null(pyBaseVersion)){
    message('We can´t check python version.')
    ## List conda after instaling
    diagnose_cola()
    stop()
  } else {
    cat( '    Using', pyBaseVersion, 'for ´base´ conda environment\n')
    (numPyVers <- as.numeric(substr(0, 4, x = gsub('Python|python| ', '', pyBaseVersion))))
    (numPyVers2 <- gsub('Python|python| ', '', pyBaseVersion))
    if ( numPyVers <= 3.11 ){
      (colaDir <- dirname(subset(condaLists, name == envName)$python))
      stop(paste0('The current ', pyBaseVersion, ' version need to be updated.\n',
                  'For this, please run the following commands and try setup_cola() again:\n',
                  '\n\treticulate::conda_remove(envname = "cola") \n',
                  '\tfile.remove("', colaDir, '")\n',
                  '\treticulate::conda_update()\n',
                  '\tsetup_cola()\n'
      ))
      # reticulate::conda_remove(envname = "cola")
      #(tryCatch( file.remove( colaDir ), error = function (e) NULL) )
    }
  }

  (numPyVers3  <- ifelse(!is.null(pyVer), pyVer, numPyVers2))

  ## Conda exists // check for cola
  if ( class(condaLists) == 'data.frame' ){

    ## Error: Env exists but empty -- prob broken uninstall
    if( envName %in% condaLists$name ){
      (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
      if (!file.exists(pyCola)){
        message(paste0(' Uninstalling corrupt previous installation'))
        tryCatch(conda_remove(envName))
        (condaLists <- tryCatch( reticulate::conda_list(), error = function (e) NULL))
      }
    }

    ## Env not existing
    if( !envName %in% condaLists$name ){

      # (ymlFile <- 'N:/Mi unidad/git/cola/inst/python/python_conda_config.yml'); read.delim(ymlFile)
      (ymlFile <- system.file('python/python_conda_config.yml', package = "cola"))
      ymlTxt <- readLines(ymlFile)
      ymlTxt <- gsub('name: .+', paste0('name: ', envName), ymlTxt) # Change env name
      ymlTxt <- gsub('^prefix: .+', paste0('prefix: ', file.path(miniPath, 'envs', envName)), ymlTxt) # Change env path
      (newYmlFile <- paste0(tempfile(), 'newYml.yml'))  # read.delim(newYmlFile)
      writeLines(text = ymlTxt, con = newYmlFile)

      if (ask){
        user_permission <- utils::askYesNo(paste0("Install '", envName, "' conda environment? Migth take some minutes"))
      }


      if ( isTRUE(user_permission) ) {

        if(!onlyIndividual){
          insCondLog <- install_cond_env(
            envName = envName, useYML = yml,
            packages = libs2Install,
            ymlFile = newYmlFile, pv = numPyVers3)

        } else {
          insCondLog <- install_cond_env(
            envName = envName, useYML = yml,
            ymlFile = newYmlFile, pv = numPyVers3)
        }

        ## Error: folder exists but empty
        if ( any( grep('Could not solve for environment specs', insCondLog) ) & isTRUE(force)){
          envDir <- file.path(miniconda_path(), 'envs', envName)
          if ( dir.exists(envDir) ){
            conda_remove(envName)
            unlink( envDir, recursive = TRUE, force = TRUE )
            insCondLog <- install_cond_env(envName = envName, useYML = yml,
                                           ymlFile = newYmlFile,
                                           pv = numPyVers3)
          }
        }
        # + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "create" "--yes" "--name" "cola" "python=3.8" "--quiet" "-c" "conda-forge"
      } else {
        message(paste0('You need to run `conda_create("', envName,', python_version = "',
                       numPyVers3,'")` before using this package'))
        stop()
      }
    } else {
      colaexe <- condaLists$python[condaLists$name %in% envName]
      cat (sep = '', '    `', envName, '` conda environment installed in ', colaexe, '\n')
      (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
      (pyColaVersion <- tryCatch(system( paste0( pyCola, ' -V') , intern = TRUE), error = function(e) NULL))
      if ( pyColaVersion != 'Python 3.12.11' ){
        warning(paste0('We are using a new Python versiom. CoLs requires 3.12.11\n',
                       'To update your CoLa version, please run the following commands and try setup_cola() again:\n',
                       '\n\treticulate::conda_remove(envname = "cola") \n',
                       '\tfile.remove("', colaDir, '")\n',
                       '\treticulate::conda_update()\n',
                       '\tsetup_cola(pyVer = "3.12.11")\n'
        ))

        if (ask){
          user_permission <- utils::askYesNo(paste0("Uninstall '", envName, "' conda environment? Migth take some minutes"))
        }

        if ( isTRUE(user_permission) ) {
          cat (sep = '', '    Updating your Cola version to 3.12.11 ... please wait\n')
          reticulate::conda_remove(envname = envName)
          file.remove(dirname(pyCola))
          reticulate::conda_update()
          setup_cola(pyVer = "3.12.11")
        }
      }
    }
  }


  ## List conda after instaling
  (condaLists <- tryCatch(reticulate::conda_list(), error = function (e) NULL))

  ## Confirm env name
  (pyCola <- tryCatch( subset(condaLists, name == envName)$python, error = function (e) NULL) )
  if ( is.null(pyCola) | length(pyCola) == 0 ){
    message(paste0('You should run `conda_create("', envName ,'")` before using this package'))
    stop()
  } else {
    (pyColaVersion <- tryCatch(system( paste0( pyCola, ' -V') , intern = TRUE), error = function(e) NULL))
    if( is.character(pyCola) & file.exists(pyCola) & !is.null(pyColaVersion) ){
      cat (sep = '', '    `', envName, '` conda environment named correctly!\n')
      cat (sep = '', '    `', envName, '` using ', pyColaVersion,'\n')

    } else {

      ## Error -- cola/python.exe saved on paths but doesn't exists
      if (!file.exists(pyCola)){
        if (!file.exists(pyCola)){
          message(paste0(' Uninstalling corrupt previous installation'))
          tryCatch(conda_remove(envName))
          insCondLog <- install_cond_env(envName = envName, useYML = yml,
                                         ymlFile = newYmlFile, python_version = numPyVers2)
        }
      }

      ## Error -- no name of conda under "conda info --envs"
      (possiblePy <- grep(paste0(envName, '/python'), condaLists$python, value = TRUE))
      if( any(length(possiblePy)) ) {
        (minibat <- file.path(miniconda_path(), 'condabin/conda.bat'))
        if( file.exists(minibat) & !(envName %in% condaLists$name) ){
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


  ## Step4. Install packages ----------------------------------------------
  cat (sep = '', '  +Step 4/', nSteps, ' Installing & checking conda modules\n')

  if (!onlyIndividual){

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
  }

  ## Try individually

  ## list packages
  (avLibs <- reticulate::py_list_packages(envname = envName))
  head(avLibs)


  ## Install conda packages
  (libs2inst <- libs2Install[
    !((libs2Install %in% avLibs$package) |
        (gsub('==', '=', libs2Install) %in% avLibs$requirement))
  ])

  if( length(libs2inst) != 0 ){
    for( l in 1:length(libs2inst)){ # l = 10
      (lib2inst <- libs2inst[l])
      (lib2 <- gsub('==.+', '', lib2inst))
      (lib3 <- sub('=', '', lib2inst))
      # Check if specific
      (versReq <- ifelse(grepl('==', lib2inst), TRUE, FALSE))
      (versOK <- FALSE)
      if(versReq){
        (instVer <- avLibs$requirement[avLibs$package %in% lib2]) # instVer <- avLibs$requirement[avLibs$package %in% 'xx']
        if(length(instVer) != 0){
          (versOK <- ifelse(instVer %in% lib3, TRUE, FALSE)) # instVer == lib3
        }
      } else {
        versOK <- TRUE
      }

      if( (! lib2 %in% avLibs$package ) | !versOK ){ #
        cat(paste0(' \n --- Installing `',  libs2inst[l], '` module\n'))

        logPkg <- tryCatch(
          reticulate::py_install(
            envname = envName, # python_version = pyCola,
            channel = "conda-forge", packages = lib2inst),
          error = function (e) {print(e); NULL})

        # reticulate::py_install(envname = envName, channel = "conda-forge", packages = 'numexpr')
        # -- Installing `gdal` module + "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" "install" "--yes" "--name" "cola" "-c" "conda-forge" "gdal"

        # reticulate::conda_remove(envname = 'cola', packages = 'zarr')
        # "C:/Users/gonza/AppData/Local/r-miniconda/condabin/conda.bat" remove --yes --name cola zarr

        ## If there's a problem installing trought conda-forge, use PIP
        if( any(!is.null(logPkg)) ){
          cat( ' Trying to use PIP to install ', lib2inst, '\n')
          logPkg2 <- tryCatch(reticulate::py_install( envname = envName,
                                                      pip = TRUE,
                                                      packages = lib2inst),
                              error = function (e) {print(e); NULL})
        }
        avLibs <- reticulate::py_list_packages(envname = envName)
      }
    }
  }

  avLibs <- reticulate::py_list_packages(envname = envName)

  ## Installed
  insLibs <- libs2Install[
    ((libs2Install %in% avLibs$package) |
       (gsub('==', '=', libs2Install) %in% avLibs$requirement))
  ]

  ## No installed
  (noInsLibs <- libs2Install[
    !((libs2Install %in% avLibs$package) |
        (gsub('==', '=', libs2Install) %in% avLibs$requirement))
  ])
  # noInsLibs <- c('a', 'f')

  # if ( length(noInsLibs) == 1){
  #   if ( grepl(pattern = 'zarr', noInsLibs[1]) ) {
  #     cat('  You are missing zarr library in a specific version.',
  #         'This package is required to run very big and specific tasks on clusters.',
  #         'If this is not your intention, you can use this cola version. Otherwise, you',
  #         'need to remove the current "cola" conda environmnent, and reinstall it again.',
  #         'for this use this instructions:',
  #         '\t reticulate::conda_remove(envname = "cola")',
  #         '\t cola::setup_cola( forece = TRUE, dss = TRUE )',
  #         sep = '\n')
  #     noInsLibs <- grep('n', 'm')
  #   }
  # }

  if( length(noInsLibs) > 0 ){
    cat (sep = '', '  `', envName, '` conda environment installed!\n')
    cat(sep = '', "   Some python libraries aren't installed. Try \n\t`reticulate::py_install(envname = '", envName,
        "', channel = 'conda-forge', packages = c('",
        paste0(noInsLibs, collapse= "', '"), "')`")
    diagnose_cola()
    stop()
  } else {
    cat (sep = '', '    All required conda modules installed!\n')
  }


  ## Step4. Set paths ----------------------------------------------
  cat (sep = '', '  +Step 5/', nSteps, ' Setting up local variables\n')
  # reticulate::py_available()

  ## Connecting lib paths
  (welcomepy <- system.file("python/welcome.py", package = "cola"))
  (cola_scripts_path <- dirname(welcomepy))

  #(welcomepy <- file.path('N:/My Drive/git/cola/inst/python/welcome.py'))
  # read.delim(welcomepy)
  if (any(grep(' ', welcomepy))){
    welcomepy <- paste0('"', welcomepy, '"')
  }

  #tryA <- tryCatch(reticulate::py_exe(system.file("python/welcome.py", package = "cola")), error = function (e) e)
  (cmd2test <- paste0( #'cd ', cola_scripts_path, '; ',
    quotepath(pyCola), ' ', quotepath(welcomepy))); #cat(tryBcmd)
  (cmdans <- tryCatch( system( cmd2test , intern = TRUE ), error = function (e) e$message)) ## error is character


  ## Try to solve issues
  # C:\Users\Admin\AppData\Local\r-miniconda\envs\cola\Lib\site-packages\osgeo\_gdal.py

  if ( any(grep('failed|ImportError: DLL load failed while importing', cmdans)) ){

    # https://stackoverflow.com/questions/47246350/conda-activate-not-working
    # (instGd <- tryCatch(conda_install(envname = envName, packages = c('gdal', 'libgdal')),
    #                     error = function (e) e))
    # conda run -n envname python -c "print('Hello!')"

    # system(
    #   paste0('conda run --cwd ',
    #          cola_scripts_path,
    #          ' -n cola python -c "import osgeo; import rasterio;import os;import cola_functions as cf;print(os.getcwd()); print(1)"') )
    #

    # https://stackoverflow.com/questions/72142036/best-way-to-activate-my-conda-environment-for-a-python-script-on-a-windows-pc
    # conda run -n my_env python your_script.py
    # conda run -p /path/to/my_env python your_script.py


    (cmdans <- tryCatch(
      system(  paste0('conda run -n ', envName,' python ', welcomepy),
               intern = TRUE ), error = function (e) e$message)) # error us simpleError

    pyCola <- paste0('conda run --cwd ', cola_scripts_path, ' -n ', envName,' python ')

    if ( any(grep("'conda' not found", cmdans)) ){

      (pyCola <- paste0(reticulate::conda_binary(), ' run --cwd ', cola_scripts_path, ' -n ', envName,' python '))

      (cmdans <- tryCatch(
        system(  paste0(pyCola, ' ', welcomepy),
                 intern = TRUE ), error = function (e) e))

    }


    #   (cmd2testA <- paste0( 'conda run -n  '));
    #   (cmdansA <- tryCatch( system( cmd2testA , intern = TRUE ), error = function (e) e))
    #
    #
    #   (cmd2testA <- paste0( 'conda config --set auto_activate_base true'));
    #   (cmdansA <- tryCatch( system( cmd2testA , intern = TRUE ), error = function (e) e))
    #
    #   (cmd2testB <- paste0( 'activate ', envName));
    #   (cmdansB <- tryCatch( system( cmd2testB , intern = TRUE ), error = function (e) e))
    #
    #   (cmd2test <- paste0( pyCola, ' ', welcomepy));
    #   (cmdans <- tryCatch( system( cmd2test , intern = TRUE ), error = function (e) e))
    #
    #   (cmd2test <- paste0( #'cd ', cola_scripts_path, '; ',
    #     pyCola, ' ', welcomepy)); #cat(tryBcmd)
    #   (cmdans <- tryCatch( system( cmd2test , intern = TRUE ), error = function (e) e))
    #   (cmd2test <- paste0( 'activate ', envName,' ; ', pyCola, ' ', welcomepy)); #cat(tryBcmd)
    #   cat (sep = '', "\n    Some errors found when running the scripts.\n\t",
    #        "    Adding the following paths to te ENVIROMENTAL PATH:\n"
    #   )
    #
    #   (current_paths <- sort(strsplit(Sys.getenv('PATH'), ';')[[1]]))
    #
    #   ## Path according to
    #   ( gdal_data_path0 <- file.path( dirname(pyCola), "Library/share/gdal") );
    #   # "C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/Library/share/gdal"
    #
    #   if( dir.exists(gdal_data_path0) ){
    #     cat (sep = '', "\t + ", (gdal_data_path0), "\n")
    #     (cmd_add_gdal <- paste0("setx /m GDAL_DATA ",
    #                             gsub(fixed = T, "/", "\\", gsub("//", "/", gdal_data_path0))
    #     ))
    #     # cat(cmd_add_gdal)
    #     logA <- tryCatch(system(cmd_add_gdal, intern = TRUE), error = function(e) e)
    #     # PS C:\WINDOWS\system32> setx /m GDAL_DATA C:\Users\Admin\AppData\Local\r-miniconda\envs\cola\Library\share\gdal
    #     # CORRECTO: se guardó el valor especificado.
    #   }
    #
    #   ## Paths with dll files
    #   (gdal_data_paths <- list.files(path = dirname(pyCola), recursive = TRUE, pattern = 'gdal.+dll$', full.names = TRUE))
    #   # C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/Library/bin/gdal.dll
    #   for(ii in  1:length(gdal_data_path)){ # ii = 1
    #     gdal_data_path <- gdal_data_paths[ii]
    #     cat (sep = '', "\t + ", dirname(gdal_data_path), "\n")
    #     (cmd_add_gdal <- paste0("setx /m GDAL_DATA ",
    #                             (gsub(fixed = T, "/", "\\", gsub("//", "/", dirname(gdal_data_path)) ))
    #     ) )
    #     # setx /m GDAL_DATA C:\Users\Admin\AppData\Local\r-miniconda\envs\cola\Library\bin
    #     # CORRECTO: se guardó el valor especificado.
    #   }
    #
    # cat (sep = '', "\n    If any error detected, try to run in the command line (or power shell) as admin the instruction: \n",
    #      "\tsetx /m GDAL_DATA C:\\path\\mentioned\\before  -- Use only one BACKSLASH for paths separators\n")

  }


  if ( any(grep('WELCOME ', cmdans)) ){
    #libP <- .libPaths()
    #cola_scripts_path <- file.path(libP, 'cola/python')

    timeMark <- gsub('[[:punct:]]| ', '',
                     format(as.POSIXct(Sys.time(), tz="CET"),
                            tz="America/Bogota",usetz=TRUE))

    pyScript <- system.file("python/s2res.py", package = "cola")
    if ( any(grep(' ', pyScript)) ){
      pyScript <- paste0('"', pyScript, '"')
    }

    outTest <- paste0(tempfile(), timeMark, '.tif')
    # pyCola <- paste0('conda run --cwd ', cola_scripts_path, ' -n ', envName,' python ')

    # Create test
    (test_suit2res <- paste0(
      quotepath(pyCola), ' ', # Python
      quotepath(pyScript), ' ', # script
      quotepath(system.file("sampledata/sampleTif.tif", package = "cola")), ' ', #in [1]
      quotepath(outTest), #out
      ' 0 1 100', # min max scale-max
      ' 1 -9999 None')) # shape nodata proj
    # cat(test_suit2res)

    ## Run tests
    intCMD <- tryCatch(system(test_suit2res, intern = TRUE), error = function(e) e$message)


    ## Define local paths if test
    if ( dir.exists(cola_scripts_path) & file.exists(outTest) ){
      ## Success in testing the system


      ## Setting cola python as environmental variable
      Sys.setenv("COLA_PYTHON_PATH" = pyCola)
      Sys.setenv("COLA_SCRIPTS_PATH" = cola_scripts_path)
      options("COLA_PYTHON_PATH" = pyCola)
      options("COLA_SCRIPTS_PATH" = cola_scripts_path)

      ## Saving paths in .Renviron
      (home <- Sys.getenv("HOME"))
      renv <- file.path(home, ".Renviron")

      if ( file.exists(renv) ) {
        # Backup original .Renviron before doing anything else here.
        file.copy(renv, file.path(home, ".Renviron_backup"), overwrite = TRUE)
      }

      if ( !file.exists(renv) ) {
        file.create(renv)
      }

      {
        # if( !any(grep('COLA_PYTHON_PATH', origRenviron)) ){
        #
        #   ## Readfile
        #   con  <- file(renv, open = "r+")
        #   lines <- as.character()
        #   ii <- 1
        #
        #   while (TRUE) {
        #     line <- readLines(con, n = 1, warn = FALSE)
        #     if (length(line) == 0) {
        #       break()
        #     }
        #     lines[ii] <- line
        #     ii <- ii + 1
        #   }
        #
        #   # Set EARTHENGINE_PYTHON in .Renviron
        #   cola_python_line <- paste0('COLA_PYTHON_PATH="', pyCola, '"')
        #   cola_scripts_line <- paste0('COLA_SCRIPTS_PATH="', cola_scripts_path, '"')
        #
        #   system_vars <- c(lines, cola_python_line, cola_scripts_line)
        #
        #   writeLines(text = system_vars, con = con)
        #
        #   on.exit(close(con), add = TRUE)
        #   invisible(TRUE)
        #
        # } else {
        # }
      }

      origRenviron <- readLines(renv)
      Renviron <- origRenviron

      posA <- grep('COLA_PYTHON_PATH', Renviron)
      (posA <- ifelse(length(posA) == 0, length(Renviron) + 1, posA))
      Renviron[posA] <- paste0('COLA_PYTHON_PATH="', pyCola, '"')

      posB <- grep('COLA_SCRIPTS_PATH', Renviron)
      (posB <- ifelse(length(posB) == 0, length(Renviron) + 1, posB))
      Renviron[posB] <- paste0('COLA_SCRIPTS_PATH="', cola_scripts_path, '"')



      pos <- grep('COLA_DATA_PATH', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_DATA_PATH='}

      pos <- grep('COLA_NCORES', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_NCORES=1'}

      pos <- grep('COLA_DSS_UPL_MB', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_DSS_UPL_MB=250'}

      pos <- grep('COLA_VIZ_THRES_PIX', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_VIZ_THRES_PIX=1000000'}

      pos <- grep('COLA_VIZ_RES_NCOL', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_VIZ_RES_NCOL=1000'}

      pos <- grep('COLA_VIZ_RES_NROW', Renviron)
      # (pos <- ifelse(length(pos) == 0, length(Renviron) + 1, pos))
      if (length(pos) == 0){Renviron[length(Renviron) + 1] <- 'COLA_VIZ_RES_NROW=1000'}


      #cat(Renviron, sep = '\n')
      writeLines(text = Renviron, con = renv)
      #gcp <- pyCola
      #options('COLA_PYTHON_PATH')

      ## Checking gdal_calc.py


      # Sys.getenv(c("COLA_PYTHON_PATH", "COLA_SCRIPTS_PATH"))
      # Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH = new)
      # on.exit(Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH = old), add = TRUE)
      cat (sep = '', '\n\n\t=== Ready to connect landscapes! ===\n\n')

      if(dss){
        ## Step4. Set paths ----------------------------------------------
        cat (sep = '', '  + Extra step   Installing DSS GUI\n')
        cola::setup_cola_dss()
      }

      cat (sep = '', '\n\n',
           '\tCustomize your local parameteres by editing the file:\n\t',
           file.path(Sys.getenv("HOME"), ".Renviron"),'\n\n',
           '\n\tOpen it on R/Rstudio with the command:\n',
           '\tfile.edit(file.path(Sys.getenv("HOME"), ".Renviron"))\n\n',
           '\tPlease restart R to update the new settings\n'
      )

    } else {
      Sys.unsetenv("COLA_SCRIPTS_PATH")
      Sys.unsetenv("COLA_PYTHON_PATH")
      cat (sep = '', "    -- Final test didn't run. System vars COLA_SCRIPTS_PATH and COLA_PYTHON_PATH removed.\n\t", intCMD)
      diagnose_cola()
      stop()
    }
  }  else {
    cat (sep = '\n', "    Error: Can't load conda modules.\n\n\t Here the error: \n\t", cmdans)
  }
}

## Not run
## setup_cola()
## reticulate::conda_remove('cola')
## create --yes --name cola2 "python=3.12.11" zarr==3.1.2 networkit==11.1
## geopandas gdal h5py numexpr rasterio pytables pandas cython numba fiona shapely kdepy scikit-image kdepy --quiet -c conda-forge
## conda_create(envname = 'cola', python_version = '3.12.11', channel = 'conda-forge', packages = c('zarr', 'networkit', 'geopandas',
## 'gdal', 'h5py', 'numexpr', 'rasterio',  'pytables', 'pandas',  'cython', 'numba', 'fiona', 'shapely',  'kdepy', 'scikit-image'))
