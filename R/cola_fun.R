#' @title Launch \emph{COLA} decision support system (DSS) dashboard
#' @description This function runs the \emph{cola} dashboard
#' @param launch.browser Run this app on a new window? Default `TRUE`
#' @note Please see the official website (\url{https://github.com/connectingLandscapes/cola/})
#' @examples
#' library(cola)
#' if(interactive()) {
#' cola_dss()
#' }
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
cola_dss <- function(launch.browser = TRUE)  {
  #app_path <- system.file("shiny", package = "wallace")
  dssLocation <- system.file('app', package = "cola")
  #knitcitations::cleanbib()
  #options("citation_format" = "pandoc")
  preexisting_objects <- ls(envir = .GlobalEnv)
  on.exit(rm(list = setdiff(ls(envir = .GlobalEnv), preexisting_objects), envir = .GlobalEnv))
  return( shiny::runApp(dssLocation, launch.browser = launch.browser) )
}


#' @title Quote file path if space is detected
#' @description This function add double quotes to the string if a space is detected
#' @path The file path to be quoted
#' @examples
#' library(cola)
#' quotepath('C:/Users/First Second Name/Documents')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
quotepath <- function(path)  {
  #path <- 'N:/My Drive/connectivity-nasa-USFSIP/01_original-data/KAZA2025resamp/Leop_Clip_500m.tif'
  (newpath <- ifelse(grepl(" ", path), yes = paste0('"', path, '"'), path))
  return( path )
}


#' @title  Makes a population structure map from CDPOP results interpolation
#' @description Takes CDPOP results and generates a raster interpolation
#' @param py Python executable location
#' @param pyscript Python script location
#' @param grids String, List of CDPOP grid.csv files containing population genetic structure
#' @param template String. Path/filename of a template raster used for interpolation
#' @param method String. Interpolation method. One of 'multiquadric', 'thin_plate_spline', 'linear', or 'idw'.
#' @param neighbors String. Number of neighbors to use for interpolation. Either a positive integer < the total number of occupied points or 'all'.
#' @param crs String. User provided CRS as EPSG or ESRI string. Can also be 'None' in which case the CRS will be extracted from the template raster.
#' @return List with three slots: a) file, with NA if no result given or a single file if function was successful, b) newFiles, string vector with resulting files, c)log, string with the message obtained from the console execution. 0 indicates success.
#' @examples
#' cdpop_mapstruct( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

cdpop_mapstruct <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                            pyscript = system.file(package = 'cola', 'python/interpolate_popstructure.py'),
                            grids, template,
                            method = 'thin_plate_spline',
                            neighbors, crs = 'None', cml = TRUE, show.result = TRUE){

  # 1. List of CDPOP grid.csv files containing population genetic structure
  # 2. Path/filename of a template raster used for interpolation
  # 3. Interpolation method. One of 'multiquadric', 'thin_plate_spline', 'linear', or 'idw'.
  # 4. Number of neighbors to use for interpolation. Either a positive integer < the total number of occupied points or 'all'.
  # 5. User provided CRS as epsg or esri string. Can also be 'None' in which case the CRS will be extracted from the template raster.

  # neighbors = 5; method = 'thin_plate_spline'; crs = 'None'; py = Sys.getenv("COLA_PYTHON_PATH"); pyscript = system.file(package = 'cola', 'python/interpolate_popstructure.py')
  # grids = list.files(path = '/mnt/c/tempRLinux/RtmpNGsi0b/colaADC2024073113022305/cdpopoutLNY__1722460599', recursive = TRUE, pattern = 'grid.+csv$', full.names = TRUE)
  # template <- '/mnt/c/tempRLinux/RtmpNGsi0b/colaADC2024073113022305/out_surface_ZRY2024073113024005.tif'
  # grids <- grids[1]; read.csv(grids)

  # allele = '/mnt/c/tempRLinux/RtmpNGsi0b/colaADC2024073113022305/cdpopoutLNY__1722460599/allele.tif'
  # hetero = '/mnt/c/tempRLinux/RtmpNGsi0b/colaADC2024073113022305/cdpopoutLNY__1722460599/hetero.tif'

  if( !method %in% c('multiquadric', 'thin_plate_spline', 'linear', 'idw')){
    stop('Not valid method')
  }

  ### Create CMD
  (cmd_inter <- paste0(
    py, ' ', pyscript, ' ',
    grids, ' ', template, ' ',
    # allele, ' ', hetero, ' ',
    method, ' ', neighbors, ' ', crs, ' 2>&1'))
  if (cml){
    cat('\n\tCMD interpol struct: \n')
    cat(cmd_inter <- gsub(fixed = TRUE, '\\', '/', cmd_inter))
    cat('\n')
  }

  prevFiles <- list.files(path = dirname(grids), full.names = TRUE)
  intCMD <- tryCatch(system( cmd_inter ,
                             intern = TRUE),
                     error = function(e) e$message)
  #newFiles <- setdiff(list.files(path = dirname(grids), full.names = TRUE), prevFiles)
  newFiles <- grep(value = TRUE, pattern = 'heterozygosity.+.tif|alleles.+.tif',
                   list.files(path = dirname(grids), full.names = TRUE))

  # show.result = TRUE;
  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(any(file.exists(grep('tif', newFiles, value = TRUE))),
                             newFiles, ''),
               newFiles = newFiles,
               #log =  c(intCMD, read.delim('cdpop_mapstruct.txt')) ) )
               log = paste0("", intCMD) ) )

}


#' @title  Makes a population density map from CDPOP results interpolation
#' @description Takes CDPOP results and generates a raster interpolation
#' @param py Python executable location
#' @param cdpopscript Python location
#' @param grids String. Vector with the CDPOP output files to be interpolated
#' @param template String. Raster template for interpolating the files
#' @param method String. Method for interpolation. Options are: 'isj', 'silvermans', 'scotts', 'average', 'cv', 'user'
#' @param bandwiths String. Bandwidth value, either 'None' (default) or a list of user supplied bandwidths that can be converted from string to float or integer
#' @param type String. Specify whether count per cell or count per km2 should be returned. Options are: 'count', 'density'
#' @param crs String. User provided CRS as EPSG or ESRI string. Can also be 'None' in which case the CRS will be extracted from the template raster.
#' @return List with three slots: a) file, with NA if no result given or a single file if function was successful, b) newFiles, string vector with resulting files, c)log, string with the message obtained from the console execution. 0 indicates success.
#' @examples
#' cdpop_mapdensity( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

cdpop_mapdensity <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                             pyscript = system.file(package = 'cola', 'python/interpolate_popdensity.py'),
                             grids, template, method = 'average', bandwidths = 'None',
                             type = 'count', crs = 'None', cml = TRUE, show.result = TRUE){

  if(! method %in% c('isj', 'silvermans', 'scotts', 'average', 'cv', 'user')){
    stop("Not valid method: 'isj', 'silvermans', 'scotts', 'average', 'cv', 'user'")
  }

  if( !type %in% c('count', 'density')){
    stop('Not valid output: count or density')
  }
  logname <- paste0(tools::file_path_sans_ext(template), '_cdpop_mapdensity.txt')
  ### Create CMD
  (cmd_inter <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(grids), ' ',
    quotepath(template), ' ',
    # allele, ' ', hetero, ' ',
    method, ' ', bandwidths, ' ', type, ' ', crs
    , ' 2>&1 ' #, logname
  ))
  if (cml){
    cat('\n\tCMD interpol density: \n ')
    cat(cmd_inter <- gsub(fixed = TRUE, '\\', '/', cmd_inter))
    cat('\n')
  }

  prevFiles <- list.files(path = dirname(grids), full.names = TRUE)
  intCMD <- tryCatch(system( cmd_inter ,
                             intern = TRUE),
                     error = function(e) e$message)
  #newFiles <- setdiff(list.files(path = dirname(grids), full.names = TRUE), prevFiles)
  newFiles <- grep(value = TRUE, pattern = paste0(type, '.+', method, '.+'),
                   list.files(path = dirname(grids), full.names = TRUE))

  # , show.result = TRUE
  if(show.result){
    print(intCMD)
  }


  return( list(file = ifelse(any(file.exists(grep('tif', newFiles, value = TRUE))), newFiles, NA),
               newFiles = newFiles,
               log = paste0("", intCMD) ) )
  #log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
}



#' @title  Run CDPOP model
#' @description Run CDPOP model
#' @param py String. Python executable location
#' @param cdpopscript String. CDOPOP Python script location
#' @param inputvars String. CSV parameters file path location. Full path or relative to the CDOPOP Python script
#' @param agevars String. CSV age file path location. Full path or relative to the CDOPOP Python script
#' @param cdmat String. CSV cost-distance matrix file path location. Full path or relative to the CDOPOP Python script
#' @param xy String. CSV XY coordinates file path location. Full path or relative to the CDOPOP Python script
#' @param tempfolder String. Folder where results will be saved
#' @param prefix String. String added to output folder
#' @param cml String. Print the cola command line?. Default TRUE
#' @return List with two slots: a) newFiles with generated results, b) cdpopPath, with the folder with the results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

cdpop_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                     cdpopscript = system.file(package = 'cola', 'CDPOP/src/CDPOP.py'),
                     inputvars = NULL,
                     agevars = NULL,
                     cdmat = NULL,
                     xy = NULL,
                     tempFolder,
                     prefix = paste0('cdpopout', sessionIDgen(only3 = TRUE)),
                     cml = TRUE, show.result = TRUE){
  #xyfilename  NO .csv required
  #agefilename .csv required
  #matecdmat	cdmats/EDcdmatrix16
  #dispcdmat	cdmats/EDcdmatrix16

  # file.copy('invars.csv', '/home/user/cola/inst/examples/invars.csv')
  # file.copy('age.csv', '/home/user/cola/inst/examples/agevars.csv')
  if ( is.null(inputvars) ){
    cat(' // Using default CDPOP invars.csv\n')
    (inputvars <- system.file(package = 'cola', 'sampledata/invars.csv'))
  } else {
    cat(' // Using provided CDPOP ', inputvars,'\n')
  }

  if ( is.null(agevars) ){
    (agevars <- system.file(package = 'cola', 'sampledata/age.csv'))
  }

  if (is.null(cdmat)){
    stop(' CDmatrix required')
  }

  if (is.null(xy)){
    stop(' XY required')
  }

  datapath <- tempFolder # datapath = tempFolder
  timeMarkCDPOP <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE))
  (cdpopPath <- paste0(prefix, '__')) # timeMarkCDPOP
  logname <- paste0(tools::file_path_sans_ext(prefix), '_cdpop.txt')

  # inputvars <- read.csv('invars.csv')
  # agevars <- read.csv('age.csv')
  # xy <- read.csv('xy.csv')
  # dim(xy)
  # cdmat <- read.csv('cdmat.csv', header = FALSE)
  # dim(cdmat)

  # (inputvars <- system.file(package = 'cola', 'CDPOP/data/inputvars.csv'))
  # (agevars <- system.file(package = 'cola', 'CDPOP/data/age.csv'))
  # python CDPOP.py %userprofile%\dockerdata\CDPOP\data inputvars.csv outAnac1
  # file.copy(inputvars, paste0(datapath, '/inputvars.csv'), overwrite = TRUE)
  # file.copy('invars.csv', paste0(system.file(package = 'cola', 'sampledata'), '/invars.csv'))
  # file.copy('invars.csv','/home/user/cola/inst/sampledata/invars.csv')

  file.copy(cdmat, paste0(datapath, '/cdmat.csv'), overwrite = TRUE)
  file.copy(inputvars, paste0(datapath, '/invars.csv'), overwrite = TRUE)
  file.copy(agevars, paste0(datapath, '/age.csv'), overwrite = TRUE)


  (cmd <- paste0(
    quotepath(py), ' ',
    quotepath(cdpopscript), ' ',
    quotepath(datapath), ' invars.csv ', quotepath(cdpopPath)
    , ' 2>&1 ' #, logname
  ))
  if (cml){
    cat('\n\tCMD CDPOP: \n')
    cat(cmd, '\n')
  }

  CMDcp <- tryCatch(system(cmd, intern = TRUE), error = function(e) NULL)

  if(show.result){
    print(CMDcp)
  }

  newFiles0 <- list.files(path = datapath, recursive = TRUE, full.names = TRUE)
  (newFiles <- grep(pattern = cdpopPath, x = newFiles0, value = TRUE))

  ans2ret <- list(newFiles = newFiles, datapath = datapath,
                  cdpopPath = cdpopPath,
                  log =  CMDcp )
  #log =  paste0(intCMD, ' -- ', read.delim(logname)) )
  return(ans2ret)

  # (gridFiles <- grep(pattern = '/grid.+csv$', x = newFiles, value = TRUE))
  # gridFiles <- gridFiles[order( as.numeric(gsub('grid|.csv', '', basename(gridFiles) ) ) )]
  # gridLast <- read.csv(tail(gridFiles, 1))

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # inputvars <- read.csv('invars.csv')
  # agevars <- read.csv('age.csv')
  # xyorig <- read.csv('xy.csv'); dim(xy)

  # outfiles:
  #CDPOP/data/out11684988478/batchrun0mcrun0/grid0.csv [0, 1, ..]
  #CDPOP/data/out11684988478/batchrun0mcrun0/output.csv
  #CDPOP/data/out11684988478/batchrun1mcrun0/XY0.csv
  #CDPOP/data/out21684988758/cdpop.log
}


#' @title  shapefile to CDPOP XY
#' @description Converts Shapefile into xy file
#' @param shapefile String. Path to the shapefile
#' @param outxy String. Path to the resulting .xy file
#' @param tempDir String. Path to the working directory
#' @param mortrast String. Optional. Path to mortality/resistance raster. If it is a valid  raster, extracts the values to populate the 'Subpop_mortperc' field. The extracted values are scaled between 0 and 100.
#' @param survrast String. Optional. Path to survival/suitability raster. This argument will be ignored if a valid mortrast is provided. If it is a valid  raster, extracts the values to populate the 'Subpop_mortperc' field. The extracted values are flipped upside down and scaled between 0 and 100.
#' @param porcEmpty Integer. percentage of locations (points) that randomly will be considered empty and possible to be colonized.
#' @return String. Path to a temporal xy.csv file. Writes the same file at the outxy path
#' @examples
#' shp2xy( shapefile = 'shapefilepathhere.shp',outxy = 'out.xy', tempDir = 'temFolder')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

shp2xy <- function(shapefile, outxy, tempDir,
                   mortrast = NULL, survrast = NULL, porcEmpty = 0){

  # The n-(x,y) grid location values. This is a comma delimited file with 5 column headings:
  # (Subpopulation)- a unique identifier for  each individual corresponding to a unique subpopulation
  # (XCOORD)-x-coordinate location, (YCOORD)-y-coordinate location (YCOORD)
  # (ID)-a string label identifier, and
  # (sex)-an initial sex assignment (use 0/1 or F/M).
  # See xyED16.csv for an example xyfilename. The column order is necessary with a header file.


  # tempDir = "C:/temp/cola/colaXAQ2024111901384105"
  # shapefile = "C:/temp/cola//colaXAQ2024111901384105/out_simpts_RKF2024111901384805.shp"
  # mortrast = "C:/temp/cola//colaXAQ2024111901384105/out_surface_DDG2024111901385505.tif"
  # outxy = 'C:/temp/cola//colaXAQ2024111901384105/xy.csv'
  #
  # df <- foreign::read.dbf(gsub('\\..+', '.dbf', shapefile))
  # xyorig <- system.file(package = 'cola', 'sampledata/xy.csv')
  # shapefile<- system.file(package = 'cola', 'sampledata/xy.csv')

  #  shapefile<- 'C:/temp/cola/colaFFF2024111800273705/np_baseline_200points.shp'
  # survrast <- 'C:/temp/cola/colaSFY2024081405054705//in_surface_fixed_XMT2024081405055805.tif'
  # outxy <- 'C:/temp/cola/colaSFY2024081405054705//outxy.xy'
  # tempDir = "C:/temp/cola/colaSFY2024081405054705/"

  xy <- terra::vect(shapefile)

  xy$ID <- 1:nrow(xy)
  xy$Subpopulation <- 1#:nrow(xy)
  if ( any(!c('XCOORD', 'YCOORD') %in% names(xy))  ){
    xy_coord <- data.frame(geom(xy))
    xy$XCOORD <- xy_coord$x
    xy$YCOORD <- xy_coord$y
  }

  xy$Subpop_mortperc <- 0
  xy$X <- xy$XCOORD
  xy$Y <- xy$YCOORD

  xy$sex <- sample(x = c(1, 0), size = nrow(xy), replace = TRUE)
  xynew <- as.data.frame(xy)
  xynew <- xynew[, c('Subpopulation', 'X', 'Y', 'Subpop_mortperc', 'ID', 'sex')]
  xynew[, paste0('Fitness_', c('AA', 'Aa', 'aa', 'AABB', 'AaBB', 'aaBB', 'AABb', 'AaBb',
                               'aaBb', 'AAbb', 'Aabb', 'aabb'))] <- 0
  xynew$Fitness_AA <- 50
  xynew$Fitness_aa <- 100
  xynew$Fitness_Aa <- 16

  # porcEmpty <- 50
  if (porcEmpty > 0){
    nEmpty <- round( (porcEmpty/100) * nrow(xy))
    posEmpty <- sample(x = 1:nrow(xy), size = nEmpty, replace = FALSE)
    xynew$ID[posEmpty] <- 'NA'
    xynew$sex[posEmpty] <- 'NA'
  }

  if (!is.null(mortrast) | !is.null(survrast) ){
    vals2add <- 0

    if (!is.null(mortrast)){
      rast_path <- mortrast
      mort <- TRUE
    } else if( !is.null(survrast)){
      rast_path <- survrast
      mort <- FALSE
    }

    if(file.exists(rast_path)){
      # rast_path <- 'C:/temp/cola/colaFFF2024111800273705/in_dist_fixed_MBU2024111800275205.tif'
      rcdpop <- tryCatch(terra::rast(rast_path), error = function(e) NULL)
      names(rcdpop) <- 'rcdpop'
      if(!is.null(rcdpop)){
        cat('  Extracing raster values for mortality\n')
        extVals <- terra::extract(rcdpop, xy[, c('XCOORD','YCOORD')])
        extVals <- extVals2 <- extVals$rcdpop # [, 2]
        #extVals <- (20:200)
        rng <- range(extVals, na.rm = TRUE)
        if( any(!is.na(rng)) & length(unique(rng)) > 1 ){
          cat('    --- At least some values OK \n')
          extVals2 <- (extVals - min(rng))/(max(rng) - min(rng))
          if(!mort){
            extVals2 <- max(extVals2) - extVals2
          }

          vals2add <- xy$Subpop_mortperc <- extVals2 * 100

          if (!is.null(mortrast)){
            if(file.exists(mortrast)){
              #rx <- terra::rast(mortrast)
              #exct <- (terra::extract(rx, xynew[, c( 'X', 'Y')]))
              #vect2add <- as.numeric(exct[, 2])
              #vect2add <- vect2add/max(vect2add)*100
              vals2add <- vals2add
              xynew$Subpop_mortperc <- vals2add
            }
          } else if ( !is.null(survrast) ){
            if ( file.exists(survrast)) {
              #rx <- terra::rast(survrast)
              #exct <- (terra::extract(rx, xynew[, c( 'X', 'Y')]))
              #vect2add <- as.numeric(exct[, 2])
              vals2add <- (max(vect2add, na.rm = TRUE) - vect2add)/max(vect2add)*100
              xynew$Subpop_mortperc <- vals2add
            }
          }
        }
      }
    }
    xynew$Subpop_mortperc <- vals2add
  }

  # print(head(xynew))

  # ## Write CSV
  write.csv(x = xynew, file = outxy, row.names = FALSE, quote = FALSE)

  ## Write CSV as xy.csv
  # xynew[1:3, 1:3]
  xy_out <- file.path(tempDir, 'xy.csv')
  write.csv(x = xynew, file = xy_out, row.names = FALSE, quote = FALSE)

  return(xy_out)

  # xy$YCOORD <- 1:nrow(xy)
  # # xyorig <- read.csv('xy.csv'); dim(xy)
  # xynew <- xyorig[sample(x = 1:nrow(xyorig), replace = TRUE, size = nrow(df)), ]
  # xynew$Subpopulation <- 1
  # xynew$XCOORD <- df$X
  # xynew$YCOORD <- df$Y
  # xynew$ID <- 1:nrow(df)
  # xynew$sex <- sample(x = c(0, 1), replace = TRUE, size = nrow(df))
  # xynew$age <- 1
  # #dim(xynew)
  # xynew[is.na(xynew)] <- 1#"NA"
}


#' @title  Get NA value from raster layer metadata
#' @description Extracts the NoData value from raster metadata
#' @param path File location
#' @return String with the NA value. No converted to number to avoid characters loss
#' @examples
#' input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' guessNoData(input_tif)
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
guessNoData <- function(path){
  # path <- input_tif
  ans <- -9999
  if (require(gdalUtilities) & file.exists(path)){
    gi <- strsplit(gdalUtilities::gdalinfo(path, quiet = TRUE), '\n')[[1]]
    ndv <- grep('NoData ', gi, value = TRUE)
    if( any(length(ndv)) ) {
      (ans <- as.numeric(gsub('.+\\=', '', ndv)))
    }
  }
  if(is.na(ans) | ans == 'nan' | ans == '' | length(ans) == 0){
    ans <- -9999
  }
  return(ans)
}



#' @title  Get MinMax value from raster layer
#' @description Extracts the min and max value from raster
#' @param rastPath File location or raster layer. Accepts file path or SpatRaster object
#' @return Numeric vector of two numbers
#' @examples
#' input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' getMnMx(input_tif)
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
getMnMx <- function(rastPath, na.rm = TRUE){
  # rastPath0 <- rastPath
  # rastPath <- rastPath0
  #rastPath = '/data/tempR/colaEXR2025012814493505/Mountainaire_prePercentCanopy30m.tif'
  # rastPath = 'C:/temp/cola/colaOCC2025011423394005//out_lcc_JCH2025011423414805.tif'
  na.rm * 1

  ra <- c(NA, NA)
  # Path path
  if(class(rastPath) == 'character'){
    invisible(ras <- sf::gdal_utils('info', rastPath,  options = c('-mm'), quiet = TRUE))
    #ra2 <- sf::gdal_utils('info', rastPath,  options = c('-stats'))
    ra <- as.numeric(strsplit(split = ',',
                              gsub('.+=', '', grep(strsplit(ras, '\n')[[1]], pattern = 'Min/Max', value = TRUE))
    )[[1]])
    #rst <- terra::rast(rastPath)

    if( !(all(is.numeric(ra)) & all(!is.infinite(ra)))  ){
      rst <- terra::rast(rastPath)
      ra <- setMinMax(rst, force=TRUE)
      ra <- minmax(rst)[1:2]
      if( !(all(is.numeric(ra)) & all(!is.infinite(ra))) ){
        ra <- global(rst, 'range' )
        if(!(all(is.numeric(ra)) & all(!is.infinite(ra))) ){
          (ra <- range(rst[], na.rm = TRUE))
        }
      }
    }

    # Raster path
  } else if (class(rastPath) == 'SpatRaster'){
    rst <- (rastPath)
    rastPath <- sources(rst)
    if(file.exists(rastPath)){ # in disk
      invisible(ras <- sf::gdal_utils('info', rastPath,  options = c('-mm'), quiet = TRUE))
      ra <- as.numeric(strsplit(split = ',',
                                gsub('.+=', '', grep(strsplit(ras, '\n')[[1]], pattern = 'Min/Max', value = TRUE))
      )[[1]])
    } else { # in memory
      ra <- minmax(rst)[1:2]
    }
  }
  return(ra)
}
#getMnMx(rastPath=rast(volcano))
#getMnMx('/home/shiny/Probability.tif')
# cola::getMnMx('/data/tempR//colaRDB2025012818465505//out_surface_IEB2025012818473305.tif')


#' @title  Adapt file path. Change backslash to slash
#' @description Fix paths so internal console recognize paths. Change double backslash to slash
#' @param path File location
#' @return File location String using slash as separator
#' @examples
#' newPath <- adaptFilePath('\\temp\\path\\to\\file\\here.tif')
#' newPath
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
adaptFilePath <- function(path){
  # path = "C:\\Users\\ig299\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/python/s2res.py C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/sampledata/sampleTif.tif C:\\Users\\ig299\\AppData\\Local\\Temp\\RtmpwrUyVu/VM2024041715525605file51c6028258//out_surface_JQ2024041715525705file51c260d2c4b.tif 0.06788435 0.9989325 100 1 -9999 None"
  return ( gsub(fixed = TRUE, '\\', '/', path) )
}



#' @title  Suitability to resistance
#' @description Transforms suitability to resistance surface
#' @param py String. Python location or executable. The string used in R command line to activate `cola`. The default versio should point to a conda environment. Might change among computers versions
#' @param pyscript String. Python script location
#' @param intif String. File path to the input raster.
#' @param outtif String. File path of the output surface resistance
#' @param minval Numeric. The lower value on the input raster to cut off. Pixels with values under the given number will be ignored. In the front end the values automatically derived from the input file.
#' @param maxval Numeric. The upper value on the input raster to cut off. Pixels with values under the given number will be ignored. In the front end the values automatically derived from the input file.
#' @param maxout Numeric. This is the maximum resistance value after transformation from suitability. Default value is 100.  Minimum value is set to 1.
#' @param shape Numeric. A statistical parameters that defines the transformation pattern between the input and output. The shape value determines the relationship between suitability and resistance. For a linear relationship, use a value close to 0, such as 0.01. Positive values result in a greater increase in resistance as suitability declines. This is appropriate for animals that are more sensitive to the matrix in between habitat. Negative values result in a lesser increase in resistance as suitability declines. This is appropriate for animals that are less sensitive to the matrix between habitat. The more positive or more negative, the greater the effect on the shape of the relationship. Values generally range between +10 and -10, 0 is not allowed.
#' @param nodata Numeric. The no data value of the input file. For GeoTiffs, this is automatically determined. For text based files, this must be input by the user. Default value is ‘None’.
#' @param prj string. Projection information in the case the input raster [1] has no spatial projection. For GeoTiffs, this is automatically determined. For text based files, this must be input by the user. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.
#' @param cml Logical. Print the back-end command line? Default  TRUE
#' @param show.result Logical. Print the command line result? Default  TRUE
#' @return Creates a raster layer with a minimum value of 1 and maximum value given the parameter 5. The internal R object is a list of two slots. The first one contains the path of the created raster, if any, and the second slot includes any function message or log, if any.
#' @examples
#' hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' srp30 <- sui2res_py(intif = hs, outtif = 'hs_p3_0.tif', #' minval = 0, maxval = 1, maxout = 10,  shape = 1,  nodata = NULL, prj = 'None')
#' hs_rast <- terra::rast(hs)
#' plot(hs_rast, main = 'Habitat suitability')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
sui2res_py <- function(intif, outtif,
                       minval, maxval, maxout, shape,
                       nodata = NULL, prj = 'None',
                       py = Sys.getenv("COLA_PYTHON_PATH"),
                       pyscript = system.file(package = 'cola', 'python/s2res.py'),
                       cml = TRUE, show.result = TRUE){
  # minval = 0
  # maxval =  100
  # maxout = 100
  # shape = 1
  # nodata = -9999
  # prj = 'None'

  ## Guess NA value if not provided
  if(is.null(nodata)){
    nodata <- guessNoData(intif)
  }

  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')

  ## Create CMD
  (cmd_s2res <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(minval, scientific=F), ' ',
    format(maxval, scientific=F), ' ',
    format(maxout, scientific=F), ' ',
    format(shape, scientific=F), ' ',
    format(nodata, scientific=F), ' ',
    prj
    , ' 2>&1 ' #, logname
  ))
  if (cml){
    cat('\n\tCMD Surface : \n')
    cat(cmd_s2res <- gsub(fixed = TRUE, '\\', '/', cmd_s2res))
    cat('\n')
  }

  intCMD <- tryCatch(system( cmd_s2res , intern = TRUE),
                     error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(file.exists(outtif), outtif, ''),
               #log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
               log = paste0("", intCMD) ) )
}

#' @title  Create random points
#' @description Create random points
#' @param rvect Raster. Raster object co be sampled
#' @param npts Numeric. Number of points
#' @param rmin Numeric. Raster minimum values o¿to consider
#' @param rmax Numeric. Raster maximun value to consider
#' @return Path with CDPOP results
#' @examples
#' hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' newPoints <- randPtsFun(terra::rast(hs), 10, 0.2, 0.8)
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
randPtsFun <- function(rvect, npts, rmin, rmax){
  # rvect <- hs_rast[]
  # npts = 50; rmin = .20; rmax = .5

  ## Only positions among range
  newRast <- which(rvect > rmin & rvect < rmax)

  ## Selected positions
  randPos <- sample(newRast, size = npts)

  ## Extracting xy from selected cells
  xyRand <- data.frame(xyFromCell(x, cell = randPos))

  ## Extracting values from selected cells
  # xyRand$val <- terra::extract(x = x, y = randPos)[, 1]
  # xyRand$val <- rvect[randPos]

  return(xyRand)
}


#' @title  Simulate spatial points
#' @description Creates a dispersal resistance matrix among a set of points. create_cdmat.py in the command line or cdmat_py( ) in R.
#' @param py String. Python location or executable. The string used in R command line to activate `cola`. The default version should point to a conda environment. Might change among computers versions
#' @param pyscript String. Python script location
#' @param inshp  String. Source points file path to the point layer with no spaces. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
#' @param intif String. Surface resistance input raster with no spaces. Requires a projected file with square pixels. Not LonLat projection allowed
#' @param outtif String. Output point layer file path, with no spaces. Written in ESRI Shapefile format.
#' @param minval Numeric. Minimum value. The lower value of the pixels in the raster to consider to simulate the points.
#' @param maxval Numeric. Maximum value. The upper value of the pixels in the raster to consider to simulate the points.
#' @param npoints Integer. Number of points. Number of points to simulate.
#' @param issuit String. Is it suitable? ‘Yes’ (default) or ‘No’. Indicates if the provided raster [1]  is suitability. If so, the script will likely sample higher value pixels. If ‘No’, will assume it is resistance and will sample more likely lower values
#' @param upcrs String. Update CRS | | upcrs| |String| |Projection information in the case the input raster [1] has no spatial projection. For GeoTiffs, this is automatically determined. For text-based files like ASCII or RSG rasters, this must be input by the user. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with the created shapefile
#' @examples
#' library(cola)
#' library(terra)
#' hs_path <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' # hs_path <- 'C:/path/to/raster.tif'
#' points_path <- system.file(package = 'cola', 'sampledata/samplePoints.shp')
#' # points_path <- 'C:/path/to/points.shp'
#' pts_result <- points_py(intif = hs_path, outshp = 'out_pts.shp', minval = 0.2, maxval = 0.9, npoints = 50, issuit = 'Yes', upcrs = 'None')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
points_py <- function(intif, outshp,
                      smin, smax, npoints, issuit = 'Yes', upcrs = 'None',
                      py = Sys.getenv("COLA_PYTHON_PATH"),
                      pyscript = system.file(package = 'cola', 'python/create_source_points.py'),
                      cml = TRUE, show.result = TRUE){
  # smin = 2
  # smax =  95
  # npoints = 50
  # pyscript <- system.file(package = 'cola', 'python/create_source_points.py'

  logname <- paste0(tools::file_path_sans_ext(outshp), '.txt')


  (cmd_pts <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(intif), ' ',
    quotepath(outshp), ' ',
    format(smin, scientific=F), ' ',
    format(smax, scientific=F), ' ',
    format(npoints, scientific=F), ' ',
    issuit, ' ', upcrs
    , ' 2>&1 '# , logname
  ))

  if (cml){
    cat('\n\tCMD Points: \n')
    cat(cmd_pts <- gsub(fixed = TRUE, '\\', '/', cmd_pts))
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_pts, intern = TRUE), error = function(e) e$message)
  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(file.exists(outshp), outshp, ''),
               # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
               log = paste0("", intCMD) ) )


}

#' @title  Creates CDmatrix for CDPOP model
#' @description Creates a dispersal resistance matrix among a set of points. create_cdmat.py in the command line or cdmat_py( ) in R.
#' @param inshp String. Source points File path to the point layer. Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
#' @param intif String. Surface resistance File path to the input raster. Requires a GeoTIFF file with square pixels
#' @param outcsv String. Output csv file name. Path of the output csv matrix
#' @param maxdist Numeric. Max. dispersal distance in cost units. This is the maximum distance to consider when calculating kernels and should correspond to the maximum dispersal distance of the focal species. Values greater than this will be converted to 0 before summing kernels. For example, if the maximum dispersal distance of the focal species is 10 km, set this value to 10000.
#' @param ncores Numeric. Number of cores. Number of CPU cores to run the analysis
#' @param crs String. Projection string. String. Projection information in the case the input raster 'intif' has no spatial projection. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.
#' @param intif String.
#' @param py Python executable location
#' @param pyscript Python script location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with the CSV matrix
#' @examples
#' library(cola)
#' hs_path <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' # hs_path <- 'C:/path/to/raster.tif'
#' points_path <- system.file(package = 'cola', 'sampledata/samplePoints.shp')
#' # points_path <- 'C:/path/to/points.shp'
#' mat_result <- cdmat_py(inshp = points_path, intif = hs_path, outtif = 'out_mat.csv')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

cdmat_py <- function(inshp, intif, outcsv,
                     maxdist,
                     ncores = Sys.getenv("COLA_NCORES"),
                     crs = 'None',
                     py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/create_cdmat.py'),
                     cml = TRUE, show.result = TRUE){
  # maxdist = 100000
  # create_cdmat.py
  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)
  # [5] cores
  # [6] Projection # User provided CRS if using ascii or other file without projection info, as epsg or esri string e.g. "ESRI:102028"Default None

  # inshp =  '/tmp/RtmpiD9uhw/colaMJJ2024070817062505/out_simpts_XYA2024070817112705.shp'
  # intif =  '/tmp/RtmpiD9uhw/colaMJJ2024070817062505/out_surface_HWP2024070817115305.tif'
  # outcsv = '/tmp/RtmpiD9uhw/colaMJJ2024070817062505/out_cdmatrix_KQI.csv'
  # maxdist = 100000; param5 = 1; param6 = 'None'
  logname <- paste0(tools::file_path_sans_ext(outcsv), '.txt')

  # pyscript <- system.file(package = 'cola', 'python/create_cdmat.py')
  (cmd_cdmat <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outcsv),
    ' ',
    format(maxdist, scientific=F), ' ',
    ncores, ' ', crs
    , ' 2>&1 '
    # , logname
  ))
  if (cml){
    cat('\n\n\tCMD cdmat: \n')
    cat(cmd_cdmat <- gsub(fixed = TRUE, '\\', '/', cmd_cdmat))
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_cdmat, intern = TRUE), error = function(e) e$message)
  #checkcsv <- read.csv(outcsv); which(is.numeric(checkcsv)) ; summary(checkcsv); sum(checkcsv, )

  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(file.exists(outcsv), outcsv, ''),
               # log =  paste0(intCMD, ' -- ', read.delim(logname) ) ) )
               log = paste0("", intCMD) ) )
}


#' @title  Corridors
#' @description Calculate corridors
#' @param py Python executable location
#' @param pyscript Python script location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lcc_py <- function(inshp, intif, outtif,
                   maxdist, smooth, tolerance,
                   ncores = as.numeric(Sys.getenv('COLA_NCORES')), crs = 'None',
                   py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/lcc.py'),
                   cml = TRUE, show.result = TRUE){
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] Max. dispersal distance (meters)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)

  # if(is.null(param8)){
  #   param8 <- guessNoData(intif)
  # }
  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')

  (cmd_lcc <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(maxdist, scientific=F), ' ',
    format(smooth, scientific=F), ' ',
    format(tolerance, scientific=F), " ",
    format(ncores, scientific=F), " ",
    crs
    , ' 2>&1 ' #, logname
  ))

  if (cml){
    cat('\n\tCMD LCC: \n')
    cat(cmd_lcc <- gsub(fixed = TRUE, '\\', '/', cmd_lcc))
    cat('\n')
  }


  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
  return( ans )
}


#' @title  Create least cost corridors for heavy rasters
#' @description Run CDPOP model
#' @param py Python executable location
#' @param pyscript Python script location
#' @param py Python executable location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lccHeavy_py <- function(inshp, intif, outtif,
                        maxdist, smooth, tolerance,
                        ncores = as.numeric(Sys.getenv('COLA_NCORES')),
                        crs = 'None', tempFolder = NULL,
                        py = Sys.getenv("COLA_PYTHON_PATH"),
                        pyscript = system.file(package = 'cola', 'python/lcc_heavy.py'),
                        cml = TRUE, show.result = TRUE){

  # "lcc_hdf5_v6.py" "pts.shp inraster.tif out.tif 10000000 0 1000 6 None first.h5 second.h5 rmlimitinGB"
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (should be in meters*)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)
  # [7] number of cores
  # [8] projection if missing
  # [9] first h5 temp file
  # [10] second h5 temp file
  # [11] Max GB ram allowed

  if (is.null(tempFolder)){
    tempFolder <- tempdir()
  }
  tempH5 <- basename(tempfile())
  h5file1 <- paste0(tempFolder, '/', tempH5, '_A.h5')
  h5file2 <- paste0(tempFolder, '/', tempH5, '_B.h5')

  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')

  (cmd_lcc <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(maxdist, scientific=F), ' ',
    format(smooth, scientific=F), ' ',
    format(tolerance, scientific=F), " ",
    format(ncores, scientific=F), " ",
    crs, " ",
    h5file1, " ",
    h5file2, " ",
    '50'
    , ' 2>&1 ' #, logname
  ))

  (cmd_lcc <- gsub(fixed = TRUE, '\\', '/', cmd_lcc))

  if (cml){
    cat('\n\tCMD LCC:\n', cmd_lcc)
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE), error = function(e) e$message)

  tryCatch(file.remove(c(h5file1, h5file2)), error = function(e) NULL)

  if(show.result){
    print(intCMD)
  }

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
  # print('ANS LCC');print(ans)
  return( ans )

}



#' @title  Create least cost corridors using parallel computing
#' @description Run CDPOP model
#' @param py Python executable location
#' @param pyscript Python script location
#' @param py Python executable location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lccJoblib_py <- function(inshp, intif, outtif,
                         maxdist, smooth, tolerance,
                         ncores = as.numeric(Sys.getenv('COLA_NCORES')),
                         crs = 'None',
                         maxram = 6,
                         tempFolder = NULL,
                         py = Sys.getenv("COLA_PYTHON_PATH"),
                         pyscript = system.file(package = 'cola', 'python/lcc_joblib.py'),
                         cml = TRUE, show.result = TRUE){

  # "lcc_hdf5_v6.py" "pts.shp inraster.tif out.tif 10000000 0 1000 6 None first.h5 second.h5 rmlimitinGB"
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (should be in meters*)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)
  # [7] number of cores
  # [8] projection if missing
  # [9] first h5 temp file
  # [10] second h5 temp file
  # [11] Max GB ram allowed

  if (is.null(tempFolder)){
    tempFolder <- tempdir()
  }
  tempH5 <- basename(tempfile())
  h5file1 <- paste0(tempFolder, '/', tempH5, '_A.h5')
  h5file2 <- paste0(tempFolder, '/', tempH5, '_B.h5')

  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')

  (cmd_lcc <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(maxdist, scientific=F), ' ',
    format(smooth, scientific=F), ' ',
    format(tolerance, scientific=F), " ",
    format(ncores, scientific=F), " ",
    crs, " ",
    h5file1, " ",
    h5file2, " ",
    maxram
    , ' 2>&1 ' #, logname

  ))
  (cmd_lcc <- gsub(fixed = TRUE, '\\', '/', cmd_lcc))

  if (cml){
    cat('\n\tCMD LCC joblib:\n', cmd_lcc)
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  tryCatch(file.remove(c(h5file1, h5file2)), error = function(e) NULL)

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
  #print('ANS LCC');print(ans)
  return( ans )
}



#' @title  Create least cost corridors using zarr
#' @description Create a corridor raster using zarr. This approach is used for big data analysis
#' @param inshp String. Path to file holding xy coordinates
#' @param intif String. Path to resistance grid
#' @param outtif String. Path to corridor GeoTIFF result
#' @param maxdist Numeric. Distance threshold
#' @param smooth Numeric. Radius for gaussian smoother (in number of cells). The size of the kernel on each side is 2*radius + 1. E.g. a radius of 2 gives a 5x5 cell kernel
#' @param tolerance Numeric. Amount to add to the least cost path (in cost distance units) in order to generate a swath of low cost pixels, termed the least cost corridor. If 0, returns the least cost path. If > 0, this amount is added to the least cost path value so that all pixels with values <= to that value will be returned. This results in a swath of pixels instead of a single pixel wide path.
#' @param ncores Numeric. Number of threads to use
#' @param maxram Numeric. Set memory size for processing corridors, i.e. set to 16 if you want to use 16GB of RAM, when processing. Make sure you have enough RAM available when setting this value. Consider the total amount of RAM available on your computer and the amount used by other programs that may be running.
#' @param crs String. User provided CRS if using ascii or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"
#' @param sci Numeric. Default is 'None'. Start corridor index. For now, these should be zero indexed python style. E.g. for a landscape with 10,000 corridors a first batch of corridors could be 0-500. Python range is such that this would process corridors 0-499. Then next batch would be 500-1000, which would process corridors 500-999. The next batch would be 1000-1500, and so on.
#' @param eci Numeric. End corridor index. Default is 'None'
#' @param tempFolder Numeric.
#' @param py Python executable location
#' @param pyscript Python script location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' lccZarr_py( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lccZarr_py <- function(inshp, intif, outtif,
                       maxdist, smooth, tolerance,
                       ncores = as.numeric(Sys.getenv('COLA_NCORES')),
                       maxram = 6, crs = 'None',
                       sci = 'None', eci = 'None',
                       tempFolder = NULL,
                       py = Sys.getenv("COLA_PYTHON_PATH"),
                       pyscriptA = system.file(package = 'cola', 'python/lcc_hpc1_zarr.py'),
                       pyscriptB = system.file(package = 'cola', 'python/lcc_hpc2_zarr.py'),
                       cml = TRUE, show.result = TRUE){

  # inshp = system.file(package = 'cola', 'sampledata/points_sabah_50.shp');
  # intif = '/home/shiny/cola/inst/sampledata/sampleSR.tif';
  # outtif = '/home/shiny/test/testzarr.tif';
  # maxdist = 1000000; smooth = 0; tolerance = 0 ;
  # ncores = 8; maxram = 6;
  # crs = 'None'; sci = 'None'; eci = 'None'
  # tempFolder = NULL;
  # tempFolder = '/home/shiny/test/'

    # A: inshp intif outtif maxdist smooth tolerance ncores crs pazarr dazarr reOrderFile nodeidsFile maxram sci eci


  if (is.null(tempFolder)){
    #tempFolder <- tempdir()
    tempFolder <- dirname(intif)
  }

  (tempH5 <- basename(tempfile()))
  (pazarr <- paste0(tempFolder, '/', tempH5, '_pazarr.zarr'))
  (dazarr <- paste0(tempFolder, '/', tempH5, '_dazarr.zarr'))
  (reOrderFile <- paste0(tempFolder, '/', tempH5, '_reOrderFile.csv'))
  (nodeidsFile <- paste0(tempFolder, '/', tempH5, '_nodeidsFile.csv'))

  # # INPUTS part A
  # # Path to file holding xy coordinates
  # xyf = sys.argv[1]  --- inshp
  #
  # # Path to resistance grid
  # rg = sys.argv[2] --- intif
  #
  # # Distance threshold
  # dThreshold = sys.argv[3]
  #
  # # Number of threads to use
  # nThreads = sys.argv[4] # Default 1
  #
  # # User provided CRS if using ascii or other file without projection info
  # # Provide as epsg or esri string e.g. "ESRI:102028"
  # upCRS = sys.argv[5] # Default None
  #
  # # Output pairwise point hdf name
  # # This is a temporary file and gets deleted at the end of the script,
  # # unless there's a script failure.
  # ppzarr = sys.argv[6]
  #
  # # Output distance array hdf name
  # # This is a temporary file and gets deleted at the end of the script,
  # # unless there's a script failure.
  # dazarr = sys.argv[7]
  #
  # # point pair output file name
  # reOrderFile = sys.argv[8]
  #
  # # Node ids output
  # nodeidsFile = sys.argv[9]
  #
  # # Set memory size for processing corridors
  # # I.e. set to 16 if you want to use 16GB of RAM
  # # when processing. Make sure you have enough RAM
  # # available when setting this value. Consider
  # # the total amount of RAM available on your computer
  # # and the amount used by other programs that may
  # # be running.
  # gbLim = sys.argv[10] # Default 6

  # INPUTS part 2 ---------------
  # # Path to resistance grid
  # rg = sys.argv[1]
  #
  # # Output file path
  # ofile = sys.argv[2]
  #
  # # Zarr file name
  # dazarr = sys.argv[3]
  #
  # # Radius for gaussian smoother (in number of cells)
  # # The size of the kernel on each side is 2*radius + 1
  # # E.g. a radius of 2 gives a 5x5 cell kernel
  # gRad = sys.argv[4]
  #
  # # Amount to add to the least cost path (in cost distance units)
  # # in order to generate a swath of low cost pixels,
  # # termed the least cost corridor.
  # # If 0, returns the least cost path.
  # # If > 0, this amount is added to the least cost path
  # # value so that all pixels with values <= to that value
  # # will be returned. This results in a swath of pixels
  # # instead of a single pixel wide path.
  # corrTolerance = sys.argv[5]
  #
  # # Number of threads to use
  # nThreads = sys.argv[6] # Default 1
  #
  # # User provided CRS if using ascii or other file without projection info
  # # Provide as epsg or esri string e.g. "ESRI:102028"
  # upCRS = sys.argv[7] # Default None
  #
  # # point pair output file name
  # reOrderFile = sys.argv[8]
  #
  # # nodeids file name
  # nodeidsFile = sys.argv[9]
  #
  # # Start corridor index
  # # For now, these should be zero indexed python style
  # # E.g. for a landscape with 10,000 corridors
  # # a first batch of corridors could be 0-500
  # # Python range is such that this would process
  # # corridors 0-499. Then next batch would be 500-1000,
  # # which would process corridors 500-999. The next
  # # batch would be 1000-1500, and so on.
  # sci = sys.argv[10] # Default is None
  #
  # # End corridor index
  # eci = sys.argv[11] # Default is None

  (logname <- paste0(tools::file_path_sans_ext(outtif), '.txt'))


  (cmd_lcc_zarrA <- paste0(
    quotepath(py), ' ',
    quotepath(pyscriptA), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    format(maxdist, scientific=F), ' ',
    format(ncores, scientific=F), " ",
    crs, " ",
    pazarr, " ",
    dazarr, " ",
    reOrderFile, " ",
    nodeidsFile, " ",
    maxram
    , ' 2>&1 ' #, logname
  ))

  # A: inshp intif maxdist ncores crs pazarr dazarr reOrderFile nodeidsFile maxram
  # B: intif outtif pazarr smooth tolerance ncores crs reOrderFile nodeidsFile sci eci
  (cmd_lcc_zarrB <- paste0(
    quotepath(py), ' ',
    quotepath(pyscriptB), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    pazarr, " ",
    format(smooth, scientific=F), ' ',
    format(tolerance, scientific=F), " ",
    format(ncores, scientific=F), " ",
    crs, " ",
    reOrderFile, " ",
    nodeidsFile, " ",
    format(sci, scientific=F), " ",
    format(eci, scientific=F), " "
    , ' 2>&1 ' #, logname
  ))

  if (cml){
    cat('\n\tCMD LCC zarr A:\n', cmd_lcc_zarrA)
    cat('\n')
  }

  intCMDA <- tryCatch(system(cmd_lcc_zarrA, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMDA)
  }

  if (cml){
    cat('\n\tCMD LCC zarr B:\n', cmd_lcc_zarrB)
    cat('\n')
  }

  intCMDB <- tryCatch(system(cmd_lcc_zarrB, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMDB)
  }

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
            log = paste0(" A: ", intCMDA, '\n',
                         " B: ", intCMDB, '\n') )
  #print('ANS LCC');print(ans)
  return( ans )
}


#' @title  Create cumulative resistance kernels
#' @description Cumulative resistant kernels
#' @param py Python executable location
#' @param pyscript Python script location
#' @param py Python executable location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
crk_py <- function(inshp, intif, outtif,
                   maxdist, shape, transform = 'no', volume,
                   ncores = as.numeric(Sys.getenv('COLA_NCORES')),
                   crs = 'None',
                   py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/crk.py'),
                   cml = TRUE, show.result = TRUE){

  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)
  # [5] kernel shape (linear, gaussian)
  # [6] kernel volume
  # [7] cores
  # [8] proj

  (cmd_crk <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(maxdist, scientific=F), ' ', # [4] distance threshold
    format(shape, scientific=F), ' ', # [5] kernel shape (linear, gaussian)
    transform, ' ', #
    format(volume, scientific=F), ' ', # [6] kernel volume
    format(ncores, scientific=F), ' ', # [7] cores
    crs
    , ' 2>&1 ' #, logname
  ) # [8] proj
  )
  (cmd_crk <- gsub(fixed = TRUE, '\\', '/', cmd_crk))
  if (cml){
    cat('\n\tCMD Kernel:\n',cmd_crk)
    cat('\n')
  }

  intCMD <- paste('', tryCatch(system(cmd_crk, intern = TRUE), error = function(e) e$message))
  #intCMD <- tryCatch(system(cmd_lcc, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  logname <- paste0(tools::file_path_sans_ext(outtif), '.metadata')
  metaFile <- c(inshp = inshp, intif = intif, outtif = outtif, maxdist = maxdist,
                shape = shape, transform = transform, volume = volume, ncores = ncores, crs = crs,
                log = paste0(intCMD, collapse = ' - '),
                done = ifelse(file.exists(outtif), 'yes', 'no'))
  write.table(metaFile, logname )

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
  # print('ANS CRK'); print(ans)
  return( ans )
}

#' @title  Create least cost corridors using parallel computing
#' @description Run CDPOP model
#' @param py Python executable location
#' @param pyscript Python script location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

crkJoblib_py <- function(
    inshp, intif, outtif,
    maxdist, shape,
    transform = 'no',
    volume,
    ncores = as.numeric(Sys.getenv('COLA_NCORES')),
    crs = 'None',
    maxram = 6,
    tempFolder = NULL,
    py = Sys.getenv("COLA_PYTHON_PATH"),
    pyscript = system.file(package = 'cola', 'python/crk_joblib.py'),
    cml = TRUE, show.result = TRUE){

  # "crk_joblib.py" "pts.shp inraster.tif out.tif 10000000 0 1000 6 None first.h5 second.h5 rmlimitinGB"
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (should be in meters*)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)
  # [7] number of cores
  # [8] projection if missing
  # [9] first h5 temp file
  # [10] second h5 temp file
  # [11] Max GB ram allowed

  if (is.null(tempFolder)){
    tempFolder <- tempdir()
  }
  tempH5 <- basename(tempfile())
  h5file <- paste0(tempFolder, '/', tempH5, '_A.h5')
  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')

  (cmd_crk <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(inshp), ' ',
    quotepath(intif), ' ',
    quotepath(outtif), ' ',
    format(maxdist, scientific=F), ' ',
    shape, ' ', transform, " ",
    format(volume, scientific=F), ' ', # kernel volume
    format(ncores, scientific=F), " ",
    crs, " ",
    h5file, " ",
    maxram
    , ' 2>&1 ' #, logname

  ))
  (cmd_crk <- gsub(fixed = TRUE, '\\', '/', cmd_crk))

  if (cml){
    cat('\n\tCMD CRK joblib:\n', cmd_crk)
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_crk, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  logname <- paste0(tools::file_path_sans_ext(outtif), '.metadata')
  metaFile <- c(inshp = inshp, intif = intif, outtif = outtif, maxdist = maxdist,
                shape = shape, transform = transform, volume = volume, ncores = ncores, crs = crs,
                maxram = maxram,
                log = paste0(intCMD, collapse = ' - '),
                done = ifelse(file.exists(outtif), 'yes', 'no'))
  write.table(metaFile, logname )


  tryCatch(file.remove(c(h5file, h5file2)), error = function(e) NULL)

  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
  #print('ANS CRK'); print(ans)
  return( ans )
}


#' @title  Runs prioritization
#' @description Run CDPOP model
#' @param py Python executable location
#' @param pyscript Python script location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
prio_py <- function(tif, incrk, inlcc,
                    maskedcsname = paste0(tempfile(), '.tif'),
                    outshppoint, outshppol, outshppatch,
                    outtifpatch, outtif,
                    threshold = 0.5, tolerance = 1000,
                    py = Sys.getenv("COLA_PYTHON_PATH"),
                    pyscript = system.file(package = 'cola', 'python/prioritize_core_conn.py'),
                    cml = TRUE, show.result = TRUE){

  # pri_py(py, incrk, inlcc, outshp, outtif, param5 = 0.5)
  # out_pri <- pri_py(py = py,
  #                    tif = rf$tif,
  #                    incrk = rv$crk ,
  #                    inlcc = rv$lcc,
  #                    maskedcsname = maskedcsname,
  #                    outshp = out_pri_shp,
  #                    outtif = out_pri_tif,
  #                    param5 = as.numeric(input$in_prio_5),
  #                    param6 = as.numeric(input$in_prio_6)) # 0.5


  # INPUTS
  # Original cost surface
  #ocsFile = sys.argv[1] # r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\size7.tif"

  # Path to resistant kernel file (output of crk function)
  #crkFile = sys.argv[2] # r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\crk100_test6.tif"

  # Path to least cost corridor file (output of lcc function)
  #lccFile = sys.argv[3] # r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\lcc100_test14.tif"

  # Output masked cost surface (could be a temp file)
  #maskedcsname = sys.argv[4] # r'C:/Users/pj276/Projects/UNICOR_arch/unicor/fix1_nodata-9999/size7_masked_cs.tif'

  # Output shapefile name
  #oshp = sys.argv[5] # r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\ptzCorr_points.shp"

  # Output raster name
  #ofile = sys.argv[6] #r"C:\Users\pj276\Projects\UNICOR_arch\unicor\fix1_nodata-9999\ptzCorr.tif"

  # Quantile threshold for identifying high value patches from resistant kernel surface
  # Should be between 0 and 1.
  #q = sys.argv[7] # 0.5
  #q = float(q)

  # Corridor tolerance. Should be the same tolerance used in the lcc script
  #corrTol = sys.argv[8]#1000

  # pyscript <- system.file(package = 'cola', 'python/prioritize_core_conn.py')

  # incrk <- '/data/temp/RY2024011519163805file176c0d742621ed/out_crk_IF2024011520212105file176c0d2f0a3994.tif'
  # inlcc <- '/data/temp/RY2024011519163805file176c0d742621ed/out_lcc_GG2024011519164505file176c0d5f4b6267.tif'
  # outshp <- '/data/temp/LA2024011522400305file1a50376cab710d//out_pri_EX2024011522411005file1a5037c8aef66.shp'
  # outtif <- '/data/temp/LA2024011522400305file1a50376cab710d//out_pri_EX2024011522411005file1a5037c8aef66.tif'

  logname <- paste0(tools::file_path_sans_ext(outtif), '.txt')


  (cmd_prio <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(tif), ' ',
    quotepath(incrk), ' ',
    quotepath(inlcc), ' ',
    quotepath(maskedcsname), ' ',
    quotepath(outshppoint), ' ',
    quotepath(outshppol), ' ',
    quotepath(outshppatch), ' ',
    quotepath(outtifpatch), ' ',
    quotepath(outtif), ' ',
    format(threshold, scientific=F), " ",
    format(tolerance, scientific=F)
    , ' 2>&1 ' #, logname
  ))

  if (cml){
    cat('\n\tCMD prio: \n')
    cat(cmd_prio <- gsub(fixed = TRUE, '\\', '/', cmd_prio))
    cat('\n')
  }


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE),
                     error = function(e) e$message)

  if(show.result){
    cat('\n\tlog CMD prio: \n')
    print(intCMD)
  }

  return(
    list(tif = ifelse(file.exists(outtif), outtif, ''),
         shp = ifelse(file.exists(outshppoint), outshppoint, ''),
         #log =  paste0(intCMD, ' -- ', read.delim(logname)) ))
         log = paste0("", intCMD) ) )
}


#' @title  Compare maps of cumulative resistance kernels
#' @description This tool compares the cumulative resistance kernels
#' @param py Python executable location
#' @param pyscript Python script location
#' @param py Python executable location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' crk_compare_py( )
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @export
crk_compare_py <- function(intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None',
                           py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/crk_compare.py'),
                           cml = TRUE, show.result = TRUE){

  # 'C:/Users/pj276/Scratch/scenario_testing/size7.tif
  # "C:/Users/pj276/Scratch/scenario_testing/size7_crk.tif,C:/Users/pj276/Scratch/scenario_testing/size7_s1_crk.tif,C:/Users/pj276/Scratch/scenario_testing/size7_s2_crk.tif"
  # C:/Users/pj276/Scratch/scenario_testing/size7_crk_abs_comp.png
  # C:/Users/pj276/Scratch/scenario_testing/size7_crk_percent_comp.png
  # C:/Users/pj276/Scratch/scenario_testing
  # C:/Users/pj276/Scratch/scenario_testing/mys_pas_ss.shp
  # WDPA_PID',
  # [1] Baseline raster path that contains original no data values
  # [2] string, coma separated raster to compare
  # [3] out csv with stats
  # [4] out PNG file
  # [5] out folder with pairwise comparisson
  # [6] shapefile to regionalize
  # [7] shapefile attribute/column name to aggregate
  logname <- paste0(tools::file_path_sans_ext(outfolder), '.txt')

  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

  (cmd_crk_comp <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(intif), ' ',
    quotepath(intifs), ' ',
    quotepath(outcsvabs), ' ',
    quotepath(outcsvrel), ' ',
    quotepath(outpngabs), ' ',
    quotepath(outpngrel), ' ',
    quotepath(outfolder), ' ',
    quotepath(inshp), ' ',
    quotepath(shpfield)
    , ' 2>&1 ' #, logname
  )
  )
  if (cml){
    cat('\n\tCMD Compare CRK: \n')
    cat(cmd_crk_comp <- gsub(fixed = TRUE, '\\', '/', cmd_crk_comp))
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_crk_comp, intern = TRUE),
                     error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(file.exists(outpngabs), outpngabs, ''),
               # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
               log = paste0("", intCMD) ) )
}

#' @title  Compare maps of least cost paths
#' @description This tool compares the least cost paths
#' @param py Python executable location
#' @param pyscript Python script location
#' @param py Python executable location
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path with CDPOP results
#' @examples
#' crk_compare_py( )
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @export
lcc_compare_py <- function(intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None',
                           py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/lcc_compare.py'),
                           cml = TRUE, show.result = TRUE){

  # 'C:/Users/pj276/Scratch/scenario_testing/size7.tif
  # "C:/size7_crk.tif,C:/size7_s1_crk.tif,C:/size7_s2_crk.tif"
  # C:/Users/pj276/Scratch/scenario_testing/size7_crk_abs_comp.png
  # C:/Users/pj276/Scratch/scenario_testing/size7_crk_percent_comp.png
  # C:/Users/pj276/Scratch/scenario_testing
  # C:/Users/pj276/Scratch/scenario_testing/mys_pas_ss.shp
  # WDPA_PID',
  # [1] Baseline raster path that contains original no data values
  # [2] string, coma separated raster to compare
  # [3] out csv with stats
  # [4] out PNG file
  # [5] out folder with pairwise comparisson
  # [6] shapefile to regionalize
  # [7] shapefile attribute/column name to aggregate
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  logname <- paste0(tools::file_path_sans_ext(outfolder), '.txt')


  (cmd_lcc_comp <- paste0(
    quotepath(py), ' ',
    quotepath(pyscript), ' ',
    quotepath(intif), ' ',
    quotepath(intifs), ' ',
    quotepath(outcsvabs), ' ',
    quotepath(outcsvrel), ' ',
    quotepath(outpngabs), ' ',
    quotepath(outpngrel), ' ',
    quotepath(outfolder), ' ',
    quotepath(inshp), ' ',
    quotepath(shpfield)
    , ' 2>&1 ' #, logname
  )
  )
  if (cml){
    cat('\n\tCMD Comp LCC: \n ')
    cat(cmd_lcc_comp <- gsub(fixed = TRUE, '\\', '/', cmd_lcc_comp))
    cat('\n')
  }

  intCMD <- tryCatch(system(cmd_lcc_comp, intern = TRUE), error = function(e) e$message)

  if(show.result){
    print(intCMD)
  }

  return( list(file = ifelse(file.exists(outpngabs), outpngabs, ''),
               #log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
               log = paste0("", intCMD) ) )
}


#' @title  Add values to  maps of least cost paths
#' @description Rasterize a polygon and sum it to an existing raster. Both layers need to be in the same projection.
#' @param polpath String. Location of the vector layer
#' @param burnval String. Value to burn. Can be a number or a attribute/column name
#' @param rastPath String. Location of the raster layer
#' @param att Logical. Should be 'all-the-touched' pixels be considered? Default TRUE
#' @param lineBuffW Number. How many pixels should be used as buffer width for line geometries?
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path of the resulting raster layer. Same as rastPath with the '_rasterized' suffix.
#' @examples
#' burnShp( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

burnShp <- function(polPath, burnval = 'val2burn',
                    colu = FALSE,
                    rastPath, rastCRS = NA,  att = FALSE, lineBuffW = 1){
  # test <- burnShp(polDraw, burnval, rastPath, rastCRS)
  # rastPath <- '/data/temp/XZ2024041911393405file9c152374e9a2/in_edit_DJ2024041911410705file9c15429f450d.tif'
  # rast <-terra::rast(rastPath)
  # burnval = -10
  # (load(file = '/data/tempR/draw.RData')) # polDraw

  #if( burnval != 0 & is.numeric(burnval) & !is.na(burnval) ){

  #(polPath <- gsub(x = rastPath, '.tif$', '_pol.shp'))
  (rasterizedPath <- gsub(x = rastPath, '.tif$', '_rasterized.tif'))

  file.copy(rastPath, rasterizedPath, overwrite = TRUE)

  #load(file = '/data/temp/2lines.RData') # polDraw #load(file = '/data/temp/4geom.RData') # polDraw
  #str(polDraw) #polDraw$type # FeatureCollection

  # rt <- terra::rast(rastPath)
  # rastRes <- res(rt)

  if( is.na(rastCRS)){
    ## rastPath <- '/home/shiny/connecting-landscapes/docs/HS_size5_nd_squared.tif'
    #gi <- rgdal::GDALinfo(rastPath)
    #gi2 <- sf::gdal_crs(rastPath)
    #gi <- strsplit(x = gdalUtilities::gdalinfo(rastPath), '\n')
    #prj <- attr(x = gi, "projection")

    # prj <- terra::crs(rt, proj = TRUE)
    # (rastCRS <- st_crs(rt))

    # if(!is.na(prj)){
    #   pol2save@proj4string@projargs <- prj
    # }
    #ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:NEW_EPSG_NUMBER -s_srs EPSG:OLD_EPSG_NUMBER output.shp input.shp
    #EPSG:4326
  }

  # polPath <- '/data/tempR//colaQMX2024101612582905//Aproj_ADB_FeasibilityAlignment.shp'
  # rasterizedPath <- '/data/tempR//colaQMX2024101612582905//in_edit_fixed_RKS2024101613001505_rasterized.tif'
  # att = TRUE; burnval = '200'

  if( colu ) {
    print(' Add vals -- column')
    err <- tryCatch(
      gdalUtilities::gdal_rasterize(
        src_datasource = polPath,
        at = att,
        dst_filename = rasterizedPath,
        add = TRUE,
        a = burnval) #as.numeric(burnval)
      , error = function(e) e$message)

  } else {
    print(' Add vals -- value')
    err <- tryCatch(
      gdalUtilities::gdal_rasterize(
        src_datasource = polPath,
        at = att,
        dst_filename = rasterizedPath,
        add = TRUE,
        burn = burnval) #as.numeric(burnval)
      , error = function(e) e$message)
  }



  #file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)
  # rasteri <- gdalUtils::gdal_rasterize(src_datasource = polPath, at = T,
  #                                      dst_filename = rastPath,
  #                                      add = TRUE, a = 'val2burn')
  # plot(raster(rastPath))
  # plot(raster(rasterizedPath), main = 'Rasterized')
  # plot(sf::read_sf(polPath), add = TRUE)
  # file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)

  return(rasterizedPath)
  # } else {
  #   return(NA)
  # }
}


#' @title  Add values to  maps of least cost paths
#' @description Rasterize a polygon and sum it to an existing raster. Both layers need to be in the same projection.
#' @param polpath String. Location of the vector layer
#' @param burnval String. Value to burn. Can be a number or a attribute/column name
#' @param rastPath String. Location of the raster layer
#' @param att Logical. Should be 'all-the-touched' pixels be considered? Default TRUE
#' @param lineBuffW Number. How many pixels should be used as buffer width for line geometries?
#' @param cml Logical. Print the back-end command line? Default TRUE
#' @param show.result Logical. Print the command line result? Default TRUE
#' @return Path of the resulting raster layer. Same as rastPath with the '_replaced ' suffix.
#' @examples
#' replacePixels( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

replacePixels <- function(polPath, burnval = 'val2burn', rastPath, colu = FALSE,
                          att = FALSE, rastCRS = NA, gdal = TRUE){

  #test <- replacePixels(polDraw, burnval, rastPath, rastCRS)

  # polPath <- '/data/tempR/colaBMJ2024101517341605/proj_ADB_FeasibilityAlignment.shp'
  # rastPath <- '/data/tempR/colaBMJ2024101517341605/in_edit_fixed_TKG2024101517383805.tif'

  #if( burnval != 0 & is.numeric(burnval) & !is.na(burnval) ){
  ## Polygon to write
  #(polPath <- gsub(x = rastPath, '.tif$', '_pol.shp'))
  ## Raster with new features
  (rasterizedPath <- gsub(x = rastPath, '.tif$', '_rasterized2replace.tif'))
  ## Raster to create
  (replacedPath <- gsub(x = rastPath, '.tif$', '_replaced.tif'))


  rtp <- terra::rast(rastPath)
  # (rastRes <- res(rt))

  # rastPath <- '/data/tempR//colaZTL2024101522171205//in_edit_fixed_ILK2024101522172305.tif'
  gi <- gdalUtilities::gdalinfo(rastPath, quiet = TRUE)
  rastRes0 <- grep('Pixel Size', strsplit(gi, '\n')[[1]], value = TRUE)
  base::options(scipen = 999)
  (rastRes <- (strsplit(x = gsub(pattern = '.+\\(|\\)|-', '', rastRes0), ',')[[1]]))
  ts0 <- grep('Size is', strsplit(gi, '\n')[[1]], value = TRUE)
  (ts <- (strsplit(x = gsub(pattern = 'Size is | ', '', ts0), ',')[[1]]))
  base::options(scipen=0, digits=7)

  if( is.na(rastCRS)){
    #prj <- terra::crs(rt, proj = TRUE)
    #(rastCRS <- st_crs(rt))
  }

  # pol2Rast <- draws2Features(polDraw, distLineBuf = min(rastRes), rastCRS = rastCRS)
  # # pol2Rastx <- st_sf(data.frame(a = 1:length(pol2Rast), pol2Rast))
  # pol2Rastx <- st_as_sf( pol2Rast)
  # # plot(pol2Rastx, add = TRUE, border = 'blue', col = NA)
  #
  # pol2Rastx$val2burn <- as.numeric(burnval)
  # sf::st_write( obj = pol2Rastx, dsn = dirname(polPath),
  #               layer = tools::file_path_sans_ext(basename(polPath)),
  #               driver = 'ESRI Shapefile',
  #               append = FALSE,
  #               overwrite_layer = TRUE)


  ## Rasterize polygons
  rastExtent <- as.character(as.vector(terra::ext(rtp))[c('xmin', 'ymin', 'xmax', 'ymax')])

  # cat(' --- GdalRasterize: \n ++ te:', rastExtent)
  # cat('\n ++ tr ', rastRes)
  # cat('\n ++ ts ', ts)
  # cat('\n ++ polPath ', polPath)
  # cat('\n ++ rasterizedPath ', rasterizedPath)
  # cat('\n ++ burnval ', burnval)
  if (colu){
    ## Use a column/attribute
    err <- tryCatch(
      gdalUtilities::gdal_rasterize(
        te = rastExtent, # -te <xmin> <ymin> <xmax> <ymax>
        #tr = rastRes,
        ts = ts, # c(terra::ncol(rt), terra::nrow(rt)), # <width> <height>
        src_datasource = polPath,
        at = att,
        dst_filename = rasterizedPath,
        add = FALSE,
        a = burnval) #as.numeric(burnval)
      , error = function(e) {
        print('Error rast pol for replacing');
        print(e);
        return(e)
      })

  } else{
    ## Use a single value
    err <- tryCatch(
      gdalUtilities::gdal_rasterize(
        te = rastExtent, # -te <xmin> <ymin> <xmax> <ymax>
        #tr = rastRes,
        ts = ts, # c(terra::ncol(rt), terra::nrow(rt)), # <width> <height>
        src_datasource = polPath,
        at = att,
        dst_filename = rasterizedPath,
        add = FALSE,
        burn = burnval) #as.numeric(burnval)
      , error = function(e) {print('Error rast pol for replacing'); print(e);e})
  }

  if (!file.exists(rasterizedPath)){
    cat (' Error at rasterizing')
  } else {

    # rasterizedPath <- 'C:/cola/colaGLO2024121820390005/out_crk_EGL2024121820455305.tif'
    g2 <- gdalUtilities::gdalinfo(rasterizedPath, quiet = TRUE)
    cat(g2)

    #rft <- rast(rasterizedPath); plot(rft)

    if (gdal){
      cat(' --- Rasterizing with gdal_calc.py')
      # Use gdal calc in linux
      (cmdCalc <- paste0('gdal_calc.py --overwrite ',
                         ' -A ', rastPath, # original
                         ' -B ', rasterizedPath,
                         ' --calc "((B == 0 ) * A ) + ((B != 0 ) * B )"',
                         ' --outfile ', replacedPath))

      runCMD <- tryCatch(system(cmdCalc, intern = TRUE), error = function(e) NA)
    } else {
      cat(' --- Rasterizing with terra')

      ## GDAL calc
      ## works on win> C:\Users\gonza>C:\Users\gonza\AppData\Local\r-miniconda\envs\cola\python.exe C:\Users\gonza\AppData\Local\r-miniconda\envs\cola\Scripts\gdal_calc.py --help


      # gdalCalc <- file.path(dirname(Sys.getenv('COLA_PYTHON_PATH')), 'Scripts', 'gdal_calc.py')
      # if (file.exists(gdalCalc)){
      #
      #   testFile <- gsub('\\\\', '/', paste0(tempfile(), '.tif'))
      #   (gc_cmd <- paste0(gdalCalc, ' -A ',
      #                     system.file(package = 'cola', 'sampledata/sampleTif.tif'),
      #                     ' --outfile=', testFile, ' --calc="(A*2)"' ))
      #   cat(gc_cmd)
      #
      #   out_gc <- (system(gc_cmd, intern = TRUE, show.output.on.console = TRUE))
      #
      #   (a <- system2(gc_cmd, timeout = 10, wait = TRUE))
      #   (b <- system(gc_cmd, timeout = 10, wait = TRUE, intern = TRUE))
      # }

      ## Run analysis on Windows
      # Use terra::rast

      # rastPath  <- '/data/tempR//colaRFY2024100813020705/out_surface_GZO2024100813024405_rasterized2replace.tif'
      # rasterizedPath <- '/data/tempR//colaRFY2024100813020705/out_surface_GZO2024100813024405_replaced.tif'
      A <- terra::rast(rastPath)
      B <- terra::rast(rasterizedPath)
      # plot(A)
      # plot(B)

      newRast <- ((B == 0 ) * A ) + ((B != 0 ) * B )
      # plot(newRast)

      terra::writeRaster(newRast, filename = replacedPath)
    }
  }

  # plot(rast(rastPath))
  # plot(rast(rasterizedPath), main = 'Rasterized')
  # plot(sf::read_sf(polPath), add = TRUE, col = NA, border = 'red')
  # plot(rast(replacedPath), main = 'Replaced')
  # file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)

  return(replacedPath)

  #} else { return(NA) }
}
