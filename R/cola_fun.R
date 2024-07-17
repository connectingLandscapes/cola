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


#' @title  Run CDPOP model
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

runCDPOP <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                     inputvars = NULL,
                     agevars = NULL,
                     cdmat = NULL,
                     xy = NULL,
                     tempFolder){
  #xyfilename  NO .csv required
  #agefilename .csv required
  #matecdmat	cdmats/EDcdmatrix16
  #dispcdmat	cdmats/EDcdmatrix16



  # file.copy('invars.csv', '/home/user/cola/inst/examples/invars.csv')
  # file.copy('age.csv', '/home/user/cola/inst/examples/agevars.csv')
  if (is.null(inputvars)){
    inputvars <- system.file(package = 'cola', 'sampledata/invars.csv')
    file.copy(from = inputvars, to = file.path(tempFolder, 'invars'), overwrite = TRUE)
  }

  if (is.null(agevars)){
    agevars <- system.file(package = 'cola', 'sampledata/age.csv')
  }

  if (is.null(cdmat)){
    stop(' CDmatrix required')
  }

  if (is.null(xy)){
    stop(' XY required')
  }

  py <- Sys.getenv("COLA_PYTHON_PATH")
  # pyscript <- '/home/shiny/connecting-landscapes/lib/CDPOP/src/CDPOP.py'
  pyscript <- '/mnt/c/temp/cdpopx/connecting-landscapes/lib/CDPOP/src/CDPOP.py'

  datapath <- tempFolder # datapath = tempFolder
  timeMarkCDPOP <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE))
  (cdpopPath <- paste0('cdpopout_', timeMarkCDPOP, '__'))
  #cdpopPath <- 'cdpopout'

  # file.copy('out_cdmatrix_KQI.csv', 'cdmat.csv', overwrite = TRUE)
  # tempFolder <- '/tmp/RtmpiD9uhw/colaMJJ2024070817062505'; setwd(tempFolder); list.files()
  # tempFolder <- '/mnt/c/temp/tempCola'; setwd(tempFolder); list.files()
  # inputvars <- read.csv('invars.csv')
  # agevars <- read.csv('age.csv')
  # xyorig <- read.csv('xy.csv'); dim(xy)
  # unique(xy$DisperseCDist)
  # cdmat <- read.csv('cdmat.csv', header = FALSE); dim(cdmat)
  #cdmat2 <- read.csv('out_cdmatrix_KQI2024070817135705.csv', header = FALSE); dim(cdmat2)
  #cdmat2[1:5, 1:5]
  # vars <- paste0('invars.csv') # Only file name read.csv(vars)

  # setwd(tempFolder)
  # file.copy('xyA.csv', 'xy.csv')
  # templ_invars <- ''
  # templ_xy <- ''
  # templ_agevars <- ''

  # python CDPOP.py %userprofile%\dockerdata\CDPOP\data inputvars.csv outAnac1
  (cmd <- paste0(py, ' ', pyscript, ' ', datapath, ' ', vars, ' ', cdpopPath))
  cat(' CMD:\n', cmd, '\n')

  intCMD <- tryCatch(system(cmd, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)

  newFiles0 <- list.files(path = tempFolder, recursive = TRUE)
  newFiles <- grep(pattern = cdpopPath, x = newFiles0, value = TRUE)
  return(list(newFiles = newFiles, cdpopPath = cdpopPath))

  # outfiles:
  #CDPOP/data/out11684988478/batchrun0mcrun0/grid0.csv [0, 1, ..]
  #CDPOP/data/out11684988478/batchrun0mcrun0/output.csv
  #CDPOP/data/out11684988478/batchrun1mcrun0/XY0.csv
  #CDPOP/data/out21684988758/cdpop.log
}


#' @title  Shapefile to CDPOP xy
#' @description Convert Shapefile into xy file
#' @param shapefile String. Path to the file
#' @examples
#' shp2xy( shapefile = 'shapefilepathhere.shp',outxy = 'out.xy', tempDir = 'temFolder')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

shp2xy <- function(shapefile, outxy, tempDir){

  # The n-(x,y) grid location values. This is a comma delimited file with 5 column headings:
  # (Subpopulation)- a unique identifier for  each individual corresponding to a unique subpopulation
  # (XCOORD)-x-coordinate location, (YCOORD)-y-coordinate location (YCOORD)
  # (ID)-a string label identifier, and
  # (sex)-an initial sex assignment (use 0/1 or F/M).
  # See xyED16.csv for an example xyfilename. The column order is necessary with a header file.


  # tempDir = "/tmp/RtmpiD9uhw/colaMJJ2024070817062505"
  # tempDir = "/mnt/c/temp/tempCola/"
  # shapefile = '/mnt/c/temp/tempCola/out_simpts_XYA2024070817112705.shp'
  # outxy = '/mnt/c/temp/tempCola/out_simpts_XYA2024070817112705.csv'
  # df <- foreign::read.dbf(gsub('\\..+', '.dbf', shapefile))
  # xyorig <- system.file(package = 'cola', 'sampledata/xy.csv')
  xy <- terra::vect(shapefile)

  xy$ID <- xy$Subpopulation <- 1:nrow(xy)
  if (! 'X' %in% names(xy) ){
    xy_coord <- data.frame(geom(xy))
    xy$XCOORD <- xy_coord$x
    xy$YCOORD <- xy_coord$y
  }
  xy$Subpop_mortperc <- 0
  xy$sex <- sample(x = c(1, 0), size = nrow(xy), replace = TRUE)
  xynew <- as.data.frame(xy)
  xynew <- xynew[, c('Subpopulation', 'X', 'Y', 'Subpop_mortperc', 'ID', 'sex')]


  # ## Write CSV
  write.csv(x = xynew, file = outxy, row.names = FALSE, quote = FALSE)

  ## Write CSV as xy.csv
  xynew[1:3, 1:3]
  write.csv(x = xynew, file = file.path(tempDir, 'xy.csv'), row.names = FALSE, quote = FALSE)


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
#' @return String with the NA value. No converted to number to avoid characters
#' @examples
#' input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' guessNoData(input_tif)
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
guessNoData <- function(path){
  # path = raster path
  ans <- NA
  if (require(gdalUtilities) & file.exists(path)){
    gi <- strsplit(gdalinfo(path), '\n')[[1]]
    ndv <- grep('NoData ', gi, value = TRUE)
    if( any(length(ndv)) ) {
      (ans <- gsub('.+\\=', '', ndv))
    }
  }
  return(ans)
}

#' @title  Adapt file path. Change backslash to slash
#' @description Fix paths so internal console recognize paths. Change backslash to backslash
#' @param path File location
#' @return File location String using slash as separator
#' @examples
#' newPath <- adaptFilePath('\temp\path\to\file\here.tif')
#' newPath
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
adaptFilePath <- function(path){
  # path = "C:\\Users\\ig299\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/python/s2res.py C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/sampledata/sampleTif.tif C:\\Users\\ig299\\AppData\\Local\\Temp\\RtmpwrUyVu/VM2024041715525605file51c6028258//out_surface_JQ2024041715525705file51c260d2c4b.tif 0.06788435 0.9989325 100 1 -9999 None"
  return ( gsub(fixed = TRUE, '\\', '/', path) )
}

#' @title  Suitability surface to resistance surface
#' @description Transforms suitability to resistance surface
#' @param py Python location or executable. The string used in R command line to activate `cola`. The default versio should point to a conda environment. Might change among computers versions
#' @param pyscript Python script location
#' @param intif String. File path to the input raster.
#' @param outtif String. File path of the output surface resistance
#' @param param3 Numeric. The lower value on the input raster to cut off. Pixels with values under the given number will be ignored. In the front end the values automatically derived from the input file.
#' @param param4 Numeric. The upper value on the input raster to cut off. Pixels with values under the given number will be ignored. In the front end the values automatically derived from the input file.
#' @param param5 Numeric. This is the maximum resistance value after transformation from suitability. Default value is 100.  Minimum value is set to 1.
#' @param param6 Numeric. A statistical parameters that defines the transformation pattern between the input and output. The shape value determines the relationship between suitability and resistance. For a linear relationship, use a value close to 0, such as 0.01. Positive values result in a greater increase in resistance as suitability declines. This is appropriate for animals that are more sensitive to the matrix in between habitat. Negative values result in a lesser increase in resistance as suitability declines. This is appropriate for animals that are less sensitive to the matrix between habitat. The more positive or more negative, the greater the effect on the shape of the relationship. Values generally range between +10 and -10, 0 is not allowed.
#' @param param7 Numeric. The no data value of the input file. For GeoTiffs, this is automatically determined. For text based files, this must be input by the user. Default value is ‘None’.
#' @param param8 Projection information in the case the input raster [1] has no spatial projection. For GeoTiffs, this is automatically determined. For text based files, this must be input by the user. Provide it as EPSG or ESRI string e.g. "ESRI:102028". Default value is ‘None’.
#' @return Creates a raster layer with a minimum value of 1 and maximum value given the parameter 5. The internal R object is a list of two slots. The first one contains the path of the created raster, if any, and the second slot includes any function message or log, if any.
#' @examples
#' hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' srp30 <- s2res_py(intif = hs, outtif = 'hs_p3_0.tif', param3 = 0, param4 = 1, param5 = 10,  param6 = 1,  param7 = NULL, param8 = 'None')
#' hs_rast <- terra::rast(hs)
#' plot(hs_rast, main = 'Habitat suitability')
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
s2res_py <- function(intif, outtif,
                     param3, param4, param5, param6,
                     param7 = NULL, param8 = 'None',
                     py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/s2res.py')){
  # param3 = 0
  # param4 =  100
  # param5 = 100
  # param6 = 1
  # param7 = -9999
  # param8 = 'None'

  ## Guess NA value if not provided
  if(is.null(param7)){
    param7 <- guessNoData(intif)
  }

  ## Create CMD
  (cmd_s2res <- paste0(py, ' ', pyscript, ' ',
                       intif, ' ', outtif, ' ',
                       format(param3, scientific=F), ' ',
                       format(param4, scientific=F), ' ',
                       format(param5, scientific=F), ' ',
                       format(param6, scientific=F), ' ',
                       format(param7, scientific=F), ' ',
                       param8))
  cat('\n\n\tCMD:')

  cat(cmd_s2res <- gsub(fixed = TRUE, '\\', '/', cmd_s2res))

  intCMD <- tryCatch(system( cmd_s2res ,
                             intern = TRUE, ignore.stdout = TRUE),
                     error = function(e) e$message)
  return( list(file = ifelse(file.exists(outtif), outtif, NA),
               log =  intCMD) )
}

#' @title  Create random points
#' @description Create random points
#' @param rvect Python location
#' @param npts Python location
#' @param rmin Python location
#' @param rmax Python location
#' @return Path with CDPOP results
#' @examples
#' hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
#' pts30 <- s2res_py(intif = hs, outtif = 'hs_p3_0.tif', param3 = 0, param4 = 1, param5 = 10,  param6 = 1,  param7 = NULL, param8 = 'None')
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


#' @title  Transforms suitability to resistance surface
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
points_py <- function(intif, outshp,
                      param3, param4, param5, param6 = 'None',
                      py = Sys.getenv("COLA_PYTHON_PATH"),
                      pyscript = system.file(package = 'cola', 'python/create_source_points.py')){
  # param3 = 2
  # param4 =  95
  # param5 = 50
  # pyscript <- system.file(package = 'cola', 'python/create_source_points.py'
  (cmd_pts <- paste0(py, ' ', pyscript, ' ', intif, ' ', outshp, ' ',
                     format(param3, scientific=F), ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     param6))
  cat('\n\n\tCMD:')

  print(cmd_pts <- gsub(fixed = TRUE, '\\', '/', cmd_pts))

  intCMD <- tryCatch(system(cmd_pts, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  return( list(file = ifelse(file.exists(outshp), outshp, NA),
               log =  intCMD) )

}

#' @title  Creates CDmatrix for CDPOP model
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

cdmat_py <- function(inshp, intif, outcsv,
                     param4,
                     param5 = Sys.getenv("COLA_NCORES"),
                     param6 = 'None',
                     py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/create_cdmat.py')){
  # param4 = 100000
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
  # param4 = 100000; param5 = 1; param6 = 'None'

  # pyscript <- system.file(package = 'cola', 'python/create_cdmat.py')
  (cmd_cdmat <- paste0(py, ' ', pyscript, ' ', inshp, ' ', intif, ' ', outcsv,
                       ' ', param4, ' ', param5, ' ', param6))

  cat('\n\n\tCMD:')
  print(cmd_cdmat <- gsub(fixed = TRUE, '\\', '/', cmd_cdmat))

  intCMD <- tryCatch(system(cmd_cdmat, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  #checkcsv <- read.csv(outcsv); which(is.numeric(checkcsv)) ; summary(checkcsv); sum(checkcsv, )


  return( list(file = ifelse(file.exists(outcsv), outcsv, NA),
               log =  intCMD) )
}


#' @title  Create least cost corridors
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lcc_py <- function(inshp, intif, outtif,
                   param4, param5, param6,
                   param7 = as.numeric(Sys.getenv('COLA_NCORES')), param8 = 'None',
                   py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/lcc.py')){
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] Max. dispersal distance (meters)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)

  # if(is.null(param7)){
  #   param7 <- guessNoData(intif)
  # }

  (cmd_lcc <- paste0(py, ' ', pyscript, ' ',
                     inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     format(param6, scientific=F), " ",
                     format(param7, scientific=F), " ",
                     param8))
  cat('\n\n\tCMD:')
  print(cmd_lcc <- gsub(fixed = TRUE, '\\', '/', cmd_lcc))


  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  return( list(file = ifelse(file.exists(outtif), outtif, NA),
               log =  intCMD) )
}


#' @title  Create least cost corridors for heavy rasters
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

lccHeavy_py <- function(inshp, intif, outtif,
                        param4, param5, param6,
                        param7 = as.numeric(Sys.getenv('COLA_NCORES')),
                        param8 = 'None', tempFolder = rootPath,
                        py = Sys.getenv("COLA_PYTHON_PATH"),
                        pyscript = system.file(package = 'cola', 'python/lcc_heavy.py')){

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

  tempH5 <- sessionIDgen()
  h5file1 <- paste0(rootPath, '/', tempFolder, '/', tempH5, '_A.h5')
  h5file2 <- paste0(rootPath, '/', tempFolder, '/', tempH5, '_B.h5')


  (cmd_lcc <- paste0(py, ' ', pyscript, ' ', inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     format(param6, scientific=F), " ",
                     format(param7, scientific=F), " ",
                     param8, " ",
                     h5file1, " ",
                     h5file2, " ",
                     '50'
  ))
  cat('\n\n\tCMD:')
  print(cmd_lcc <- gsub(fixed = TRUE, '\\', '/', cmd_lcc))

  file.remove(c(h5file1, h5file2))

  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  return( list(file = ifelse(file.exists(outtif), outtif, NA),
               log =  intCMD) )
}

#' @title  Create cumulative resistance kernels
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#'
crk_py <- function(inshp, intif, outtif,
                   param4, param5, param6,
                   param7 = as.numeric(Sys.getenv('COLA_NCORES')), param8 = 'None',
                   py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/crk.py')){

  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)
  # [5] kernel shape (linear, gaussian)
  # [6] kernel volume
  # [7] cores
  # [8] proj

  (cmd_crk <- paste0(py, ' ', pyscript, ' ', inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ', # [4] distance threshold
                     format(param5, scientific=F), ' ', # [5] kernel shape (linear, gaussian)
                     format(param6, scientific=F), ' ', # [6] kernel volume
                     format(param7, scientific=F), ' ', # [7] cores
                     param8) # [8] proj
  )
  cat('\n\n\tCMD:')
  print(cmd_crk <- gsub(fixed = TRUE, '\\', '/', cmd_crk))

  intCMD <- tryCatch(system(cmd_crk, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  return( list(file = ifelse(file.exists(outtif), outtif, NA),
               log =  intCMD) )
}


#' @title  Runs prioritization
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#'
pri_py <- function(tif, incrk, inlcc,
                   maskedcsname = paste0(tempfile(), '.tif'),
                   outshppoint, outshppol, outshppatch,
                   outtifpatch, outtif,
                   param7 = 0.5, param8 = 1000,
                   py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/prioritize_core_conn.py')){

  # pri_py(py, incrk, inlcc, outshp, outif, param5 = 0.5)
  # out_pri <- pri_py(py = py,
  #                    tif = rf$tif,
  #                    incrk = rv$crk ,
  #                    inlcc = rv$lcc,
  #                    maskedcsname = maskedcsname,
  #                    outshp = out_pri_shp,
  #                    outif = out_pri_tif,
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

  (cmd_prio <- paste0(py, ' ',
                      pyscript, ' ',
                      tif, ' ',
                      incrk, ' ', inlcc, ' ',
                      maskedcsname, ' ',
                      outshppoint, ' ',
                      outshppol, ' ',
                      outshppatch, ' ',
                      outtifpatch, ' ',
                      outtif, ' ',
                      format(param7, scientific=F), " ",
                      format(param8, scientific=F)))

  print(' CMD prio:')
  print(cmd_prio <- gsub(fixed = TRUE, '\\', '/', cmd_prio))


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  print(intCMD)
  return(
    list(tif = ifelse(file.exists(outtif), outtif, NA),
         shp = ifelse(file.exists(outshppoint), outshppoint, NA),
         log = intCMD)
  )
}

#' @title  Runs CDOPOP2
#' @description Run CDPOP model
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
cdpop_py <- function(outshp, outtif,
                     param5 = 0.5, param6 = 1000,
                     tif, incrk, inlcc,
                     maskedcsname = paste0(tempfile(), '.tif'),
                     py = Sys.getenv("COLA_PYTHON_PATH")){


  #(base) C:\Users\Admin\dockerdata\CDPOP\pyscript>
  #  python CDPOP.py %userprofile%\dockerdata\CDPOP\data invars.csv outLinux
  ### Should be Y in the column output_unicor


  # pri_py(py, incrk, inlcc, outshp, outif, param5 = 0.5)
  # out_pri <- pri_py(py = py,
  #                    tif = rf$tif,
  #                    incrk = rv$crk ,
  #                    inlcc = rv$lcc,
  #                    maskedcsname = maskedcsname,
  #                    outshp = out_pri_shp,
  #                    outif = out_pri_tif,
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

  pyscript <- system.file(package = 'cola', 'python/prioritize_core_conn.py')


  (cmd_prio <- paste0(py, ' ',
                      pyscript, ' ',
                      tif, ' ',
                      incrk, ' ', inlcc, ' ',
                      maskedcsname, ' ',
                      outshp, ' ', outtif, ' ',
                      format(param5, scientific=F), " ",
                      format(param6, scientific=F)))

  print(cmd_prio <- gsub(fixed = TRUE, '\\', '/', cmd_prio))


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  print(intCMD)
  return(list(tif = ifelse(file.exists(outtif), outtif, NA),
              shp = ifelse(file.exists(outshp), outshp, NA),
              log = intCMD) )
}


#' @title  Compare maps of cumulative resistance kernels
#' @description This tool compares the cumulative resistance kernels
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' crk_compare_py( )
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @author Ivan Gonzalez <ig299@@nau.edu>
#'
crk_compare_py <- function(intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None',
                           py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/crk_compare.py')){

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
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

  (cmd_crk_comp <- paste0(py, ' ', pyscript, ' ',
                          intif, ' ', intifs, ' ',
                          outcsvabs, ' ',
                          outcsvrel, ' ',
                          outpngabs, ' ',
                          outpngrel, ' ',
                          outfolder, ' ',
                          inshp, ' ', shpfield)
  )
  cat('\n\n\tCMD:')
  print(cmd_crk_comp <- gsub(fixed = TRUE, '\\', '/', cmd_crk_comp))

  intCMD <- tryCatch(system(cmd_crk_comp, intern = TRUE, ignore.stdout = TRUE),
                     error = function(e) e$message)
  return( list(file = ifelse(file.exists(outpngabs), outpngabs, NA),
               log =  intCMD) )
}

#' @title  Compare maps of least cost paths
#' @description This tool compares the least cost paths
#' @param py Python location
#' @param py Python location
#' @return Path with CDPOP results
#' @examples
#' crk_compare_py( )
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @author Ivan Gonzalez <ig299@@nau.edu>
#'
lcc_compare_py <- function(intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None',
                           py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/lcc_compare.py')){

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

  (cmd_lcc_comp <- paste0(py, ' ', pyscript, ' ',
                          intif, ' ', intifs, ' ',
                          outcsvabs, ' ',
                          outcsvrel, ' ',
                          outpngabs, ' ',
                          outpngrel, ' ',
                          outfolder, ' ',
                          inshp, ' ', shpfield)
  )
  cat('\n\n\tCMD:')
  print(cmd_lcc_comp <- gsub(fixed = TRUE, '\\', '/', cmd_lcc_comp))

  intCMD <- tryCatch(system(cmd_lcc_comp, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  return( list(file = ifelse(file.exists(outpngabs), outpngabs, NA),
               log =  intCMD) )
}

