#' @title Launch \emph{COLA} decision support system (DSS) dashboard
#' @description This function runs the \emph{cola} dashboard
#' @param launch.browser Run this app on a new window? Default `TRUE`
#' @note Please see the official website (\url{https://wallaceecomod.github.io/})
#' for more details. If you have questions about the application,
#' please participate in the \href{https://groups.google.com/forum/#!forum/wallaceecomod}{Google Group},
#' or email the team directly: \email{wallaceEcoMod@@gmail.com}.
#'
#' @examples
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

runCDPOP <- function(py = Sys.getenv("COLA_PYTHON_PATH"), datapath = tempFolder){
  # outfiles:
  #CDPOP/data/out11684988478/batchrun0mcrun0/grid0.csv [0, 1, ..]
  #CDPOP/data/out11684988478/batchrun0mcrun0/output.csv
  #CDPOP/data/out11684988478/batchrun1mcrun0/XY0.csv
  #CDPOP/data/out21684988758/cdpop.log

  #xyfilename  requires NO .csv
  #agefilename requires .csv
  #matecdmat	cdmats/EDcdmatrix16
  #dispcdmat	cdmats/EDcdmatrix16

  # python CDPOP.py %userprofile%\dockerdata\CDPOP\data inputvars.csv outAnac1
  pyscript <- '/home/shiny/connecting-landscapes/lib/CDPOP/src/CDPOP.py'
  datapath <- tempFolder # datapath = tempFolder
  vars <- paste0('invars.csv') # Only file name
  timeMarkCDPOP <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE))
  cdpopPath <- paste0('cdpopout_', timeMarkCDPOP, '__')
  cdpopPath <- 'cdpopout'
  (cmd <- paste0(py, ' ', pyscript, ' ', datapath, ' ', vars, ' ', cdpopPath))

  # setwd(tempFolder)
  #file.copy('inputvars.csv', 'in.csv')
  #file.copy('xyA.csv', 'xy.csv')
  #intCMD <- tryCatch(system(cmd, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)

  newFiles <- list.files(path = tempFolder, recursive = TRUE)
  newFiles <- grep(pattern = cdpopPath, newFiles, value = TRUE)
  return(list(newFiles = newFiles, cdpopPath = cdpopPath))
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
    gi <- strsplit(gdalinfo(intif), '\n')[[1]]
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

#' @title  Transforms suitability to resistance surface
#' @description Run CDPOP model
#' @param py Python location or executable. The string used in R command line to activate `cola`
#' conda environment. Might change among computers versions
#' @param pyscript Python script location
#' @param intif Path to the Suitability surface raster layer (TIF format)
#' @param outtif Path to the resulting Surface Resistance (SR) raster layer (TIF format)
#' @param param3 Suitability grid min value. Numeric value to cutoff the original layer.
#' @param param4 Suitability grid max value. Numeric value to cutoff the original layer.
#' @param param5 Maximum resistance value to have the new raster. Min value is set to 1.
#' @param param6 Shape Parameter. Between 1 and X
#' @param param7 No data value of suitability raster. Default NULL. To be guessed from metadata if not provided
#' @param param8 CRS if using ASCII|RSG or other file without projection info. Provide as EPSG or ESRI string e.g. "ESRI:102028"
#' @return Path with CDPOP results
#' @examples
#' runCDPOP( )
#' @author Ivan Gonzalez <ig299@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export
s2res_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/s2res.py'),
                     intif, outtif,
                     param3, param4, param5, param6,
                     param7 = NULL, param8 = 'None'){
  # param3 = 0
  # param4 =  100
  # param5 = 100
  # param6 = 1
  # param7 = -9999
  # param8 = 'None

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
points_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                      pyscript = system.file(package = 'cola', 'python/create_source_points.py'),
                      intif, outshp,
                      param3, param4, param5, param6 = 'None'){
  # param3 = 2
  # param4 =  95
  # param5 = 50
  #pyscript <- system.file(package = 'cola', 'python/create_source_points.py'
  datapath <- tempFolder # datapath = tempFolder
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

cdmat_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/create_cdmat.py'),
                     inshp, intif, outcsv,
                     param3, param4,
                     param5 = 1, param6 = 'None'){
  # param3 = 25000
  # create_cdmat.py
  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)

  # pyscript <- system.file(package = 'cola', 'python/create_cdmat.py')
  (cmd_cdmat <- paste0(py, ' ', pyscript, ' ', inshp, ' ', intif, ' ', outcsv,
                       ' ', format(param3, scientific=F),
                       ' ', param4, ' ', param5))

  cat('\n\n\tCMD:')

  print(cmd_cdmat <- gsub(fixed = TRUE, '\\', '/', cmd_cdmat))

  intCMD <- tryCatch(system(cmd_cdmat, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
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

lcc_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/lcc.py'),
                   inshp, intif, outtif,
                   param4, param5, param6,
                   param7 = as.numeric(Sys.getenv('COLA_NCORES')), param8 = 'None'){
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (should be in meters*)
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

lccHeavy_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                     pyscript = system.file(package = 'cola', 'python/lcc_heavy.py'),
                     inshp, intif, outtif,
                     param4, param5, param6,
                    param7 = as.numeric(Sys.getenv('COLA_NCORES')),
                    param8 = 'None', tempFolder = rootPath){

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
crk_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                pyscript = system.file(package = 'cola', 'python/crk.py'),
                inshp, intif, outtif,
                param4, param5, param6,
                param7 = as.numeric(Sys.getenv('COLA_NCORES')), param8 = 'None'){

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
pri_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                   pyscript = system.file(package = 'cola', 'python/prioritize_core_conn.py'),
                   tif, incrk, inlcc,
                   maskedcsname = paste0(tempfile(), '.tif'),
                   outshppoint, outshppol, outshppatch,
                   outtifpatch, outtif,
                   param7 = 0.5,
                   param8 = 1000){

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

  print(cmd_prio <- gsub(fixed = TRUE, '\\', '/', cmd_prio))


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE, ignore.stdout = TRUE), error = function(e) e$message)
  print(intCMD)
  return(list(tif = ifelse(file.exists(outtif), outtif, NA),
              shp = ifelse(file.exists(outshp), outshp, NA),
              log = intCMD) )
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
cdpop_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"), tif, incrk, inlcc,
                     maskedcsname = paste0(tempfile(), '.tif'),
                     outshp, outtif,
                     param5 = 0.5, param6 = 1000){


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
crk_compare_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/crk_compare.py'),
                           intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None'){

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
lcc_compare_py <- function(py = Sys.getenv("COLA_PYTHON_PATH"),
                           pyscript = system.file(package = 'cola', 'python/lcc_compare.py'),
                           intif, intifs,
                           outcsvabs, outcsvrel,
                           outpngabs, outpngrel,
                           outfolder,
                           inshp = 'None',
                           shpfield = 'None'){

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

