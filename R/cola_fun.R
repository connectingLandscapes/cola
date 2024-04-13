cola_dss <- function(launch.browser = TRUE) ) {
  #app_path <- system.file("shiny", package = "wallace")
  dssLocation <- system.file('app', package = "cola")
  #knitcitations::cleanbib()
  #options("citation_format" = "pandoc")
  #preexisting_objects <- ls(envir = .GlobalEnv)
  on.exit(rm(list = setdiff(ls(envir = .GlobalEnv), preexisting_objects), envir = .GlobalEnv))
  return( shiny::runApp(dssLocation, launch.browser = launch.browser) )
}


#' Run CDPOP model
#'
#' Run CDPOP model
#' @return The temperature in degrees Celsius
#' @examples
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export

runCDPOP <- function(py, datapath = tempFolder){
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
  src <- '/home/shiny/connecting-landscapes/lib/CDPOP/src/CDPOP.py'
  datapath <- tempFolder # datapath = tempFolder
  vars <- paste0('invars.csv') # Only file name
  timeMarkCDPOP <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE))
  cdpopPath <- paste0('cdpopout_', timeMarkCDPOP, '__')
  cdpopPath <- 'cdpopout'
  (cmd <- paste0(py, ' ', src, ' ', datapath, ' ', vars, ' ', cdpopPath))

  # setwd(tempFolder)
  #file.copy('inputvars.csv', 'in.csv')
  #file.copy('xyA.csv', 'xy.csv')
  #intCMD <- tryCatch(system(cmd, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)

  newFiles <- list.files(path = tempFolder, recursive = TRUE)
  newFiles <- grep(pattern = cdpopPath, newFiles, value = TRUE)
  return(list(newFiles = newFiles, cdpopPath = cdpopPath))
}


runS2RES <- function(py, src,
                     intif, outtif,
                     param3, param4, param5, param6 = nCores, param7, param8 = 'None'){
  # param3 = 0
  # param4 =  100
  # param5 = 100
  # param6 = 1
  # param7 = -9999

  src <- '/home/shiny/connecting-landscapes/src/s2res.py'
  datapath <- tempFolder # datapath = tempFolder

  timeMarkCDPOP <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE))

  (cmd_s2res <- paste0(py, ' ', src, ' ', intif, ' ', outtif, ' ',
                       format(param3, scientific=F), ' ', format(param4, scientific=F),
                       ' ', format(param5, scientific=F), ' ', format(param6, scientific=F),
                       ' ', format(param7, scientific=F), ' ', param8))

  intCMD <- tryCatch(system(cmd_s2res, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  return(file = ifelse(file.exists(outtif), outtif, NA))
}


points_shp <- function(py, intif, outshp, param3, param4, param5, param6 = 'None'){
  # param3 = 2
  # param4 =  95
  # param5 = 50
  src <- '/home/shiny/connecting-landscapes/src/create_source_points.py'
  datapath <- tempFolder # datapath = tempFolder
  (cmd_pts <- paste0(py, ' ', src, ' ', intif, ' ', outshp, ' ',
                     format(param3, scientific=F), ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     param6))

  print(cmd_pts)
  intCMD <- tryCatch(system(cmd_pts, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  return(file = ifelse(file.exists(outshp), outshp, NA))
}

cdmat_py <- function(py, inshp, intif, outcsv, param3, param4 = nCores, param5 = 'None'){
  # param3 = 25000
  # create_cdmat.py
  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)

  src <- '/home/shiny/connecting-landscapes/src/create_cdmat.py'
  (cmd_cdmat <- paste0(py, ' ', src, ' ', inshp, ' ', intif, ' ', outcsv,
                       ' ', format(param3, scientific=F),
                       ' ', param4, ' ', param5))

  intCMD <- tryCatch(system(cmd_cdmat, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  return(file = ifelse(file.exists(outcsv), outcsv, NA))
}



lcc_py <- function(py, inshp, intif, outtif, param4, param5,
                   param6, param7 = nCores, param8 = 'None'){
  # param3 = 25000
  # [1] source points: Spatial point layer (any ORG driver), CSV (X, Y files), or *.xy file
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (should be in meters*)
  # [5] corridor smoothing factor (in number of cells)
  # [6] corridor tolerance (in cost distance units)

  src <- '/home/shiny/connecting-landscapes/src/lcc.py'

  if(require('gdalUtilities')){
    gi <- capture.output(gdalUtilities::gdalinfo(intif))
    nd0 <- as.numeric(
      gsub("[^-[:alnum:]]", "",
           gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)
           )
      )
    )
    # } else if(require('gdalUtils')){
    #   gi <- capture.output(gdalUtils::gdalinfo(intif))
    #   nd0 <- as.numeric(
    #     gsub("[^-[:alnum:]]", "",
    #          gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)
    #          )
    #     )
    #   )
  }


  (cmd_lcc <- paste0(py, ' ', src, ' ', inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     format(param6, scientific=F), " ",
                     format(param7, scientific=F), " ",
                     param8))

  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  return(file = ifelse(file.exists(outtif), outtif, NA))
}


lcc_py2 <- function(py, inshp, intif, outtif, param4, param5,
                    param6, param7 = nCores, param8 = 'None', tempFolder = rootPath){

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

  src <- '/home/shiny/connecting-landscapes/src/lcc_hdf5_v7.py'

  ## get NoData

  if(require('gdalUtilities')){
    gi <- capture.output(gdalUtilities::gdalinfo(intif))
    nd0 <- as.numeric(
      gsub("[^-[:alnum:]]", "",
           gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)
           )
      )
    )
    # } else if(require('gdalUtils')){
    #   gi <- capture.output(gdalUtils::gdalinfo(intif))
    #   nd0 <- as.numeric(
    #     gsub("[^-[:alnum:]]", "",
    #          gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)
    #          )
    #     )
    #   )
  }



  (cmd_lcc <- paste0(py, ' ', src, ' ', inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     format(param6, scientific=F), " ",
                     format(param7, scientific=F), " ",
                     param8, " ",
                     h5file1, " ",
                     h5file2, " ",
                     '50'
  ))

  print(cmd_lcc)

  intCMD <- tryCatch(system(cmd_lcc, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  file.remove(c(h5file1, h5file2))
  return(file = ifelse(file.exists(outtif), outtif, NA))
}

crk_py <- function(py, inshp, intif, outtif, param4,
                   param5, param6, param7 = nCores, param8 = 'None'){
  # [1] source points
  # [2] resistance surface
  # [3] output file name
  # [4] distance threshold (in cost distance units)
  # [5] kernel shape (linear, gaussian)
  # [5] kernel volume

  src <- '/home/shiny/connecting-landscapes/src/crk.py'
  (cmd_crk <- paste0(py, ' ', src, ' ', inshp, ' ', intif, ' ', outtif, ' ',
                     format(param4, scientific=F), ' ',
                     format(param5, scientific=F), ' ',
                     format(param6, scientific=F), ' ',
                     format(param7, scientific=F), ' ', param8))

  intCMD <- tryCatch(system(cmd_crk, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  return(file = ifelse(file.exists(outtif), outtif, NA))
}



pri_py <- function(py, tif, incrk, inlcc,
                   maskedcsname = paste0(tempfile(), '.tif'),
                   outshp, outtif,
                   param5 = 0.5, param6 = 1000){

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

  src <- '/home/shiny/connecting-landscapes/src/prioritize_core_conn.py'

  # incrk <- '/data/temp/RY2024011519163805file176c0d742621ed/out_crk_IF2024011520212105file176c0d2f0a3994.tif'
  # inlcc <- '/data/temp/RY2024011519163805file176c0d742621ed/out_lcc_GG2024011519164505file176c0d5f4b6267.tif'
  # outshp <- '/data/temp/LA2024011522400305file1a50376cab710d//out_pri_EX2024011522411005file1a5037c8aef66.shp'
  # outtif <- '/data/temp/LA2024011522400305file1a50376cab710d//out_pri_EX2024011522411005file1a5037c8aef66.tif'

  (cmd_prio <- paste0(py, ' ',
                      src, ' ',
                      tif, ' ',
                      incrk, ' ', inlcc, ' ',
                      maskedcsname, ' ',
                      outshp, ' ', outtif, ' ',
                      format(param5, scientific=F), " ",
                      format(param6, scientific=F)))


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  print(intCMD)
  return(list(tif = ifelse(file.exists(outtif), outtif, NA),
              shp = ifelse(file.exists(outshp), outshp, NA)) )
}


cdpop_py <- function(py, tif, incrk, inlcc,
                     maskedcsname = paste0(tempfile(), '.tif'),
                     outshp, outtif,
                     param5 = 0.5, param6 = 1000){


  #(base) C:\Users\Admin\dockerdata\CDPOP\src>
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

  src <- '/home/shiny/connecting-landscapes/src/prioritize_core_conn.py'


  (cmd_prio <- paste0(py, ' ',
                      src, ' ',
                      tif, ' ',
                      incrk, ' ', inlcc, ' ',
                      maskedcsname, ' ',
                      outshp, ' ', outtif, ' ',
                      format(param5, scientific=F), " ",
                      format(param6, scientific=F)))


  intCMD <- tryCatch(system(cmd_prio, intern = TRUE, ignore.stdout = TRUE), error = function(e) NULL)
  print(intCMD)
  return(list(tif = ifelse(file.exists(outtif), outtif, NA),
              shp = ifelse(file.exists(outshp), outshp, NA)) )
}
