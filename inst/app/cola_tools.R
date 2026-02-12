## Cola utils

# options(scipen = 999)
# options(scipen = 9)

## Errors -- only works for shiny server
if (TRUE){

  dataFolder <- ifelse(exists('dataFolder'), dataFolder, '')
  path_error <- ifelse(exists('path_error'), path_error, '')
  logFilePath <<- base::paste0(dataFolder, '/cola_logFolders.txt')
  allLogs <- base::list.files(path = path_error, pattern = 'cola|connec')

  if( length(allLogs) != 0 ){
    # cleanMemory(logFilePath)

    validLogs <- unname(
      na.omit(
        sapply(allLogs, simplify = TRUE, function(x){
          # x = allLogs[1]
          y <- tryCatch(read.delim(file.path(path_error, x), quote = "")[, 1], error = function(e) NULL)
          if (any(grep('[E|e]rror', y)) & !is.null(y) ){
            return (x)
          } else {
            return (NA)
          }
        })
      ))
  } else {
    ## Temporal variable --- not to use outside the server
    validLogs <- c('')
  }
}


## time stamp -----
## Create temporal ID
sessionIDgen <- function(letter = TRUE, sep = '', short = TRUE, folder = FALSE, only3 = FALSE){
  (tempID <- basename(tempfile()))
  (timeMark <- gsub('[[:punct:]]| ', '', format(as.POSIXct(Sys.time(), tz="CET"), tz="America/Bogota",usetz=TRUE)))

  if( short){
    (sessionID <- timeMark)
  } else {
    (sessionID <- paste0( timeMark, sep, tempID ))
  }

  if (letter){
    sessionID <- paste0(sample(LETTERS, 1), sample(LETTERS, 1), sample(LETTERS, 1),
                        sep, sessionID)
  }

  if(folder){
    sessionID <- paste0('cola', sessionID)
  }

  if(only3){
    sessionID <- paste0(sample(LETTERS, 1), sample(LETTERS, 1), sample(LETTERS, 1))
  }
  return(sessionID)
}


## Clean files
cleanMemory <- function(logFilePath){
  dfm <- data.frame(Rtmp = tempdir(), tempFolder = tempFolder)

  if(file.exists(logFilePath)){
    logDF <- tryCatch(read.csv(logFilePath), error = function(e) NULL)
  } else {
    logDF <- NULL
  }

  logDF <- rbind(logDF, dfm)

  openFolders <- dir.exists(logDF$Rtmp)
  sapply(logDF$tempFolder[!openFolders], unlink, recursive = TRUE)
  logDF <- logDF[dir.exists(logDF$Rtmp), ]
  write.csv(x = logDF, logFilePath, row.names = FALSE)
}


## try deleting files
delFiles <- function(...){
  invisible(suppressWarnings(
    tryCatch(file.remove(c(...), recursive = TRUE),
             error = function(e) NULL)
  ))
  invisible(suppressWarnings(
    tryCatch(unlink(c(...), force = TRUE, recursive = TRUE),
             error = function(e) NULL)
  ))
}

## Evals if the raster is GEO or PROJ
isProjected <- function(rastPath, details = FALSE){
  #rastPath = 'C:/temp/cola/colaRKW2024081218272505//in_points_CHL2024081218292405.tif'
  cat(' +++ Evaluaring if is proj: ', rastPath, '\n')

  if (require('gdalUtilities')){
    gi <- gdalUtilities::gdalinfo(rastPath, quiet = TRUE)
    #gi <- capture.output(gdalUtilities::gdalinfo(rastPath))

    if (details){
      print(gi)
    }

    g2 <- strsplit(x = gi, split = '\n')[[1]]
    (isProj <- (length( grep('^GEOGCRS', g2) ) == 0) & any( grep('PROJCRS', g2) ))
  } else{
    isGEO <- terra::is.lonlat(rastPath);
    (isProj <- !isGEO)
  }
  return(isProj)
}


## Resample raster layer
resampIfNeeded <- function(rastPath){
  # rastPath <- '/data/temp/ZS2023111311113105file4f6823209882/in_surface_fixed_QK2023111311142905file4f687fa39935.tif'
  # rastPath <- '/data/temp/LI2024011611295205file1a5d191fa71355/in_lcc_CQ2024011611304405file1a5d19779aa294.tif'
  # rastPath <- '/home/shiny/preCanopyClass30mCost.tif'

  if(!require(terra)){
    stop('Need terra library')
  }

  r <- terra::rast(rastPath)
  (totpixels <- terra::ncol(r) * terra::nrow(r))
  #names(r)
  #rasname <- r@pnt@.xData$filenames()
  (resamPath <- gsub(x = rastPath,
                     '.tif$', '_resam.tif'))
  if(totpixels > options('COLA_VIZ_THREs_PIX')){
    #if(file.exists()){
    if (!file.exists(resamPath)){

      gdalUtilities::gdalwarp(srcfile = rastPath,
                              ts = c(options('COLA_VIZ_RES_NCOL'), options('COLA_VIZ_RES_NROW')),
                              dstfile = resamPath)
      cat(paste(' ---- >>>> Resampling raster to ',
                options('COLA_VIZ_RES_NCOL'),' - ',
                options('COLA_VIZ_RES_NROW'), '\n'))
    }
    #}
    return(resamPath)
  } else {
    return(rastPath)
  }
}

canLeafletPlotIt <- function(x){
  # x is rast file path or SpatRast
  epsg4326 <- "+proj=longlat +datum=WGS84 +no_defs"
  epsg3857 <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs" # nolint
  require(terra)
  # x <- system.file(package = 'cola', 'sampledata/sampleTif.tif')

  if(class(x) == 'character'){
    if(file.exists(x)){
      x <- terra::rast(x)
    } else{
      stop(' Raster file not existing. Failed at check leaflet compatibility')
    }
  }

  bounds <- tryCatch(terra::ext(
    terra::project(
      terra::project(
        terra::as.points(terra::ext(x), crs = terra::crs(x)),
        epsg3857
      ),
      epsg4326
    )), error = function(e) {print(e); NULL})

  if (is.null(bounds)){
    # Can't plot in dss
    return(FALSE)
  } else {
    # Works
    return(TRUE)
  }
}

draws2Features <- function(polDraw, distLineBuf = NULL, rastCRS, crs2assign = 4326){

  ## Create list of empty features types
  polList <- list(
    line = list(),
    polygon = list(),
    rectangle = list(),
    circle = list(),
    marker= list()
  )

  #(ss <- load(file = '/data/temp/draw.RData')); polDraw$geometry$coordinates  # polDraw
  # unlist(lapply(polDraw$features, function(x) x$geometry$type))
  # unlist(lapply(polDraw$features, function(x) x$properties$feature_type))

  # Iterate
  pol2Rast <- NULL;
  for(l in 1:length(polDraw$features) ){
    # l = 4
    feat <- polDraw$features[[l]]
    (featType <- feat$properties$feature_type)
    # print(featType)

    if(featType == 'polyline'){
      #lapply(polDraw$geometry$coordinates, length)
      #coords <- do.call(rbind, polDraw$geometry$coordinates[[1]])
      coords <- matrix(unlist(feat$geometry$coordinates),
                       ncol = 2, byrow = TRUE)
      # coords <- do.call(c, feat$geometry$coordinates)
      # pts_geo <- SpatialPoints(coords, proj4string = CRS('EPSG:4326'))
      # pts_geo <-  sf::st_as_sf(as.data.frame(coords), coords = c(1,2))
      # st_crs(pts_geo) <- crs('EPSG:4326')
      # #pts <- spTransform(pts_geo, CRSobj = rastCRS)
      # pts <- sf::st_transform(pts_geo, to = rastCRS)
      # line_ <- as(pts,"SpatialLines")

      line_ <-  (st_linestring(coords))
      #st_crs(line_) <- crs('EPSG:4326')

      # plot(line_, axes = TRUE)
      lnbuf <- st_sfc(st_polygon(st_buffer(line_, 0.000001)), crs = crs2assign)
      # lnbuf <- SpatialPolygonsDataFrame(gBuffer(line_, width = 1),
      #                                   data = data.frame(ID = l),
      #                                   match.ID = FALSE)
      lnbuf <- st_buffer(sf::st_transform(lnbuf, crs = rastCRS), dist = distLineBuf)

      pol2Add <- lnbuf
      # plot(pol2Add, axes = TRUE); plot(rast(rastPath), add = TRUE); plot(pol2Add, add = TRUE)
    }


    if(featType == 'polygon'){
      #(ss <- load(file = '/data/temp/2pols.RData'));
      coordx <- feat$geometry$coordinates[[1]]  # polDraw
      # str(coordx)
      cx <- (as.matrix.data.frame(do.call(rbind, coordx)))
      #cx <- rbind(cx, cx[nrow(cx), ])

      pol2save_ <- (st_polygon(list(cx)))
      pol2save <- sf::st_transform(st_sfc(pol2save_, crs = crs2assign),
                                   crs = rastCRS)

      # pol2save <- SpatialPolygonsDataFrame(
      #   data = data.frame(ID = l),
      #   SpatialPolygons(
      #     Srl = list( Polygons(
      #       list(Polygon(as.matrix(cx))), 1) ),
      #     proj4string = CRS('EPSG:4326'))
      # ) # plot(pol2save)

      #pol2save <- spTransform(pol2save, CRSobj = rastCRS)
      #pol2save <- sf::st_transform(pol2save, to = rastCRS)

      #pol2save <- Polygon(as.matrix(cx))
      # if(length(polList$polygon) == 0){
      #   polList$polygon <- pol2save
      # } else {
      #   polList$polygon <- bind(polList$polygon, pol2save)
      # }
      pol2Add <- pol2save
    }


    if(featType == 'rectangle'){
      coordx <- feat$geometry$coordinates[[1]]  # polDraw
      # str(coordx)
      cx <- (as.matrix.data.frame(do.call(rbind, coordx)))

      pol2save_ <- (st_polygon(list(cx)))
      pol2save <- st_sfc(pol2save_, crs = crs2assign)
      pol2save <- sf::st_transform(pol2save, crs = rastCRS)

      # pol2save <- SpatialPolygonsDataFrame(
      #   data = data.frame(ID = l),
      #   SpatialPolygons(
      #     Srl = list( Polygons(
      #       list(Polygon(as.matrix(cx))), 1) ),
      #     proj4string = CRS('EPSG:4326'))
      # )
      # # plot(pol2save)
      # pol2save <- sf::st_transform(pol2save, crs = rastCRS)

      #pol2save <- Polygon(as.matrix(cx))
      # if(length(polList$rectangle) == 0){
      #   polList$rectangle <- pol2save
      # } else {
      #   polList$rectangle <- bind(polList$rectangle, pol2save)
      # }

      pol2Add <- pol2save
    }


    if(featType == 'circle'){
      # coordx <- feat$geometry$coordinates  # polDraw
      # cx <- (as.matrix.data.frame(do.call(cbind, coordx)))
      #
      # pts_geo <- SpatialPoints(cx, proj4string = CRS('EPSG:4326'))
      # pts <- sf::st_transform(pts_geo, CRSobj = rastCRS)
      # # plot(pts, axes = TRUE)
      # circ <- SpatialPolygonsDataFrame(gBuffer(pts, width = feat$properties$radius),
      #                                  match.ID = FALSE,
      #                                  data = data.frame(ID = l))

      coords <- matrix(unlist(feat$geometry$coordinates),
                       ncol = 2, byrow = TRUE)

      pt_ <-  st_sfc(st_point(coords), crs = crs2assign)
      pt_ <- sf::st_transform(pt_, crs = rastCRS)

      circ <- st_sfc((st_buffer(pt_, dist = feat$properties$radius)))

      # if(length(polList$circle) == 0){
      #   polList$circle <- circ
      # } else {
      #   polList$circle <- bind(polList$circle, circ)
      # }
      pol2Add <- circ
    }

    if(featType == 'marker'){
      coords <- do.call(cbind, feat$geometry$coordinates)
      # pts_geo <- SpatialPoints(coords, proj4string = CRS('EPSG:4326'))

      # pts_geo <-  sf::st_as_sf(as.data.frame(coords), coords = c(1,2))
      # st_crs(pts_geo) <- crs('EPSG:4326')
      # pts <- sf::st_transform(pts_geo, CRSobj = rastCRS)
      #
      # # plot(pts, axes = TRUE)
      # ptbuf <- SpatialPolygonsDataFrame(gBuffer(pts, width = 1),
      #                                   match.ID = FALSE,
      #                                   data = data.frame(ID = l))

      pt_ <-  st_sfc(st_point(coords), crs = crs2assign) # coords = c("x","y")
      pt_ <- sf::st_transform(pt_, crs = rastCRS)

      ptbuf <- st_sfc((st_buffer(pt_, dist = 1)))

      pol2Add <- ptbuf
    }

    #plot(pol2Add, main = l, axes = TRUE)

    ## Compile
    if(length(pol2Rast) == 0){
      pol2Rast <- pol2Add
    } else {
      #pol2Rast <- bind(pol2Rast, pol2Add)
      pol2Rast <- c(pol2Rast, pol2Add)
    }
    # plot(raster(rastPath), main = l); plot(pol2Rast, add = TRUE)
    #plot(pol2Rast, main = l, axes = TRUE)
  }


  # pol2Rastx <- st_sf(data.frame(a = 1:length(pol2Rast), pol2Rast))
  pol2Rastx <- st_as_sf( pol2Rast)
  return(pol2Rastx)
}




tif2rsg <- function(path, outdir = NULL){
  ## Convert TIF files into ASC, then modify header to match RSG format.
  ## Needs the raster path
  ## If no oudir provided, input folder is used
  ## Returns the RSG filename


  # path = 'N:/My Drive/connectivity-nasa-USFSIP/01_original-data/Sabah_example_CDPOP+UNICOR/roads.tif'
  # path = 'N:/My Drive/connectivity-nasa-USFSIP/01_original-data/Sabah_example_CDPOP+UNICOR/'
  # outdir <- path
  if (file.exists(path)){
    tif <- raster(path)

    (outdir <- ifelse(is.null(outdir), yes = dirname(path), no = gsub(pattern = basename(outdir), x = outdir, replacement = '') ))
    (outdir <- ifelse(outdir == '', '.', outdir))
    fname <- paste0(outdir, '/', tools::file_path_sans_ext( basename(path)), '_temp.asc')
    fname_rsg <- paste0(outdir, '/', tools::file_path_sans_ext( basename(path)), '.rsg')

    terra::writeRaster(x = tif, format = 'ascii', overwrite = TRUE, filename = fname )

    asc <- read.delim(fname, header = FALSE)
    asc[1:5, 1] <- gsub('[[:blank:]]', '\t', gsub(' $', '', tolower(asc[1:5, 1])))
    write.table(asc, file = fname_rsg, col.names = FALSE, row.names = FALSE, quote = FALSE)
    return(fname_rsg)
  } else{
    stop(' Error: file not found')
  }
}




#case0 <- read.delim('C:/Users/gonza/dockerdata/small_test.rip', header = FALSE, sep = '|')
#setwd('C:/Users/Admin/dockerdata/UNICOR/unicor')
origRip <- "
Session_label    case1A
Grid_Filename    size1_side3px_10totpix.rsg
XY_Filename    sabah_10.xy
Use_Direction	FALSE
Type_Direction	FlowAcc
Use_Resistance	TRUE
Barrier_or_U_Filename	completebarr_test2.txt
Direction_or_V_Filename	flowaccu_test2.txt
Speed_To_Resistance_Scale	0;10
Use_ED_threshold	False
ED_Distance   10000
Edge_Type	normal
Transform_function	linear
Const_kernal_vol	False
Kernel_volume   10000
Edge_Distance    1000000
Number_of_Processes    4
KDE_Function    Gaussian
KDE_GridSize    2
Number_of_Categories    5
Save_Path_Output    TRUE
Save_IndividualPaths_Output    FALSE
Save_GraphMetrics_Output    FALSE
Save_KDE_Output    TRUE
"
RIP <- c('Session_label    SESSIONLABEL',
         'Grid_Filename    RSGFILE',
         'XY_Filename    XYFILE',
         'Use_Direction	FALSE',
         'Type_Direction	TYPEDIRECTION',
         'Use_Resistance	USERESISTANCE',
         'Barrier_or_U_Filename	BARRIER',
         'Direction_or_V_Filename	DIRECTION',
         'Speed_To_Resistance_Scale	SPEED',
         'Use_ED_threshold    USEEDTHRESHOLD',
         'ED_Distance   EDDISTANCE',
         'Edge_Type	EDGETYPE',
         'Transform_function	TRANSFUNCTION',
         'Const_kernal_vol	CONSTKERNALVOL',
         'Kernel_volume   KERNELVOLUME',
         'Edge_Distance    EDGEDISTANCE',
         'Number_of_Processes    NPROCESSES',
         'KDE_Function    KDEFUNCTION',
         'KDE_GridSize    KDEGRIDSIZE',
         'Number_of_Categories    NOFCATEGORIES',
         'Save_Path_Output    SAVEPATHOUTPUT',
         'Save_IndividualPaths_Output    SAVEINDIVIDUALPATHSOUTPUTS',
         'Save_GraphMetrics_Output    SAVEGRAPHS',
         'Save_KDE_Output    SAVEKDEOUTPUT',
         'Save_Category_Output    SAVECATEGORYOUTPUT',
         'Save_CDmatrix_Output    SAVECDMATRIXOUTPUT')


createrip <- function(caseName = NULL, Grid_Filename = NULL, XY_Filename = NULL,
                      Use_Direction = NULL, Type_Direction = NULL, Use_Resistance = NULL,
                      Barrier_or_U_Filename = NULL, Direction_or_V_Filename = NULL,
                      Speed_To_Resistance_Scale= NULL, Use_ED_threshold = NULL, ED_Distance = NULL,
                      Edge_Type = NULL, Transform_function = NULL, Const_kernal_vol = NULL,
                      Kernel_volume = NULL, Edge_Distance = NULL, Number_of_Processes = NULL,
                      KDE_Function = NULL, KDE_GridSize = NULL, Number_of_Categories = NULL,
                      Save_Path_Output = NULL, Save_IndividualPaths_Output = NULL,
                      Save_GraphMetrics_Output = NULL, Save_KDE_Output = NULL,
                      Save_Category_Output = NULL, Save_CDmatrix_Output = NULL){

  ripFile <- ripTemplate

  if (is.null(caseName)){
    ripFile['Session_label', 1] <- gsub('SESSIONLABEL', caseName)
  }
  if (is.null(Grid_Filename)){
    ripFile['Grid_Filename', 1] <- gsub('RSGFILE', Grid_Filename)
  }
  if (is.null(XY_Filename)){
    ripFile['XY_Filename', 1] <- gsub('XYFILE',XY_Filename)
  }
  if (is.null(Use_Direction)){
    ripFile['Use_Direction', 1] <- gsub('FALSE', Use_Direction)
  }
  if (is.null(Type_Direction)){
    ripFile['Type_Direction', 1] <- gsub('TYPEDIRECTION',  Type_Direction)
  }
  if (is.null(Use_Resistance)){
    ripFile['Use_Resistance', 1] <- gsub('USERESISTANCE',  Use_Resistance)
  }
  if (is.null(Barrier_or_U_Filename)){
    ripFile['Barrier_or_U_Filename', 1] <- gsub('BARRIER',  Barrier_or_U_Filename)
  }
  if (is.null(Direction_or_V_Filename)){
    ripFile['Direction_or_V_Filename', 1] <- gsub('DIRECTION', Direction_or_V_Filename )
  }
  if (is.null(Speed_To_Resistance_Scale)){
    ripFile['Speed_To_Resistance_Scale', 1] <- gsub('SPEED',  Speed_To_Resistance_Scale)
  }
  if (is.null(Use_ED_threshold)){
    ripFile['Use_ED_threshold', 1] <- gsub('USEEDTHRESHOLD', Use_ED_threshold )
  }
  if (is.null(ED_Distance)){
    ripFile['ED_Distance', 1] <- gsub('EDDISTANCE',  ED_Distance)
  }
  if (is.null(Edge_Type)){
    ripFile['Edge_Type', 1] <- gsub('EDGETYPE',  Edge_Type)
  }
  if (is.null(Transform_function)){
    ripFile['Transform_function', 1] <- gsub('TRANSFUNCTION', Transform_function )
  }
  if (is.null(Const_kernal_vol)){
    ripFile['Const_kernal_vol', 1] <- gsub('CONSTKERNALVOL',  )
  }
  if (is.null(Kernel_volume)){
    ripFile['Kernel_volume', 1] <- gsub('KERNELVOLUME',  Const_kernal_vol)
  }
  if (is.null(Edge_Distance)){
    ripFile['Edge_Distance', 1] <- gsub('EDGEDISTANCE', Edge_Distance )
  }
  if (is.null(Number_of_Processes)){
    ripFile['Number_of_Processes', 1] <- gsub('NPROCESSES', Number_of_Processes )
  }
  if (is.null(KDE_Function)){
    ripFile['KDE_Function', 1] <- gsub('KDEFUNCTION', KDE_Function)
  }
  if (is.null(KDE_GridSize)){
    ripFile['KDE_GridSize', 1] <- gsub('KDEGRIDSIZE', KDE_GridSize )
  }
  if (is.null(Number_of_Categories)){
    ripFile['Number_of_Categories', 1] <- gsub('NOFCATEGORIES',  Number_of_Categories)
  }
  if (is.null(Save_Path_Output)){
    ripFile['Save_Path_Output', 1] <- gsub('SAVEPATHOUTPUT', Save_Path_Output )
  }
  if (is.null(Save_IndividualPaths_Output)){
    ripFile['Save_IndividualPaths_Output', 1] <- gsub('SAVEINDIVIDUALPATHSOUTPUTS', Save_IndividualPaths_Output )
  }
  if (is.null(Save_GraphMetrics_Output)){
    ripFile['Save_GraphMetrics_Output', 1] <- gsub('SAVEGRAPHS', Save_GraphMetrics_Output )
  }
  if (is.null(Save_KDE_Output)){
    ripFile['Save_KDE_Output', 1] <- gsub('SAVEKDEOUTPUT', Save_KDE_Output )
  }

  if (is.null(Save_Category_Output)){
    ripFile['Save_Category_Output', 1] <- gsub('SAVECATEGORYOUTPUT', Save_Category_Output )
  }
  if (is.null(Save_CDmatrix_Output)){
    ripFile['Save_CDmatrix_Output', 1] <- gsub('SAVECDMATRIXOUTPUT', Save_CDmatrix_Output)
  }
  write.table(ripFile, file = 'caseName.rip', col.names = FALSE, row.names = FALSE, quote = FALSE)
}


ripTemplate <- data.frame(RIP, row.names = gsub(' .+|\t.+', '', RIP) )



#write.table(ripFile, file = 'caseName.rip', col.names = FALSE, row.names = FALSE, quote = FALSE)



### Orig fitRaster2cola, cell size and nodata
fitRaster2cola <- function(inrasterpath, outrasterpath = NULL){
  # setwd('/home/shiny')
  # inrasterpath = 'orig_tifs/size6.tif'
  # outrasterpath = 'size6.tif'
  # inrasterpath = '/data/tempR/colaELU2024080412561105//in_points_JLT2024080412564005.tif'
  # outrasterpath = '/data/tempR/colaELU2024080412561105//in_points_fixed_JLT2024080412564005.tif'

  inraster <- inrasterpath
  outraster <- outrasterpath

  outraster0 <- NA

  if( ! (file.exists(inraster) | is.na(inraster) | is.null(inraster)) ){
    stop( print(paste('  >>> Infile not found - ', inraster)))
  }

  if( is.null(outraster) ){
    outraster <- paste0(tools::file_path_sans_ext(inraster), 'out',
                        basename(tempfile()) ,'.tif')
  } else {
    if( !file.exists(inraster)){
      stop( print(paste('  >>> Outfile not found - ', inraster)))
    }
  }


  if (!isProjected(inraster)){
    return(NA)
    print(' Raster MUST have proceted coordinate system. Please update your layer.')
  } else {

    if( require('gdalUtilities') ){
      # options(scipen = 999)
      # options(scipen = 9)

      gi <- capture.output(gdalUtilities::gdalinfo(inraster))
      (nd0 <- as.numeric(gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE))))
      if (length(nd0) == 0) { nd0 <- 0}
      nd <- (length(nd0) == 1) & (nd0 == -9999)
      if ( is.nan(nd0) | is.na(nd0)){ nd <- FALSE }

      pixsize0 <- unlist(strsplit(split = ',', gsub('Pixel Size = \\(|\\)', '',
                                                    grep('Pixel Size', gi, value = TRUE))))
      options(digits = max(nchar(pixsize0)))
      (pixsize <- abs(as.numeric(pixsize0 )))
      (ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] ))

      #if( !( nd & ps ) ) {
      if( !( ps ) ) {
        gdalUtilities::gdalwarp(srcfile = inraster, dstfile = outraster,
                                # ot = 'Float64',
                                #srcband  = 1,
                                tr = rep(max(pixsize), 2)
                                #, dstnodata = -9999
        )
        outraster0 <- outraster
        # }
        # else if( nd & !ps ){
        #   cat ( ' \n Converting raster: Changing no data \n' )
        #   gdalUtilities::gdalwarp(srcfile = inraster, dstfile = outraster,
        #                           #srcband  = 1,
        #                           #srcnodata = nd0,
        #                           ot = 'Float64',
        #                           tr = rep(max(pixsize), 2)
        #                           )
        #   outraster0 <- outraster
        #
        # } else if (!nd & ps ){
        #   cat ( ' \n Converting raster: Changing no data \n' )
        #   gdalUtilities::gdalwarp(srcfile = inraster,
        #                           dstfile = outraster,
        #                           #b = 1,
        #                           ot = 'Float64'
        #                           #, srcnodata = nd0 ,
        #                           #dstnodata = -9999
        #   )
        #   outraster0 <- outraster
        # } else if (!ps ){
        #   cat ( '  Making pixel squared' )
        #   gdalUtilities::gdalwarp(srcfile = inraster,
        #                           dstfile = outraster,
        #                           #b = 1,
        #                           ot = 'Float64'
        #                           #, srcnodata = nd0 ,
        #                           #dstnodata = -9999
        #   )
        #   outraster0 <- outraster
      } else{
        outraster0 <- inraster
      }

    } else if(require('raster')){
      # print(3)

      options(digits = 20)
      rx <- raster(inraster)
      nd0 <- rx@file@nodatavalue
      if (length(nd0) == 0) { nd0 <- 0}
      nd0 <- nd0 == -9999

      (pixsize <- res(rx))
      ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )
      if( !( nd & ps ) ) {
        templ <- raster(  crs = rx@crs, res = rep(max(res(rx)), 2), ext = extent(rx))
        raster::resample(x = rx, y = templ, filename = outraster,
                         # NAflag = -9999,
                         overwrite = FALSE)
        outraster0 <- outraster
      } else{
        outraster0 <- inraster
      }
      # } else if(require('raster')){
      #     # print(3)
      #
      #     options(digits = 20)
      #     rx <- raster(inraster)
      #     nd <- rx@file@nodatavalue == -9999
      #
      #     (pixsize <- res(rx))
      #     ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )
      #     if( !( nd & ps ) ) {
      #       templ <- raster(  crs = rx@crs, res = rep(max(res(rx)), 2), ext = extent(rx))
      #       raster::resample(x = rx, y = templ, filename = outraster, NAflag = -9999, overwrite = FALSE)
      #       outraster0 <- outraster
      #     } else{
      #       outraster0 <- inraster
      #     }
      #   }
    }

    options(digits = 5)
    return(outraster0)

  }

}

### Only pix size
fitRaster2colaOnlyPxSize <- function(inrasterpath, outrasterpath = NULL){
  # setwd('N:/Mi unidad/git/connecting-landscapes/performance-tests/inputs')

  # inrasterpath = 'orig_tifs/size6.tif'
  # outrasterpath = 'size6.tif'

  # setwd('/home/shiny')
  # inrasterpath = 'preCanopyClass30mCost.tif'
  # outrasterpath = 'preCanopyClass30mCost_OUT.tif'

  inraster <- inrasterpath
  outraster <- outrasterpath

  outraster0 <- NA

  if( ! (file.exists(inraster) | is.na(inraster) | is.null(inraster)) ){
    stop( print('  >>> Infile not found - ', inraster))
  }

  if( is.null(outraster) ){
    outraster <- paste0(tools::file_path_sans_ext(inraster), 'out',
                        basename(tempfile()) ,'.tif')
  } else {
    if( !file.exists(inraster)){
      stop( print('  >>> Outfile not found - ', inraster))
    }
  }

  # if (require('gdalUtils')){
  #   # print(1)
  #
  #   gi <- gdalUtils::gdalinfo(inraster)
  #
  #   nd0 <- as.numeric(gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)))
  #   nd <- length(nd0) == 1 & nd0 == -9999
  #
  #   pixsize0 <- unlist(strsplit(split = ',', gsub('Pixel Size = \\(|\\)', '', grep('Pixel Size', gi, value = TRUE))))
  #   options(digits = max(nchar(pixsize0)))
  #   (pixsize <- abs(as.numeric(pixsize0 )))
  #   ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )
  #
  #   if( !( ps ) ) {
  #     gdalUtils::gdalwarp(srcfile = inraster, dstfile = outraster,
  #                         tr = rep(max(pixsize), 2))
  #     outraster0 <- outraster
  #   } else {
  #     outraster0 <- inraster
  #   }
  # } else
  if(require('gdalUtilities')){
    # print(2)

    # options(scipen = 999)
    # options(scipen = 9)

    gi <- capture.output(gdalUtilities::gdalinfo(inraster))

    nd0 <- as.numeric(gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)))
    nd <- length(nd0) == 1 & nd0 == -9999

    pixsize0 <- unlist(strsplit(split = ',',
                                gsub('Pixel Size = \\(|\\)', '',
                                     grep('Pixel Size', gi, value = TRUE))))
    options(digits = max(nchar(pixsize0)))
    (pixsize <- abs(as.numeric(pixsize0 )))
    ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )

    if( !( ps ) ) {
      gdalUtilities::gdalwarp(srcfile = inraster,
                              dstfile = outraster,
                              #ot = 'Float64',
                              tr = rep(max(pixsize), 2)
                              #, dstnodata = -9999
      )
      outraster0 <- outraster
    } else{
      outraster0 <- inraster
    }


  } else if(require('raster')){
    # print(3)

    options(digits = 20)
    rx <- raster(inraster)
    nd <- rx@file@nodatavalue == -9999

    (pixsize <- res(rx))
    ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )
    if( !( nd & ps ) ) {
      templ <- raster(  crs = rx@crs, res = rep(max(res(rx)), 2), ext = extent(rx))
      raster::resample(x = rx, y = templ, filename = outraster, overwrite = FALSE)
      outraster0 <- outraster
    } else{
      outraster0 <- inraster
    }
  }

  options(digits = 5)
  return(outraster0)
}

## Orig function
fitRaster2cola0 <- function(inrasterpath, outrasterpath = NULL){
  # setwd('N:/Mi unidad/git/connecting-landscapes/performance-tests/inputs')
  # inrasterpath = 'orig_tifs/size6.tif'
  # outrasterpath = 'size6.tif'
  # inrasterpath = intif
  # outrasterpath = outtif

  inraster <- inrasterpath
  outraster <- outrasterpath
  outraster0 <- inrasterpath


  if( ! (file.exists(inraster) | is.na(inraster) | is.null(inraster)) ){
    stop( print('  >>> Infile not found - ', inraster))
  }

  if( is.null(outraster) ){
    outraster <- paste0(tools::file_path_sans_ext(inraster), 'out',
                        basename(tempfile()) ,'.tif')
  } else {
    if( !file.exists(inraster)){
      stop( print('  >>> Outfile not found - ', inraster))
    }
  }

  if (require('gdalUtils')){
    # print(1)

    gi <- gdalUtils::gdalinfo(inraster)

    nd0 <- as.numeric(gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)))
    nd <- length(nd0) == 1 & nd0 == -9999

    pixsize0 <- unlist(strsplit(split = ',', gsub('Pixel Size = \\(|\\)', '', grep('Pixel Size', gi, value = TRUE))))
    options(digits = max(nchar(pixsize0)))
    (pixsize <- abs(as.numeric(pixsize0 )))
    ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )

    if( !( nd & ps ) ) {
      gdalUtils::gdalwarp(srcfile = inraster, dstfile = outraster,
                          tr = rep(max(pixsize), 2), dstnodata = -9999)
      outraster0 <- outraster
    }

  } else if(require('gdalUtilities')){
    # print(2)

    # options(scipen = 999)
    # options(scipen = 9)

    gi <- capture.output(gdalUtilities::gdalinfo(inraster))

    nd0 <- as.numeric(gsub('^.+\\=', '', grep('NoData Value=', gi, value = TRUE)))
    nd <- length(nd0) == 1 & nd0 == -9999

    pixsize0 <- unlist(strsplit(split = ',', gsub('Pixel Size = \\(|\\)', '', grep('Pixel Size', gi, value = TRUE))))
    options(digits = max(nchar(pixsize0)))
    (pixsize <- abs(as.numeric(pixsize0 )))
    ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )

    if( !( nd & ps ) ) {
      gdalUtilities::gdalwarp(srcfile = inraster, dstfile = outraster,
                              tr = rep(max(pixsize), 2),
                              #                      srcnodata = nd0, dstnodata = nd0,
                              overwrite = TRUE)
      # gdalUtilities::gdalinfo(outraster)
      outraster0 <- outraster
    }

  } else if(require('raster')){
    # print(3)

    options(digits = 20)
    rx <- raster(inraster)
    nd <- rx@file@nodatavalue == -9999

    (pixsize <- res(rx))
    ps <- (length(pixsize) == 2 & pixsize[1]==pixsize[2] )
    if( !( nd & ps ) ) {
      templ <- raster(  crs = rx@crs, res = rep(max(res(rx)), 2), ext = extent(rx))
      raster::resample(x = rx, y = templ, filename = outraster, NAflag = -9999, overwrite = FALSE)
      outraster0 <- outraster
    }
  }
  options(digits = 5)
  return(outraster0)
}


loadRast <- function(inFiles, tempFolder, sessID){
  # inFile <- input$shapefile

}


# xx <- loadShp(inFiles, tempFolder, sessID, rastTemp = NULL)
loadShp <- function(inFiles, tempFolder, sessID, rastTemp = NULL){ # inFiles <- input$shapefile
  #inFiles: data frame with files. column 'name' with file name, and 'newFile' with full new name path

  ## setwd('/data/tempR/')
  # sessID <- 'colaJOU2025072917142505/'
  # tempFolder <- paste0('/data/tempR/', sessID)
  # inFiles <- data.frame(newFile = list.files(
  # full.names = TRUE, path = tempFolder, pattern = 'Anoa'))
  # inFiles$name <- basename(inFiles$newFile)
  # rv$inDistSessID
  # sessID <- (inDistSessID <- sessionIDgen())
  # tempFolder <- "/data/temp/SL2023112814374705file138b6c29cf34/"
  # save(inFiles, file = paste0(rootPath, '/inFiles_input_shp.RData'))
  # sss <- load(paste0('/data/tempR/inFiles_input_shp.RData'))
  # load(paste0(rootPath, '/inFiles_input_shp.RData'))
  # sss <- load(paste0(tempFolder, '/shpfiles.RData')) # sss

  #print(inFiles)
  outshp <- list()
  saveShp <- FALSE
  if ( class(inFiles) == "NULL" ){
    outshp$mssg <- 'Error in loading shapefile'
  } else {
    if ( class(inFiles) == 'data.frame' & nrow(inFiles) == 1){
      ## zip files
      if(grepl('*\\.zip', inFiles$name)){
        outZip <- paste0(tempFolder, '/shp_', sessID);
        dir.create(outZip)
        unzip(zipfile = inFiles$newFile, exdir = outZip)
        uZ <- list.files(outZip)
        #x0 <- uZ; print(x0); print(class(x0)); print(str(x0))
        outshp$shp <- tryCatch(sf::read_sf(outZip, layer = tools::file_path_sans_ext(uZ[1])), error = function (e) NULL)
        outshp$files <- uZ
        outshp$layer <- grep(pattern = '.shp', outshp$files, value = TRUE)
      } else if (grepl('\\.SQLite|\\.gpkg|\\.GeoJSON', inFiles$name)){ ## single
        #save(inFiles, file = 'inFileSingle.RData');
        outshp$shp <- tryCatch(sf::read_sf(inFiles$newFile[1]), error = function (e) NULL)
        outshp$files <- inFiles$newFile
        outshp$layer <- inFiles$newFile[1]
      } else if (grepl('\\.csv', inFiles$name)){ ## single CSV
        #save(inFiles, file = 'inFileSingle.RData'); load(paste0(rootPath, '/inFiles_input_shp.RData')) # sss
        # inFiles <- inFiles[5, ]
        outshp$shp <- tryCatch(sf::read_sf(inFiles$newFile[1]), error = function (e) NULL)
        outshp$files <- inFiles$newFile
        outshp$layer <- inFiles$newFile[1]
      } else if (grepl('\\.xy', inFiles$name)){ ## single
        #save(inFiles, file = 'inFileSingle.RData'); load(paste0(rootPath, '/inFiles_input_shp.RData')) # sss
        # inFiles <- inFiles[5, ]
        outshp$shp <- tryCatch(read.csv(inFiles$newFile[1]), error = function (e) NULL)
        outshp$shp$x1 <- outshp$shp[, grep('X|x', colnames(outshp$shp), value = TRUE)[1]]
        outshp$shp$y1 <- outshp$shp[, grep('Y|y', colnames(outshp$shp), value = TRUE)[1]]
        #coordinates(outshp$shp) =~ x1 + y1
        outshp$shp <- sf::st_as_sf(outshp$shp, coords = c("x1","y1"))
        outshp$files <- inFiles$newFile
        outshp$layer <- inFiles$newFile[1]
      }
    } else if ( nrow(inFiles) >= 3  & all(sapply(c('\\.shp', '\\.shx', '\\.dbf'), grep, inFiles$name)) ){ ## shp several

      #inFiles$datapath2 <- gsub('\\/[0-9]\\.', '/1.', inFiles$newFile)
      #sapply(inFiles$newFile, USE.NAMES = F, function(x){file.rename(x,  gsub('\\/[0-9]\\.', '/1.', x) ) })

      # setwd("/data/temp/PH2023100311442505file8513323368416")
      # outshp <- inFiles <- list(newFile = '/data/temp/KI2023100313381405file8513368894f1c/sp50_multi.shp')
      # inFiles <- data.frame(newFile = 'C:/cola/Anoa/Anoa_present_ardianti.shp')

      outshp$shp <- tryCatch(sf::read_sf(dirname(inFiles$newFile[1]),
                                         basename(tools::file_path_sans_ext(inFiles$newFile[1]))),
                             error = function  (e) {print(e); return(e)})
      # [1] "sf"         "tbl_df"     "tbl"        "data.frame"
      cat('\n Vectorial layer loaded. Class: ', paste0(class(outshp$shp), collapse = '-'), '\n')

    }

    ## Second try
    if( is.null(outshp$shp) ){
      outshp$shp <- tryCatch(sf::read_sf(
        grep('shp$', inFiles$newFile, value = TRUE)),
        error = function (e) {print(e); retunr(e)})
    }


    # v <- vect(f)
    #   s <- sf::st_as_sf(v) # sf object from a SpatVector
    #   vv <- vect(s) # SpatVector from an sf:
    #
    # library(raster)
    # x <- as(v, "Spatial") # To create an sp object from a SpatVector:
    # vs <- vect(x) # SpatVector from a Spatial* vector type object with


    #if worked
    if( !is.null(outshp$shp) ){

      # if (any(class(outshp$shp$geometry) == 'sfc_LINESTRING')){
      #   outshp$shp <- st_cast(outshp$shp, "LINE")
      #   outshp$shp <- as(outshp$shp, 'Spatial')
      # }

      #print( ' ............. '); print(outshp$shp)
      if ( any(class(outshp$shp$geometry) == 'sfc_MULTIPOINT') ){
        #cat (' Converting multipoint to singlepoint \n')
        outshp$shp <- st_cast(outshp$shp, "POINT")
        outshp$shp <- as(outshp$shp, 'Spatial')
        saveShp <- TRUE
      }

      if ( any(class(outshp$shp) %in% c('sf', 'SpatialPointsDataFrame')) ){
        #if (class(outshp$shp) == 'SpatialPointsDataFrame'){
        outshp$shp$sortID <- 1:nrow(outshp$shp)
      }

      if ( is.na(st_crs(outshp$shp)) ){
        # if (is.na(outshp$shp@proj4string@projargs)){
        #cat(' is.na(st_crs(outshp$shp)) \n')

        outshp$mssg <- 'No projection in shapefile'
      } else {
        #cat(' !is.na(st_crs(outshp$shp)) \n')

        # Project as raster
        if (! is.null(rastTemp) ){
          #cat(' ! is.null(rastTemp) \n')

          if ( class(rastTemp) %in% 'SpatRaster' ){
            #cat(' class(rastTemp) %in% SpatRaster \n')
            # rastTemp <- terra::rast('/data/tempR/colaJOU2025072917142505/in_crk_SHX2025072917144105.tif')
            # rastTemp <- terra::rast('C:/cola/Anoa/sur_colaMMR2025072905584605.tif')
            # outshp <- list(shp = sf::read_sf('C:/cola/Anoa/Anoa_present_ardianti.shp'),
            #                                  inFiles = 'C:/cola/Anoa/Anoa_present_ardianti.shp'
            # )

            #crsRst <- terra::crs(rastTemp)
            crsVct <- terra::vect(outshp$shp[1, 1])
            #cat(' Checking equal vector and raster CRS\n')
            #print(terra::same.crs(crsVct, rastTemp))
            if(! terra::same.crs(crsVct, rastTemp) ) {
              # cat(' Reprojecting points to raster CRS\n')
              # cat(' Old Projection: ', sf::st_crs(outshp$shp)$wkt,'\n')
              if ( any(class(outshp$shp) %in% c('sf')) ){
                outshp$shp <- st_transform(outshp$shp, crs = sf::st_crs(rastTemp))
              } else if ( any(class(outshp$shp) %in% c('SpatialPointsDataFrame')) ){
                outshp$shp <- terra::project(
                  terra::vect(outshp$shp),
                  rastTemp)
              }

              #shp2 <- terra::vect(outshp$shp) # SpatVector
              #outshp$shp <- terra::project(shp2, crs(rastTemp))
              # cat(' New Projection: ', sf::st_crs(outshp$shp)$wkt,'\n')
              saveShp <- TRUE
            }
          } else {

          }
        }
      }

      outshp$files <- inFiles$newFile
      outshp$layer <- grep(pattern = '.shp', outshp$files, value = TRUE)

      if (saveShp){
        # cat(' Vectorial layer loaded. Class: ', paste0(class(outshp$shp), collapse = '-'), '\n')
        # file.remove(inFiles$newFile)
        # inFiles <- data.frame(
        #   newFile = list.files(
        #   path = 'C:/cola/colaBKE2025073122594805/', pattern = 'Anoa')
        #   )
        outshp$files <- gsub(x = inFiles$newFile, '\\.', '_proj.')
        #print(outshp$files)
        #inFiles$newFile <- gsub(x = inFiles$newFile, '\\..', '_proj.')
        outshp$layer <- grep(pattern = '.shp', x = outshp$files, value = TRUE)
        #print(outshp$layer)
        if ( any(class(outshp$shp) %in% c('sf')) ){
          #cat (' Saving SF points into a new layer with raster projection \n',outshp$layer, '\n')
          sf::write_sf(outshp$shp,  overwrite=TRUE,
                       outshp$layer
                       #gsub(x = inFiles$newFile[1], '\\..+', '_proj.shp')
          )
        } else if ( any(class(outshp$shp) %in% c('SpatialPointsDataFrame')) ){
          #cat (' Saving TERRA points into a new layer with raster projection \n',outshp$layer, '\n')
          terra::writeVector(terra::vect(outshp$shp),  overwrite=TRUE,
                             outshp$layer
                             #gsub(x = inFiles$newFile[1], '\\..+', '_proj.shp')
          )
        } else if ( any(class(outshp$shp) %in% c('SpatVector')) ){
          terra::writeVector((outshp$shp),  overwrite=TRUE,
                             outshp$layer
                             #gsub(x = inFiles$newFile[1], '\\..+', '_proj.shp')
          )
        }
        #cat (' Writing projection: \n', sf::st_crs(outshp$shp)$wkt, '\n')
        outshp$shp <- tryCatch(sf::read_sf(
          outshp$layer
          #inFiles$newFile[1], '\\..+', '_proj.shp'
        ),
        error = function  (e) {print(e); return(e)})
      }

      # updateSelectInput(session, 'aoi_sur', choices = c('Dibujar', 'Capa'), selected = 'Capa')
      # pdebug("is.null(py)", 'inSurSessID', sep = '\n', pre = ' -- ')
      # try(file.remove(inFiles$name))
      # tempFolder <- 'C:/cola/colaJAA2025073022473805'
      save(outshp, file = paste0(tempFolder, '/outshp_loaded.RData'))
      # ss <- load(paste0(tempFolder, '/outshp_loaded.RData'))
      #cat('\n Vectorial layer loaded. Class: ', paste0(class(outshp$shp), collapse = '-'), '\n')
    }
  }
  #print(outshp)
  return(outshp)
}


pdebug <- function(devug, pre = '\n --\n', sep = '\n-',  ...){
  if (devug){
    x. = c(...)
    # x. = c('is.null(rv$newtifPath_dist)', 'rv$newtifPath_dist')
    # x. = c('grps'); grps = c('A', 'B', 'C')
    # print(x.)
    cat('\n', pre)
    invisible(sapply(x., function(x){
      # x = x.[2]
      y <- tryCatch(expr = eval(parse(text = x)), error = function(e) '-err-')
      if(is.null(y)){ y <- 'NULL' }
      tryCatch(cat(sep, x, ": ", y), error = function(e) e)
    }))
    cat('\n\t \n')
  }
} # pdebug("is.null(py)", 'inSurSessID', sep = '\n', pre = ' -- ')

# pdebug(devug=devug,sep='\n',pre='--', 'is.null(rv$newtifPath_dist)', 'rv$newtifPath_dist')
# pdebug(devug=devug,sep='\n',pre='--','print(inFiles)')

rowx <- function(...) {
  tags$div(class="row", ...)
}

colx <- function(width, ...) {
  tags$div(class=paste0("span", width), ...)
}


##function for removing the colum

# removecolumn <- function(df, nameofthecolumn){
#   df[ , -which(names(df) %in% nameofthecolumn)]
# }

removecolumn <- function(df, nameofthecolumn = NULL){
  id <- ifelse(is.null(nameofthecolumn),
               ncol(df),
               which(names(df) %in% nameofthecolumn))
  if(ncol(df) == 2){
    df[ , -id, drop = FALSE]
  } else if (ncol(df) == 1){
    df
  } else if (ncol(df) > 2){
    df[ , -id]
  }
}

addcolumn <- function(df, nameofthecolumn = NULL){
  id <- ifelse(is.null(nameofthecolumn), ncol(df), nameofthecolumn)
  cbind(df, df[ , id])
}






# for( x in 1:nrow(gd)){ # x = 1
#   case <- case0
#
#   xy <- gd[x, 1]
#   rs <- gd[x, 2]
#
#   (caseName <- paste0(outDir, '/case_', LETTERS[x], '.rip'))
#   (case[1, 1] <- gsub('small_test', paste0('case1', LETTERS[x]), case[1, 1]))
#   (case[2, 1] <- gsub('small_test.rsg', rs, case[2, 1]))
#   (case[3, 1] <- gsub('small_test_10pts.xy', xy, case[3, 1]))
#
#   (case[21, 1] <- "Save_Path_Output    TRUE") # Save_Path_Output
#   (case[24, 1] <- "Save_KDE_Output    TRUE") # Save_KDE_Output
#   (case[25, 1] <- "Save_Category_Output    TRUE") # Save_Category_Output
#   (case[26, 1] <- "Save_CDmatrix_Output    TRUE") # Save_CDmatrix_Output
#
#   if(file.exists(xy) & file.exists(rs)){
#     write.table(case, file = caseName, col.names = FALSE, row.names = FALSE, quote = FALSE)
#   }
# }
#
# outConf <- expand.grid(o1 = c(FALSE, TRUE), o2 = c(FALSE, TRUE), o3 = c(FALSE, TRUE), o4 = c(FALSE, TRUE))
# colnames(outConf) <- c('Save_Path_Output', 'KDE_Output', 'Category_Output', 'CDmatrix_Output')
#
# write.csv(outConf, 'caseO_output_config_time.csv')
#
# (x <- which(gd$case == 'R'))
# xy <- gd[x, 1]
# rs <- gd[x, 2]
# for( h in 1:nrow(outConf)){ # h = 1
#
#   case <- case0
#   (caseName <- paste0(outDir, '/case',LETTERS[x],'_', h,'.rip'))
#   print(caseName)
#   (case[1, 1] <- gsub('small_test', paste0( 'case',LETTERS[x], '_', h ), case[1, 1]))
#   (case[2, 1] <- gsub('small_test.rsg', rs, case[2, 1]))
#   (case[3, 1] <- gsub('small_test_10pts.xy', xy, case[3, 1]))
#
#   (case[21, 1] <- paste0("Save_Path_Output    ", outConf[h, 1])) # Save_Path_Output
#   (case[24, 1] <- paste0("Save_KDE_Output    ", outConf[h, 2])) # Save_KDE_Output
#   (case[25, 1] <- paste0("Save_Category_Output    ", outConf[h, 3])) # Save_Category_Output
#   (case[26, 1] <- paste0("Save_CDmatrix_Output    ", outConf[h, 4])) # Save_CDmatrix_Output
#
#   if(file.exists(xy) & file.exists(rs) ){
#     write.table(case, file = caseName, col.names = FALSE, row.names = FALSE, quote = FALSE)
#   }
# }


# cd %userprofile%\dockerdata\UNICOR\unicor
# python UNICOR.py caseA.rip

#for %a in (4 5 6 7 8 9 10 11 12 13 14 15) do python UNICOR.py caseR_%a.rip
#for /l %a in (1, 1, 16) do python UNICOR.py caseR_%a.rip

# >python UNICOR.py sabah.rip
# Log output directed to     : (sabah_test.log)
# The path list is empty, try increasing the path threshold, caclulating path as Nan
# Total UNICOR program run-time: 1:03:13.909943

#
# #### Plot results
# result <- read.csv('N:/Mi unidad/connectivity-nasa/01_original-data/Sabah_example_CDPOP+UNICOR/cases_config_time_A-X.csv')[1:24, ]
# result$nXWpts <- as.numeric(gsub('[a-zA-Z]|\\.|_', '', result$xy))
# result$npix <- gsub('size.+px_|tot.+', '', result$rsg)
#
# library(ggplot2)
# ggplot(result, aes(x = nXWpts, y = time, group = npix, color = npix)) +
#   geom_point() + geom_line() +
#   labs(title = 'Elapsed computing time in seconds',
#        x = 'Number of XW points', y = 'Time in seconds', color = 'Raster size in tot pixels') +
#   theme(legend.position="bottom")
