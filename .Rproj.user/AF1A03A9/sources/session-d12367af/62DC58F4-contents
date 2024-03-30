## Cola utils

# library(dismo)
# library(raster)
# 
# setwd('N:/My Drive/connectivity-nasa-USFSIP/01_original-data/Sabah_example_CDPOP+UNICOR')
# r <- raster('N:/My Drive/connectivity-nasa-USFSIP/01_original-data/Sabah_example_CDPOP+UNICOR/roads.tif')
# outDir <- 'C:/Users/Admin/dockerdata/UNICOR/unicor/'


# options(scipen = 999)
# options(scipen = 9)

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
    
    (outdir <- ifelse(is.null(outdir), basename(path), gsub(pattern = basename(outdir), x = outdir, replacement = '') ))
    fname <- paste0(outdir, '/', tools::file_path_sans_ext( basename(path)), '.asc')
    fname_rsg <- paste0(outdir, '/', tools::file_path_sans_ext( basename(path)), '.rsg')
    
    writeRaster(x = tif, format = 'ascii', overwrite = TRUE, filename = fname )
    
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




fitRaster2cola0 <- function(inrasterpath, outrasterpath = NULL){
  # setwd('N:/Mi unidad/git/connecting-landscapes/performance-tests/inputs')
  # inrasterpath = 'orig_tifs/size6.tif'
  # outrasterpath = 'size6.tif'
  
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
                              tr = rep(max(pixsize), 2), dstnodata = -9999)
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

fitRaster2cola <- function(inrasterpath, outrasterpath = NULL){
  # setwd('N:/Mi unidad/git/connecting-landscapes/performance-tests/inputs')
  # inrasterpath = 'orig_tifs/size6.tif'
  # outrasterpath = 'size6.tif'
  
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
    } else{
      outraster0 <- inraster
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
                              tr = rep(max(pixsize), 2), dstnodata = -9999)
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
      raster::resample(x = rx, y = templ, filename = outraster, NAflag = -9999, overwrite = FALSE)
      outraster0 <- outraster
    } else{
      outraster0 <- inraster
    }
  }
  
  options(digits = 5)
  return(outraster0)
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
