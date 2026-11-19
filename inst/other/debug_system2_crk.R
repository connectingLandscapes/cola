if(FALSE){
  
  inshp = 'C:/cola/andy/out_simpts_QYI2026050114545005.shp'
  intif = 'C:/cola/andy/scl1000m.tif'
  outtif = 'C:/cola/andy/out_crk.tif'
  maxdist = 10000
  shape = 'linear'
  transform = 'no'
  volume = 1
  ncores = as.numeric(Sys.getenv('COLA_NCORES'))
  crs = 'None'
  maxram = as.numeric(Sys.getenv('COLA_RAMGB'))
  tempFolder = NULL
  py = Sys.getenv("COLA_PYTHON_PATH")
  pyscript = system.file(package = 'cola', 'python/crk_joblib.py')
  cml = TRUE; show.result = TRUE
  
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
    cat('\n\n')
  }
  
  intCMD <- tryCatch(system(cmd_crk, intern = TRUE), error = function(e) e$message)
  # intCMD2 <- tryCatch(system2(cmd_crk, stdout = TRUE), error = function(e) e$message)
  # output <- system2(py, args = c (
  #   quotepath(pyscript),
  #   quotepath(inshp),
  #   quotepath(intif), quotepath(outtif),
  #   format(maxdist, scientific=F),
  #   shape, transform,
  #   format(volume, scientific=F),
  #   format(ncores, scientific=F),
  #   crs, h5file, maxram), stdout = TRUE)
  
  
  if(show.result){
    print(intCMD)
  }
  
  logname <- paste0(tools::file_path_sans_ext(outtif), '.metadata')
  metaFile <- c(inshp = inshp, intif = intif, outtif = outtif, maxdist = maxdist,
                shape = shape, transform = transform, volume = volume,
                ncores = ncores, crs = crs,
                maxram = maxram,
                log = paste0(intCMD, collapse = ' - '),
                done = ifelse(file.exists(outtif), 'yes', 'no'))
  write.table(metaFile, logname )
  
  
  tryCatch(file.remove(c(h5file, h5file2)), error = function(e) NULL)
  
  ans <- list(file = ifelse(file.exists(outtif), outtif, ''),
              # log =  paste0(intCMD, ' -- ', read.delim(logname)) ) )
              log = paste0("", intCMD) )
}