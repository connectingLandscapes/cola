if(FALSE){
  julia --threads 4
  
  Threads.nthreads()
  using Omniscape
  cd("/home/shiny/paper/ini+scripts/")
  readdir()
  file = "omni10.ini"
  content = read(file, String)
  println(content)
  run_omniscape(file)
  
  # [Required]
  # resistance_file = /home/shiny/paper/inputs/size10_85m_24Mpix.tif
  # radius = 5429
  # block_size = 1
  # project_name = /home/shiny/paper/output/test10_5429
  # 
  # [General options]
  # source_from_resistance = false
  # source_file = /home/shiny/paper/inputs/points100_as_size10.tif
  # r_cutoff = 1
  # calc_normalized_current = true
  # 
  # parallelize = true
  # parallel_batch_size = 1
  # 
  # [Output options]
  # write_raw_currmap = true
}

setwd('/home/shiny/paper/inputs/')
list.files()
library(terra)
pts <- vect('points_sabah_50.shp')
r <- rast('size5.tif')
rp <- terra::rasterize(pts, r, background  = 0)
writeRaster(rp, 'points50_as_size5.tif')
plot(rp)
plot(r)
minmax(rp)
plot(pts, add = TRUE)

plot(rast('../output/test5_1/cum_currmap.tif'))
plot(rast('../output/test7_2465/cum_currmap.tif'), main = '5 2465')
plot(rast('../output/test7_2465/normalized_cum_currmap.tif'), main = '5 2465')

par(mfrow = c(2, 2))
stk5 <- rast(c(paste0('../output/test5_',c(100, 200, 390, 500) ,'/normalized_cum_currmap.tif')))
global(stk5, 'max')
setMinMax(stk5)

plot(rast('../output/test5_100/normalized_cum_currmap.tif'), main = '5 100')
plot(rast('../output/test5_200/normalized_cum_currmap.tif'), main = '5 200')
plot(rast('../output/test5_390/normalized_cum_currmap.tif'), main = '5 390')
plot(rast('../output/test5_500/normalized_cum_currmap.tif'), main = '5 500')

plot(rast('../output/test7_2465_1/normalized_cum_currmap.tif'), main = '7')
plot(rast('../output/test6_1232/normalized_cum_currmap.tif'), main = '6')
plot(rast('../output/test5_390/normalized_cum_currmap.tif'), main = '5')
plot(rast('../output/'), main = '5')



gdalUtilities::gdalinfo('size2.tif') # 12
gdalUtilities::gdalinfo('size3.tif') # 38
gdalUtilities::gdalinfo('size4.tif') # 123
gdalUtilities::gdalinfo('size5.tif') # 389
gdalUtilities::gdalinfo('size6.tif') # 1232
gdalUtilities::gdalinfo('size7.tif') # 2465
gdalUtilities::gdalinfo('size8_170m_6Mpix.tif') # 2715
gdalUtilities::gdalinfo('size9.tif') # 3846 * 3121 12M
gdalUtilities::gdalinfo('size10_85m_24Mpix.tif') # 5429, 4406



library(cola)
inshp <- '/home/shiny/paper/inputs/points_sabah_100.shp'

system.time(
  crk6 <- cola::crkJoblib_py(inshp = inshp, intif = '/home/shiny/paper/inputs/size6.tif', 
                             outtif = '/home/shiny/paper/output/cola_crk6.tif', volume = 0,
                             maxdist = 10000, shape = 'linear', ncores = 6, maxram = 50)
)
# plot(rast('/home/shiny/paper/output/cola_crk6.tif'))
# user  system elapsed 
# 38.926   3.960  24.742

system.time(
  crk7 <- cola::crkJoblib_py(inshp = inshp, intif = '/home/shiny/paper/inputs/size7.tif', 
                             outtif = '/home/shiny/paper/output/cola_crk7.tif', volume = 0,
                             maxdist = 10000, shape = 'linear', ncores = 6, maxram = 50)
)

# user  system elapsed 
# 284.355  27.610  84.834
system.time(
  crk8 <- cola::crkJoblib_py(inshp = inshp, intif = '/home/shiny/paper/inputs/size8_170m_6Mpix.tif', 
                             outtif = '/home/shiny/paper/output/cola_crk8.tif', volume = 0,
                             maxdist = 10000, shape = 'linear', ncores = 6, maxram = 50)
)

# user  system elapsed 
# 709.854  39.879 177.393 
system.time(
  crk9 <- cola::crkJoblib_py(inshp = inshp, intif = '/home/shiny/paper/inputs/size9.tif', 
                              outtif = '/home/shiny/paper/output/cola_crk9.tif', volume = 0,
                              maxdist = 10000, shape = 'linear', ncores = 6, maxram = 50)
)

# user  system elapsed 
# 609.930  35.551 153.317 

system.time(
  crk10 <- cola::crkJoblib_py(inshp = inshp, intif = '/home/shiny/paper/inputs/size10_85m_24Mpix.tif', 
                              outtif = '/home/shiny/paper/output/cola_crk10.tif', volume = 0,
                              maxdist = 10000, shape = 'linear', ncores = 6, maxram = 50)
)
# user   system  elapsed 
# 1486.027   38.212  332.557



library(terra)
plot(rast(crk10$file))
plot(rast('/home/shiny/paper/output/test9_3846_2/normalized_cum_currmap.tif'))

if (FALSE){
  
##
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python /home/shiny/cola/inst/UNICOR/UNICOR.py /home/shiny/paper/ini+scripts/crk0.rip &> /home/shiny/paper/logUNcrk.txt
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk0.ini

time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk0.ini

sudo su
cd /home/shiny/cola/inst/UNICOR/unicor
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk8_50.ini &> /home/shiny/paper/logcrk_850.txt
#2m10
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk8_100.ini &> /home/shiny/paper/logcrk_8100.txt
#2m24
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk9_50.ini &> /home/shiny/paper/logcrk_950.txt
# 4m09s
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk9_100.ini &> /home/shiny/paper/logcrk_9100.txt
# 4m54
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk10_50.ini 
# 8m38s
time /home/shiny/.local/share/r-miniconda/envs/cola/bin/python UNICOR.py /home/shiny/paper/ini+scripts/crk10_100.ini 
# 9m34s
pwd
pwd

plot(rast('/size8.asc_pts50.kdepaths'))
plot(rast('/size8.asc_pts50.addedpaths.txt'))

# Session_label    resistant_kernel
# Grid_Filename    /home/shiny/paper/inputs/size9.asc
# XY_Filename    /home/shiny/paper/inputs/pts100.xy
# Use_Direction	FALSE
# Type_Direction	FlowAcc
# Use_Resistance	TRUE
# Barrier_or_U_Filename	
# Direction_or_V_Filename	
# Speed_To_Resistance_Scale	
# Use_ED_threshold    False
# ED_Distance   10000
# Edge_Type	all_paths
# Transform_function	linear
# Const_kernal_vol	True
# Kernel_volume   10000
# Edge_Distance    10000
# Number_of_Processes    6
# KDE_Function    Gaussian
# KDE_GridSize    2
# Number_of_Categories    5
# Save_Path_Output    TRUE
# Save_IndividualPaths_Output    FALSE
# Save_GraphMetrics_Output    FALSE
# Save_KDE_Output    FALSE
# Save_Category_Output    FALSE
# Save_CDmatrix_Output    FALSE
}