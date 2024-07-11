library(cola)
library(terra)

#wd <- '/home/shiny/cola/inst/examples'
wd <- '/home/user/cola_examples'
dir.create(wd)
setwd(wd)


hs <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
hsRast <- rast(hs)
nd <- guessNoData(hs)
plot(terra::rast(hs), main = 'Habitat suitability')


## S2R ---------------------
## Param 3 -- min val
srp30 <- s2res_py(intif = hs,
                  outtif = paste0(wd, '/hs_p3_0.tif'),
                  param3 = 0, # min val
                  param4 = 1, # max val
                  param5 = 10,
                  param6 = 1,
                  param7 = NULL, param8 = 'None')

srp305 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p3_05.tif'),
                   param3 = 0.5, # min val
                   param4 = 1, # max val
                   param5 = 10, # max output val
                   param6 = 1, # Shape
                   param7 = NULL, param8 = 'None')

srp309 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p3_05.tif'),
                   param3 = 0.9, # min val
                   param4 = 1, # max val
                   param5 = 10, # max output val
                   param6 = 1, # Shape
                   param7 = NULL, param8 = 'None')

plot(rast(srp30$file), main = 'Surface resistance. Min. val = 0')
plot(rast(srp305$file), main = 'Surface resistance. Min. val = 0.5')
plot(rast(srp309$file), main = 'Surface resistance. Min. val = 0.9')

range(rast(srp30$file))
range(rast(srp305$file))
range(rast(srp309$file))

## Param 4 -- max val
srp40 <- s2res_py(intif = hs,
                  outtif = paste0(wd, '/hs_p4_1.tif'),
                  param3 = 0, # min val
                  param4 = 1, # max val
                  param5 = 10,
                  param6 = 1,
                  param7 = NULL, param8 = 'None')

srp409 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p4_09.tif'),
                   param3 = 0, # min val
                   param4 = 0.9, # max val
                   param5 = 10, # max output val
                   param6 = 1, # Shape
                   param7 = NULL, param8 = 'None')

srp405 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p4_05.tif'),
                   param3 = 0, # min val
                   param4 = 0.5, # max val
                   param5 = 10, # max output val
                   param6 = 1, # Shape
                   param7 = NULL, param8 = 'None')

plot(rast(srp40$file), main = 'Surface resistance. Max. val = 1')
plot(rast(srp409$file), main = 'Surface resistance. Max. val = 0.9')
plot(rast(srp405$file), main = 'Surface resistance. Max. val = 0.5')

plot( rast(srp40$file) - rast(srp409$file) )



## Param 5 -- max val
srp510 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p5_10.tif'),
                   param3 = 0, # min val
                   param4 = 1, # max val
                   param5 = 10,
                   param6 = 1,
                   param7 = NULL, param8 = 'None')

srp550 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p5_50.tif'),
                   param3 = 0, # min val
                   param4 = 1, # max val
                   param5 = 50, # max output val
                   param6 = 1, # Shape
                   param7 = NULL, param8 = 'None')

srp5100 <- s2res_py(intif = hs,
                    outtif = paste0(wd, '/hs_p5_100.tif'),
                    param3 = 0, # min val
                    param4 = 1, # max val
                    param5 = 100, # max output val
                    param6 = 1, # Shape
                    param7 = NULL, param8 = 'None')

plot(rast(srp510$file), main = 'Surface resistance. Max. ouput val = 10')
plot(rast(srp550$file), main = 'Surface resistance. Max. ouput val = 50')
plot(rast(srp5100$file), main = 'Surface resistance. Max. ouput val = 100')

## Param 6 -- max out val
srp60 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p6_0.tif'),
                   param3 = 0, # min val
                   param4 = 1, # max val
                   param5 = 100,
                   param6 = 0,
                   param7 = NULL, param8 = 'None')


srp6m1 <- s2res_py(intif = hs,
                  outtif = paste0(wd, '/hs_p6_m1.tif'),
                  param3 = 0, # min val
                  param4 = 1, # max val
                  param5 = 100,
                  param6 = -1,
                  param7 = NULL, param8 = 'None')


srp61 <- s2res_py(intif = hs,
                  outtif = paste0(wd, '/hs_p6_1.tif'),
                  param3 = 0, # min val
                  param4 = 1, # max val
                  param5 = 100,
                  param6 = 1,
                  param7 = NULL, param8 = 'None')

srp6m5 <- s2res_py(intif = hs,
                    outtif = paste0(wd, '/hs_p6_m5.tif'),
                    param3 = 0, # min val
                    param4 = 1, # max val
                    param5 = 100,
                    param6 = -5,
                    param7 = NULL, param8 = 'None')

srp65 <- s2res_py(intif = hs,
                   outtif = paste0(wd, '/hs_p6_5.tif'),
                   param3 = 0, # min val
                   param4 = 1, # max val
                   param5 = 100,
                   param6 = 5,
                   param7 = NULL, param8 = 'None')

srp6m10 <- s2res_py(intif = hs,
                    outtif = paste0(wd, '/hs_p6_m10.tif'),
                    param3 = 0, # min val
                    param4 = 1, # max val
                    param5 = 100,
                    param6 = -10,
                    param7 = NULL, param8 = 'None')

srp610 <- s2res_py(intif = hs, outtif = paste0(wd, '/hs_p6_10.tif'),
                   param3 = 0, param4 = 1, param5 = 100,
                   param7 = NULL, param8 = 'None',
                   param6 = 10)


## plot results

breakss <- seq(0, 100, by = 10)
colx <- rev(terrain.colors(length(breaks)))

plot(rast(srp6m1$file), main = 'Surface resistance. Shape val = -1', breaks = breakss, col = colx)
plot(rast(srp61$file), main = 'Surface resistance. Shape val = 1', breaks = breakss, col = colx)
plot(rast(srp6m5$file), main = 'Surface resistance. Shape val = -5', breaks = breakss, col = colx)
plot(rast(srp65$file), main = 'Surface resistance. Shape val = 5', breaks = breakss, col = colx)
plot(rast(srp6m10$file), main = 'Surface resistance. Shape val = -10', breaks = breakss, col = colx)
plot(rast(srp610$file), main = 'Surface resistance. Shape val = 10', breaks = breakss, col = colx)


plot(rast(srp6m1$file), main = 'Surface resistance. Shape val = -1')
plot(rast(srp6m1$file), main = 'Surface resistance. Shape val = -1', breaks = breakss, col = colx)

plot(rast(srp6m1$file), breaks = breakss, col = colx)

plot(hsRast[], rast(srp60$file)[], pch = 20, ylim = c(0, 100))
points(hsRast[], rast(srp6m1$file)[], pch = 20, col = 2)
points(hsRast[], rast(srp61$file)[], pch = 20, col = 3)
points(hsRast[], rast(srp6m5$file)[], pch = 20, col = 4)
points(hsRast[], rast(srp65$file)[], pch = 20, col = 5)
points(hsRast[], rast(srp6m10$file)[], pch = 20, col = 6)
points(hsRast[], rast(srp610$file)[], pch = 20, col = 7)



## Points ---------------------


param3 = '20' # Max rast value
param4 =  95 # Min rast value
param5 = 50 # Number of points

## Change min - max
pts_20_95_50 <- points_py(intif = hs, outshp = '/home/shiny/cola/inst/sampledata/points_sabahA_50.shp',
                          param3 = 0.20, param4 = 0.95, param5 = 50, param6 = 'None')

pts_20_95_100 <- points_py(intif = hs, outshp = '/home/shiny/cola/inst/sampledata/points_sabahA_100.shp',
                          param3 = 0.20, param4 = 0.95, param5 = 100, param6 = 'None')

## Change number of points
pts_10_100_100 <- points_py(intif = hs, outshp = '/home/shiny/cola/inst/sampledata/points_sabahB_100.shp',
                            param3 = 0.10, param4 = 0.5, param5 = 50, param6 = 'None')


## Plot results
hsRast <- rast(hs)
exA1 <- extract(hsRast, vect(pts_20_95_50$file))
range(exA1$sampleTif)
nrow(exA1)
exA2 <- extract(hsRast, vect(pts_20_95_100$file))
range(exA2$sampleTif)
nrow(exA2)
exB1 <- extract(hsRast, vect(pts_10_100_100$file))
range(exB1$sampleTif)
nrow(exB1)


plot(hsRast, main = '50 points [0.2 - 0.9]')
points(vect(pts_20_95_50$file))
plot(hsRast, main = '100 points [0.2 - 0.9]')
points(vect(pts_20_95_100$file))
plot(hsRast, main = '50 points [0.1 - 0.5]')
points(vect(pts_10_100_100$file))


## Kernels ---------------------

# [4] distance threshold (in cost distance units)
# [5] kernel shape (linear, gaussian)
# [6] kernel volume

pts_path <- '/home/shiny/cola/inst/sampledata/points_sabahA_50.shp'
pts_vect <- vect(pts_path)

## Param 4 -- distance threshold

crk_p4A <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p4A.tif',
                   param4 = 125000, param5 = 'linear', param6 = 1)

crk_p4B <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p4B.tif',
                   param4 = 100000, param5 = 'linear', param6 = 1)

crk_p4C <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p4C.tif',
                   param4 = 50000, param5 = 'linear', param6 = 1)

plot(rast(crk_p4A$file), main = 'Corridor tolerance: 125k')
points(pts_vect, col = '1', cex = .4)

plot(rast(crk_p4B$file), main = 'Corridor tolerance: 100k')
points(pts_vect, col = '1', cex = .4)

plot(rast(crk_p4C$file), main = 'Corridor tolerance: 50k')
points(pts_vect, col = '1', cex = .4)


## Param 5 -- kernel shape (linear, gaussian)

crk_p5A <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p5A.tif',
                   param4 = 20000, param5 = 'linear', param6 = 1)

crk_p5B <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p5B.tif',
                   param4 = 20000, param5 = 'gaussian', param6 = 1)

plot(rast(crk_p5A$file), main = 'Kernel shape: linear')
points(pts_vect, col = '1', cex = .4)

plot(rast(crk_p5B$file), main = 'Kernel shape: gausian')
points(pts_vect, col = '1', cex = .4)

## Param 6 -- kernel volume

crk_p6A <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p6A.tif',
                   param4 = 20000, param5 = 'linear', param6 = 1)

crk_p6B <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p6B.tif',
                   param4 = 20000, param5 = 'linear', param6 = 2)

crk_p6C <-  crk_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/crk_p6C.tif',
                   param4 = 20000, param5 = 'linear', param6 = 10)

plot(rast(crk_p6A$file), main = 'kernel volume: 1')
points(pts_vect, col = '1', cex = .4)

plot(rast(crk_p6B$file), main = 'kernel volume: 5')
points(pts_vect, col = '1', cex = .4)

plot(rast(crk_p6C$file), main = 'kernel volume: 10')
points(pts_vect, col = '1', cex = .4)




## Corridors ---------------------

# Param [4] distance threshold (should be in meters*)
# Param [5] corridor smoothing factor (in number of cells)
# Param [6] corridor tolerance (in cost distance units)

pts_path <- '/home/shiny/cola/inst/sampledata/points_sabahA_50.shp'
pts_path <- '/home/user/cola/inst/sampledata/points_sabahA_50.shp'

hs_path <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
points_path <- system.file(package = 'cola', 'sampledata/samplePoints.shp')

## Param 4 -- distance threshold

p4_150k <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p4150k.tif',
                  param4 = 150000, param5 = 0, param6 = 0)

p4_100k <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p4100k.tif',
                  param4 = 100000, param5 = 0, param6 = 0)

p4_50k <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p450k.tif',
                 param4 = 50000, param5 = 0, param6 = 0)

plot(rast(p4_150k$file), main = 'LCC: 150k')
points(pts_vect, col = '1', cex = .4)

plot(rast(p4_100k$file), main = 'LCC: 100k')
points(pts_vect, col = '1', cex = .4)

plot(rast(p4_50k$file), main = 'LCC: 50k')
points(pts_vect, col = '1', cex = .4)


## Param 5 -- corridor smoothing factor

p5_0 <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p0.tif',
               param4 = 150000, param5 = 0, param6 = 0)

p5_2 <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p51.tif',
               param4 = 150000, param5 = 1, param6 = 0)

p5_10 <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p510.tif',
                param4 = 150000, param5 = 10, param6 = 0)

plot(rast(p5_0$file), main = 'LCC: 0')
plot(rast(p5_2$file), main = 'LCC: 1')
plot(rast(p5_10$file), main = 'LCC: 10')


## Param 6 -- corridor tolerance

p6_0 <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p60.tif',
               param4 = 150000, param5 = 1, param6 = 0)

p6_b <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p6b.tif',
               param4 = 150000, param5 = 1, param6 = 5)

p6_c <- lcc_py(inshp = pts_path, intif = hs, outtif = '/home/shiny/cola/inst/sampledata/lcc_p6c.tif',
                param4 = 150000, param5 = 1, param6 = 100)

plot(rast(p6_0$file), main = 'Corridor tolerance: 0')
plot(rast(p6_b$file), main = 'Corridor tolerance: 5')
plot(rast(p6_c$file), main = 'Corridor tolerance: 100')



## Prioritization  ---------------------

pri_py(tif, incrk, inlcc,
       maskedcsname = paste0(tempfile(), '.tif'),
       outshppoint, outshppol, outshppatch,
       outtifpatch, outtif,
       param7 = 0.5,
       param8 = 1000)
