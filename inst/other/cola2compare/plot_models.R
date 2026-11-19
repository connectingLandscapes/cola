library(broman)
library(dismo)
library(ggplot2)
library(gridExtra)
library(maptools)
library(raster)
library(rasterVis)
library(rgdal)

setwd('X:/CI')

(r <- list.files(pattern = '.+tif'))
df <- data.frame(t(sapply(tools::file_path_sans_ext(r), function(x){
  strsplit(x, '_')[[1]]
})))


head(df)
table(df$X1, df$X2)
table(df$X1, df$X2, df$X3)

df2 <- subset(df, !grepl(pattern = 'obs', X2))
table(df2$X1, df2$X2)

head(df)
subset(df, grepl(pattern = 'bs', X2))

outDir <- 'X:/Google Drive/CI/Entregables_Retamo/05_Modelos-distribucion/'
load('X:/Google Drive/CI/occ/occ602_final.RData')
my.at <- seq(0, 1000, 100)


## Full Ulex ----

png(filename = paste0(outDir, 'Ulex_obs_full.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Ulex_obs_Full_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()


## Full Genista ----

png(filename = paste0(outDir, 'Genista_obs_full.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Genista_obs_Full_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()




## Small Ulex ----

png(filename = paste0(outDir, 'Ulex_obs_full.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Ulex_obs_Small_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()


## Full Ulex ----

png(filename = paste0(outDir, 'Genista_obs_full.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Genista_obs_Small_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()




## PT Ulex ----

png(filename = paste0(outDir, 'Ulex_obs_TP.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Ulex_obs_TP_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()


##  PT Genista ----

png(filename = paste0(outDir, 'Genista_obs_TP.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Genista_obs_TP_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p12, p13, p14, p15, p16, ncol=3)
dev.off()



## P Ulex ----

png(filename = paste0(outDir, 'Ulex_obs_P.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Ulex_obs_P_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
# p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
#                  margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p13, p14, p15, p16, ncol=3)
dev.off()


## P Genista  ----

png(filename = paste0(outDir, 'Genista_obs_P.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Genista_obs_P_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
# p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
#                  margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11,  p13, p14, p15, p16, ncol=3)
dev.off()


## T Ulex ----

png(filename = paste0(outDir, 'Ulex_obs_T.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Ulex_obs_T_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
# p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
#                  margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11, p13, p14, p15, p16, ncol=3)
dev.off()



## T Genista  ----

png(filename = paste0(outDir, 'Genista_obs_T.png'),
    width = 739, height = 629, units = 'px', pointsize = 1/100) #par(mfrow = c(1, 6), bg = NA)
prefix <- 'Genista_obs_T_'
p11 <- levelplot(raster(paste0(prefix, 'ANN.tif')), main = paste0('ANN'),  at = my.at,# colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
# p12 <- levelplot(raster(paste0(prefix, 'GBM.tif')), main = paste0('GBM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
#                  margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p13 <- levelplot(raster(paste0(prefix, 'GLM.tif')), main = paste0('GLM'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90)))
p14 <- levelplot(raster(paste0(prefix, 'MAX.tif')), main = paste0('MAX'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p15 <- levelplot(raster(paste0(prefix, 'RF.tif')), main = paste0('RF'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
p16 <- levelplot(raster(paste0(prefix, 'SRE.tif')), main = paste0('SRE'), at = my.at,# colorkey=list(space="bottom"),colorkey=list(space="bottom"),
                 margin=FALSE, xlab = '', ylab = '',  scales=list(y=list(rot=90))) #+
grid.arrange(p11,  p13, p14, p15, p16, ncol=3)
dev.off()


############### Scenarios ---



df2 <- subset(df, !grepl(pattern = 'obs', X2))
table(df2$X1, df2$X2)

unique(df$X4)
head(df)
subset(df, grepl(pattern = 'bs', X2))
subset(df, grepl(pattern = 'SRE', X4))

outDir <- 'X:/Google Drive/CI/Entregables_Retamo/05_Modelos-distribucion/'
my.at <- seq(0, 1000, 100)

U_ANN <- stack(c('Ulex_obs_Small_ANN.tif', 'Ulex_P1RCP26_Small_ANN.tif', 'Ulex_P3RCP85_Small_ANN.tif'))
U_GBM <- stack(c('Ulex_obs_Small_GBM.tif', 'Ulex_P1RCP26_Small_GBM.tif', 'Ulex_P3RCP85_Small_GBM.tif'))
U_GLM <- stack(c('Ulex_obs_Small_GLM.tif', 'Ulex_P1RCP26_Small_GLM.tif', 'Ulex_P3RCP85_Small_GLM.tif'))
U_MAX <- stack(c('Ulex_obs_Small_MAX.tif', 'Ulex_P1RCP26_Small_MAX.tif', 'Ulex_P3RCP85_Small_MAX.tif'))
U_RF <- stack(c('Ulex_obs_Small_RF.tif', 'Ulex_P1RCP26_Small_RF.tif', 'Ulex_P3RCP85_Small_RF.tif'))
U_SRE <- stack(c('Ulex_obs_Small_SRE.tif', 'Ulex_P1RCP26_Small_SRE.tif', 'Ulex_P3RCP85_Small_SRE.tif'))


levelplot(U_ANN)
levelplot(U_GBM) # NO
levelplot(U_GLM)
levelplot(U_MAX)
levelplot(U_RF)
levelplot(U_SRE)



G_ANN <- stack(c('Genista_obs_Small_ANN.tif', 'Geni_P1RCP26_Small_ANN.tif', 'Geni_P3RCP85_Small_ANN.tif'))
G_GBM <- stack(c('Genista_obs_Small_GBM.tif', 'Geni_P1RCP26_Small_GBM.tif', 'Geni_P3RCP85_Small_GBM.tif'))
G_GLM <- stack(c('Genista_obs_Small_GLM.tif', 'Geni_P1RCP26_Small_GLM.tif', 'Geni_P3RCP85_Small_GLM.tif'))
G_MAX <- stack(c('Genista_obs_Small_MAX.tif', 'Geni_P1RCP26_Small_MAX.tif', 'Geni_P3RCP85_Small_MAX.tif'))
G_RF <- stack(c('Genista_obs_Small_RF.tif', 'Geni_P1RCP26_Small_RF.tif', 'Geni_P3RCP85_Small_RF.tif'))
G_SRE <- stack(c('Genista_obs_Small_SRE.tif', 'Geni_P1RCP26_Small_SRE.tif', 'Geni_P3RCP85_Small_SRE.tif'))


levelplot(G_ANN)
levelplot(G_GBM) # NO
levelplot(G_GLM)
levelplot(G_MAX)
levelplot(G_RF)
levelplot(G_SRE)

#Genista monspessulana


