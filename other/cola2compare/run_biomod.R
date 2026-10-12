#install.packages('biomod2')
#install.packages('rJava')
library(biomod2)
library(dismo)
library(ggplot2)
library('rJava')
library(rgdal)
maxent()

outDir <- 'X:/Google Drive/CI/Entregables_Retamo/'

##
aoi <- readOGR('H:/CI/Corredor_gcs.shp')
pnn <- spTransform(readOGR('H:/CI/pnn.shp'), CRSobj = aoi@proj4string)
load('C:/GoogleDrive/CI/Entregables_Retamo/03_Registros-especies/Registros_especies_puntos/occ602_final.RData')
dim(occ)
dim(rrbb)
table(rrbb$species, useNA = 'alw')

coords <- coordinates(occ)
occ$x_ <- occ$x
occ$y_ <- occ$y




# trainPos <- tapply(rrbb$pa, rrbb$species, function(x){
#   n <- round(length(x) * prTrain)
#   sample(length(x), n)
# })
##
# occ$train <- rrbb$train <- 0
# occ$train[occ$species == 'Genista monspessulana'][trainPos[[1]]] <- rrbb$train[occ$species == 'Genista monspessulana'][trainPos[[1]]] <- 1
# occ$train[occ$species == 'Ulex europaeus'][trainPos[[2]]] <- rrbb$train[occ$species == 'Ulex europaeus'][trainPos[[2]]] <- 1

#### Clean records ------

prTrain <-  .7
maskk <- stkObs[[1]]
maskk[!is.na(maskk[])] <- 1:length(which(!is.na(maskk[])))
plot(maskk)

##

occ$ID <- rrbb$ID <- raster::extract(maskk, occ)
extValsID <- paste0(occ$species, occ$ID)
pos <- duplicated(extValsID)
cbind(extValsID, pos)
occ$dup <- rrbb$dup <- 0
occ$dup[pos] <- rrbb$dup[pos] <- 1
posNA <- grep('NA$', extValsID)
occ$dup[posNA] <- rrbb$dup[posNA] <- 1
table(occ$dup, occ$species)
table(rrbb$dup, rrbb$species)

seedNum <- 1000+500+747
set.seed(seedNum+u1)


## valid records
rrbbG <- rrbb[rrbb$species == 'Genista monspessulana' & rrbb$dup == 0, ]
occG  <-  occ[occ$species  == 'Genista monspessulana' & occ$dup == 0, ]

rrbbU <- rrbb[rrbb$species == 'Ulex europaeus' & rrbb$dup == 0, ]
occU  <-  occ[occ$species  == 'Ulex europaeus' & occ$dup == 0,  ]

## Generate pseudo absences
set.seed(nrow(rrbbG))
backgrG <- randomPoints(maskk, nrow(rrbbG))
set.seed(nrow(rrbbU))
backgrU <- randomPoints(maskk, nrow(rrbbU))

## Assign tags for train / test
trainPreU <- rep(0, nrow(occU))
testPreU <-  rep(0, nrow(occU))

nTrainRecU <- round(nrow(occU) * prTrain)
nTrainRecG <- round(nrow(occG) * prTrain)

## Sp 1: Ulex
set.seed(stringi::stri_reverse(nrow(rrbbU)))
pos_trainPresU <- sample(nrow(occU), nTrainRecU) # Pres train
set.seed(stringi::stri_reverse(nrow(rrbbU) + 1))
pos_trainAbsU <- sample(nrow(occU), nTrainRecU) # Abs train

## Sp 2
set.seed(stringi::stri_reverse(nrow(rrbbG)))
pos_trainPresG <- sample(nrow(occG), nTrainRecG) # Pres train
set.seed(stringi::stri_reverse(nrow(rrbbG) + 1))
pos_trainAbsG  <- sample(nrow(occG), nTrainRecG) # Abs train


pos_train_presG <- pos_train_absG <- rep(0, nrow(occG)) # Store in vector
pos_train_presU <- pos_train_absU <- rep(0, nrow(occU)) # Store in vector

pos_train_presG[pos_trainPresG] <- 1
pos_train_absG[pos_trainAbsG] <- 1

pos_train_presU[pos_trainPresU] <- 1
pos_train_absU[pos_trainAbsU] <- 1

tagsG <- c(pos_train_presG, pos_train_absG) # Store in vector
tagsU <- c(pos_train_presU, pos_train_absU) # Store in vector
rec$train <- pos_train_pres

oc <- rb <- rbind.data.frame(cbind.data.frame(rrbbG[, c('x', 'y', 'source', 'species', 'pa', 'one')],
                                        train = pos_train_presG),
                       cbind.data.frame(backgrG, source = 'pseudoAbs', species = 'Genista monspessulana',
                                        pa = 0, one = 1, train = pos_train_absG),

                       cbind.data.frame(rrbbU[, c('x', 'y', 'source', 'species', 'pa', 'one')],
                                        train = pos_train_presU),
                       cbind.data.frame(backgrU, source = 'pseudoAbs', species = 'Ulex europaeus',
                                        pa = 0, one = 1, train = pos_train_absU)
                       )

table( rb$train, rb$pa, rb$species)

coordinates(oc) =~ x + y
oc$x_ <- rb$x
oc$y_ <- rb$y

### Save
save(oc, rb, file = paste0(outDir, '03_Registros-especies/oc_rb_694_readyForModel.RData'))
save(oc, rb, file = paste0(outDir, '05_Modelos-distribucion/oc_rb_694_readyForModel.RData'))
writeOGR(oc, dsn = paste0(outDir, '03_Registros-especies'), layer = '694-presen-ausen_para_modelos',
         driver = 'ESRI Shapefile')

load( paste0(outDir, '03_Registros-especies/oc_rb_694_readyForModel.RData'))

## C:/Users/Admin/Documents/R/win-library/4.0/dismo/java

### Make stacks -----

load('H:/CI/Entregables_capas/stkObs.RData')
namesStk <- c('PpAnual',
              'TmAnual', 'TmDEF', 'TmMAM', 'TmJJA', 'TmSON',
              'TnAnual', 'TnDEF', 'TnMAM', 'TnJJA', 'TnSON',
              'TxAnual', 'TxDEF', 'TxMAM', 'TxJJA', 'TxSON')
names(stkObs) <- namesStk
stkObsS <- stkObs[[c('PpAnual', 'TmAnual', 'TnAnual', 'TxAnual')]]
stkObsTP <- stkObs[[c('PpAnual', 'TmAnual')]]
stkObsP <- stkObs[[c('PpAnual')]]
stkObsT <- stkObs[[c('TmAnual')]]
names(stkObsS)



###  -----
U_train <- oc[oc$species == 'Ulex europaeus' & oc$train == 1, 'pa']
U_test <- oc[oc$species == 'Ulex europaeus' & oc$train == 0, 'pa']

G_train <- oc[oc$species == 'Genista monspessulana' & oc$train == 1, 'pa']
G_test <- oc[oc$species == 'Genista monspessulana' & oc$train == 0, 'pa']


## Prepare data for models -----

#setwd('X:/Google Drive/CI/Entregables_Retamo/05_Modelos-distribucion/Models_R/')
setwd('X:/CI/')

system.time(UF.mod.dat <- BIOMOD_FormatingData(resp.var = U_train, eval.resp.var = U_test,
                                               expl.var = stkObs, eval.expl.var = stkObs,
                                               resp.name = "Ulex_full"))

system.time(US.mod.dat <- BIOMOD_FormatingData(resp.var = U_train, eval.resp.var = U_test,
                                               expl.var = stkObsS, eval.expl.var = stkObsS,
                                               resp.name = "Ulex_small"))

system.time(UTP.mod.dat <- BIOMOD_FormatingData(resp.var = U_train, eval.resp.var = U_test,
                                                expl.var = stkObsTP, eval.expl.var = stkObsTP,
                                                resp.name = "Ulex_PT"))

system.time(UP.mod.dat <- BIOMOD_FormatingData(resp.var = U_train, eval.resp.var = U_test,
                                                expl.var = stkObsP, eval.expl.var = stkObsP,
                                                resp.name = "Ulex_P"))

system.time(UT.mod.dat <- BIOMOD_FormatingData(resp.var = U_train, eval.resp.var = U_test,
                                                expl.var = stkObsT, eval.expl.var = stkObsT,
                                                resp.name = "Ulex_T"))

 ## Genista
system.time(GF.mod.dat <- BIOMOD_FormatingData(resp.var = G_train, eval.resp.var = G_test,
                                               expl.var = stkObs, eval.expl.var = stkObs,
                                               resp.name = "Genista_full"))

system.time(GS.mod.dat <- BIOMOD_FormatingData(resp.var = G_train, eval.resp.var = G_test,
                                               expl.var = stkObsS, eval.expl.var = stkObsS,
                                               resp.name = "Genista_small"))

system.time(GTP.mod.dat <- BIOMOD_FormatingData(resp.var = G_train, eval.resp.var = G_test,
                                                expl.var = stkObsTP, eval.expl.var = stkObsTP,
                                                resp.name = "Genista_PT"))

system.time(GP.mod.dat <- BIOMOD_FormatingData(resp.var = G_train, eval.resp.var = G_test,
                                               expl.var = stkObsP, eval.expl.var = stkObsP,
                                               resp.name = "Genista_P"))

system.time(GT.mod.dat <- BIOMOD_FormatingData(resp.var = G_train, eval.resp.var = G_test,
                                               expl.var = stkObsT, eval.expl.var = stkObsT,
                                               resp.name = "Genista_T"))



save(UF.mod.dat, US.mod.dat, UTP.mod.dat, UP.mod.dat, UT.mod.dat,
     GF.mod.dat, GS.mod.dat, GTP.mod.dat, GP.mod.dat, GT.mod.dat,
     file = 'mod.Data.Models_2spp-6Meth-4VarSets.RData')
load('mod.Data.Models_2spp-6Meth-4VarSets.RData')



###### Run models ----------

models2run <- c('GLM','ANN', 'GBM', 'SRE', 'RF', 'MAXENT.Phillips') # c('MAXENT.Phillips', 'MAXENT.Phillips.2')

system.time(UF.mod <- BIOMOD_Modeling(data = UF.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(US.mod <- BIOMOD_Modeling(data = US.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(UTP.mod <- BIOMOD_Modeling(data = UTP.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(UP.mod <- BIOMOD_Modeling(data = UP.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(UT.mod <- BIOMOD_Modeling(data = UT.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))


system.time(GF.mod <- BIOMOD_Modeling(data = GF.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(GS.mod <- BIOMOD_Modeling(data = GS.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(GTP.mod <- BIOMOD_Modeling(data = GTP.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(GP.mod <- BIOMOD_Modeling(data = GP.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))
system.time(GT.mod <- BIOMOD_Modeling(data = GT.mod.dat, models = models2run, SaveObj = TRUE, VarImport = 5))


# save(UF.mod, US.mod, UTP.mod, UP.mod, UT.mod,
#      GF.mod, GS.mod, GTP.mod, GP.mod, GT.mod,
#      file = 'mod.Models_2spp-6Meth-4VarSets.RData')
load('mod.Models_2spp-6Meth-4VarSets.RData')

#### Gather all metrics ------------
UF.mod.eval <- get_evaluations(UF.mod)
US.mod.eval <- get_evaluations(US.mod)
UTP.mod.eval <- get_evaluations(UTP.mod)
UP.mod.eval <- get_evaluations(UP.mod)
UT.mod.eval <- get_evaluations(UT.mod)

GF.mod.eval <- get_evaluations(GF.mod)
GS.mod.eval <- get_evaluations(GS.mod)
GTP.mod.eval <- get_evaluations(GTP.mod)
GP.mod.eval <- get_evaluations(GP.mod)
GT.mod.eval <- get_evaluations(GT.mod)

# UF.mod.eval["TSS","Testing.data",,,]
# UF.mod.eval["ROC","Testing.data",,,]
# UF.mod.eval["KAPPA","Testing.data",,,]

U.mod.eval.melt <- rbind(cbind.data.frame(reshape2::melt(UF.mod.eval), vars = 'Full'),
                         cbind.data.frame(reshape2::melt(US.mod.eval), vars = 'PpTmTxTn'),
                         cbind.data.frame(reshape2::melt(UTP.mod.eval), vars = 'Pp-Tm'),
                         cbind.data.frame(reshape2::melt(UP.mod.eval), vars = 'Pp'),
                         cbind.data.frame(reshape2::melt(UT.mod.eval), vars = 'Tm'))

G.mod.eval.melt <- rbind(cbind.data.frame(reshape2::melt(GF.mod.eval), vars = 'Full'),
                         cbind.data.frame(reshape2::melt(GS.mod.eval), vars = 'PpTmTxTn'),
                         cbind.data.frame(reshape2::melt(GTP.mod.eval), vars = 'Pp-Tm'),
                         cbind.data.frame(reshape2::melt(GP.mod.eval), vars = 'Pp'),
                         cbind.data.frame(reshape2::melt(GT.mod.eval), vars = 'Tm'))

U.mod.eval.melt$Var3 <- gsub("MAXENT.Phillips.2|MAXENT.Phillips", 'MAX', U.mod.eval.melt$Var3)
G.mod.eval.melt$Var3 <- gsub("MAXENT.Phillips.2|MAXENT.Phillips", 'MAX', G.mod.eval.melt$Var3)

colnames(U.mod.eval.melt) <-
  colnames(G.mod.eval.melt) <- c('Metric', 'Object', 'Model', 'x', 'y', 'Value', 'Vars')

head(U.mod.eval.melt)

U.mod.eval.melt$Vars <- factor(U.mod.eval.melt$Vars,
                               levels = c("Full", "PpTmTxTn",  "Pp-Tm", "Pp", "Tm" ))

G.mod.eval.melt$Vars <- factor(G.mod.eval.melt$Vars,
                               levels = c("Full", "PpTmTxTn",  "Pp-Tm", "Pp", "Tm" ))


mod.var.imp0 <- rbind.data.frame(
  ## Genista
  cbind.data.frame(reshape2::melt(drop(get_variables_importance(GF.mod))), vars = 'Full', sp = 'Genista monspessulana'),
                     cbind.data.frame(reshape2::melt(drop(get_variables_importance(GS.mod))), vars = 'PpTmTxTn', sp = 'Genista monspessulana'),
                     cbind.data.frame(reshape2::melt(drop(get_variables_importance(GTP.mod))), vars = 'Pp-Tm', sp = 'Genista monspessulana'),
                     cbind.data.frame(reshape2::melt((get_variables_importance(GP.mod)))[, c(1,2,5)], vars = 'Pp', sp = 'Genista monspessulana'),
                     cbind.data.frame((reshape2::melt((get_variables_importance(GT.mod))))[, c(1,2,5)], vars = 'Tm', sp = 'Genista monspessulana'),
    ## Ulex
                     cbind.data.frame(reshape2::melt(drop(get_variables_importance(UF.mod))), vars = 'Full', sp = 'Ulex europaeus'),
                     cbind.data.frame(reshape2::melt(drop(get_variables_importance(US.mod))), vars = 'PpTmTxTn', sp = 'Ulex europaeus'),
                     cbind.data.frame(reshape2::melt(drop(get_variables_importance(UTP.mod))), vars = 'Pp-Tm', sp = 'Ulex europaeus'),
                     cbind.data.frame((reshape2::melt((get_variables_importance(UP.mod))))[, c(1,2,5)], vars = 'Pp', sp = 'Ulex europaeus')
)

UT.mod.VI <- get_variables_importance(UT.mod); rownames(UT.mod.VI) <- 'TmAnual'

mod.var.imp <- rbind(mod.var.imp0,
                     cbind.data.frame(reshape2::melt((UT.mod.VI))[, c(1,2,5)], vars = 'Tm', sp = 'Ulex europaeus'))

# head(mod.var.imp)
# head(cbind.data.frame(reshape2::melt((get_variables_importance(UT.mod)))[, c(1,2,5)], vars = 'Tm', sp = 'Ulex europaeus'))

colnames(mod.var.imp) <- c('Var', 'Model', 'Value', 'Set', 'Sp')
mod.var.imp$Set <- factor(mod.var.imp$Set,
                               levels = c("Full", "PpTmTxTn",  "Pp-Tm", "Pp", "Tm" ))


table(mod.var.imp$Var, mod.var.imp$Set, mod.var.imp$Sp)

## Metrics
save(mod.var.imp, U.mod.eval.melt, G.mod.eval.melt,
     file = 'summary.Models_2spp-6Meth-4VarSets.RData')

# save.image('E:/CI/models_Tests_biomod.RData')



#### Projections -----

#setwd('X:/Google Drive/CI/Entregables_Retamo/05_Modelos-distribucion/Models_R/')
#setwd('X:/CI/')

#  2403.93 / 40

system.time(UF.pred <- BIOMOD_Projection(modeling.output = UF.mod,
                                         new.env = stkObs, # modern environment
                                         proj.name = 'current' ,
                                         selected.models = 'all' , binary.meth = 'TSS' ,
                                         compress = 'xz' , clamping.mask = F, output.format = '.RData' ))

save(UF.pred, file = 'UF.pred_obs.RData')
origNames <- UF.pred@models.projected
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', UF.pred@models.projected)))
tifBase <- paste0('Ulex_obs_Full_')
writeRaster(UF.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(UF.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(UF.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(UF.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(UF.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(UF.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(UF.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

UF.pred@models.projected <- namesTif
UF.pred@proj.names <- 'Obs'
UF.pred@sp.name <- 'Ulex - Full -'
plot(UF.pred)


system.time(GF.pred <- BIOMOD_Projection(modeling.output = GF.mod, # 2517.17 /60 = [1] 41.95283
                                         new.env = stkObs, # modern environment
                                         proj.name = 'current' ,
                                         selected.models = 'all' , binary.meth = 'TSS' ,
                                         build.clamping.mask = F,
                                         output.format = '.RData' ))

#Error in res[i] <- readBin(x@file@con, what = dtype, n = 1, size = dsize,  :

save(GF.pred, file = 'GF.pred_obs.RData')
origNames <- GF.pred@models.projected
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', GF.pred@models.projected)))
tifBase <- paste0('Geni_obs_Full_')
writeRaster(GF.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(GF.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(GF.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(GF.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(GF.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(GF.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(GF.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

GF.pred@models.projected <- namesTif
GF.pred@proj.names <- 'Obs'
GF.pred@sp.name <- 'Genista - Full -'
#plot(GF.pred)


## Small ---------------------------------------

system.time(US.pred <- BIOMOD_Projection(modeling.output = US.mod, #  1177.14 /60 = [1] 19.619
                                         new.env = stkObsS,  #
                                         proj.name = 'current' , #
                                         selected.models = 'all' , #
                                         binary.meth = 'TSS' ,
                                         build.clamping.mask = T,
                                         output.format = '.RData' ))

save(US.pred, file = 'US.pred_obs.RData') #
origNames <- US.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', US.pred@models.projected))) #
tifBase <- paste0('Ulex_obs_Small_') #
writeRaster(US.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(US.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(US.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(US.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(US.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(US.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(US.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

US.pred@models.projected <- namesTif
US.pred@proj.names <- 'Obs'
US.pred@sp.name <- 'Ulex - Small -'
#plot(US.pred)

## --------------------------------------
#tempdir()
system.time(GS.pred <- BIOMOD_Projection(modeling.output = GS.mod, # 1224.11/60 [1] 20.40183
                                         new.env = stkObsS,  #
                                         proj.name = 'current' , #
                                         selected.models = 'all' , #
                                         binary.meth = 'TSS' ,
                                         build.clamping.mask = F,
                                         output.format = '.RData' ))

save(GS.pred, file = 'GS.pred_obs.RData') #
origNames <- GS.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', GS.pred@models.projected))) #
tifBase <- paste0('Geni_obs_Small_') #
writeRaster(GS.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(GS.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(GS.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(GS.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(GS.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(GS.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(GS.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

GS.pred@models.projected <- namesTif
GS.pred@proj.names <- 'Obs'
GS.pred@sp.name <- 'Genista - Small -'
#plot(GS.pred)




## PT ---------------------------------------

system.time(UTP.pred <- BIOMOD_Projection(modeling.output = UTP.mod, #
                                          new.env = stkObsTP,  #
                                          proj.name = 'current' , #
                                          selected.models = 'all' , #
                                          binary.meth = 'TSS' ,
                                          build.clamping.mask = F,
                                          output.format = '.RData' ))

save(UTP.pred, file = 'UTP.pred_obs.RData') #
origNames <- UTP.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', UTP.pred@models.projected))) #
tifBase <- paste0('Ulex_obs_TP_') #
writeRaster(UTP.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(UTP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

UTP.pred@models.projected <- namesTif
UTP.pred@proj.names <- 'Obs'
UTP.pred@sp.name <- 'Ulex - TP -'
plot(UTP.pred)

## PT --------------------------------------

system.time(GTP.pred <- BIOMOD_Projection(modeling.output = GTP.mod, #  1107.56 /60 [1] 18.45933
                                          new.env = stkObsTP,  #
                                          proj.name = 'current' , #
                                          selected.models = 'all' , #
                                          binary.meth = 'TSS' ,
                                          build.clamping.mask = F,
                                          output.format = '.RData' ))

save(GTP.pred, file = 'GTP.pred_obs.RData') #
origNames <- GTP.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', GTP.pred@models.projected))) #
tifBase <- paste0('Geni_obs_TP_') #
writeRaster(GTP.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(GTP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

GTP.pred@models.projected <- namesTif
GTP.pred@proj.names <- 'Obs'
GTP.pred@sp.name <- 'Genista - Small -'
plot(GS.pred)


## P ---------------------------------------

system.time(UP.pred <- BIOMOD_Projection(modeling.output = UP.mod, #
                                          new.env = stack(stkObsP),  #
                                          proj.name = 'current' , #
                                          selected.models = 'all' , #
                                          binary.meth = 'TSS' ,
                                          build.clamping.mask = F,
                                          output.format = '.RData' ))

save(UP.pred, file = 'UP.pred_obs.RData') #
origNames <- UP.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', UP.pred@models.projected))) #
tifBase <- paste0('Ulex_obs_P_') #
writeRaster(UP.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(UP.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(UP.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(UP.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(UP.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(UP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(UP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

UP.pred@models.projected <- namesTif
UP.pred@proj.names <- 'Obs'
UP.pred@sp.name <- 'Ulex - TP -'
#plot(UP.pred)

## P --------------------------------------

system.time(GP.pred <- BIOMOD_Projection(modeling.output = GP.mod, #  1107.56 /60 [1] 18.45933
                                          new.env = stack(stkObsP),  #
                                          proj.name = 'current' , #
                                          selected.models = 'all' , #
                                          binary.meth = 'TSS' ,
                                          build.clamping.mask = F,
                                          output.format = '.RData' ))

save(GP.pred, file = 'GP.pred_obs.RData') #
origNames <- GP.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', GP.pred@models.projected))) #
tifBase <- paste0('Geni_obs_P_') #
writeRaster(GP.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(GP.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(GP.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(GP.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(GP.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(GP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(GP.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

GP.pred@models.projected <- namesTif
GP.pred@proj.names <- 'Obs'
GP.pred@sp.name <- 'Genista - Small -'
#plot(GS.pred)

## T ---------------------------------------

system.time(UT.pred <- BIOMOD_Projection(modeling.output = UT.mod, #
                                         new.env = stack(stkObsT),  #
                                         proj.name = 'current' , #
                                         selected.models = 'all' , #
                                         binary.meth = 'TSS' ,
                                         build.clamping.mask = F,
                                         output.format = '.RData' ))

save(UT.pred, file = 'UT.pred_obs.RData') #
origNames <- UT.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', UT.pred@models.projected))) #
tifBase <- paste0('Ulex_obs_T_') #
writeRaster(UT.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(UT.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(UT.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(UT.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(UT.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(UT.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(UT.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

UT.pred@models.projected <- namesTif
UT.pred@proj.names <- 'Obs'
UT.pred@sp.name <- 'Ulex - TP -'
#plot(UP.pred)

## --------------------------------------

system.time(GT.pred <- BIOMOD_Projection(modeling.output = GT.mod, #  1107.56 /60 [1] 18.45933
                                         new.env = stack(stkObsT),  #
                                         proj.name = 'current' , #
                                         selected.models = 'all' , #
                                         binary.meth = 'TSS' ,
                                         build.clamping.mask = F,
                                         output.format = '.RData' ))

save(GT.pred, file = 'GT.pred_obs.RData') #
origNames <- GT.pred@models.projected #
(namesTif <- gsub('MAXENT.Phillips|MAXENT.Phillips2', 'MAX', gsub('.+Full_', '', GT.pred@models.projected))) #
tifBase <- paste0('Geni_obs_T_') #
writeRaster(GT.pred@proj@val@layers[[1]], filename = paste0(tifBase, namesTif[1], '.tif'))
writeRaster(GT.pred@proj@val@layers[[2]], filename = paste0(tifBase, namesTif[2], '.tif'))
writeRaster(GT.pred@proj@val@layers[[3]], filename = paste0(tifBase, namesTif[3], '.tif'))
writeRaster(GT.pred@proj@val@layers[[4]], filename = paste0(tifBase, namesTif[4], '.tif'))
writeRaster(GT.pred@proj@val@layers[[5]], filename = paste0(tifBase, namesTif[5], '.tif'))
writeRaster(GT.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))
writeRaster(GT.pred@proj@val@layers[[6]], filename = paste0(tifBase, namesTif[6], '.tif'))

GT.pred@models.projected <- namesTif
GT.pred@proj.names <- 'Obs'
GT.pred@sp.name <- 'Genista - Small -'
#plot(GS.pred)

