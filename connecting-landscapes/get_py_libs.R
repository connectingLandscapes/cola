
## Your git path
setwd('N:/Mi unidad/git/connecting-landscapes/src/')
pyFiles <- list.files(pattern = '.py$')

impLines <- NULL
for (i in 1:length(pyFiles)){ # i = 1
  py.i <- readLines(pyFiles[i], ok = TRUE)
  line.i <- grep('import', py.i, value = TRUE)
  impLines <- rbind.data.frame(impLines, 
                    cbind(file = pyFiles[i], 
                          orig = line.i))
}

uCmd <- sort(unique(impLines$orig))
uCmd <- grep('streamlit|ipysheet',  uCmd, value = TRUE, invert = TRUE)

(uCmd <- c('import osgeo', uCmd, 'print ("WELCOME -- libraries loaded succesfully")'))
write.table(x = uCmd, file = 'welcome.py', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = uCmd, quote = FALSE, row.names = FALSE, col.names = FALSE,
            file = 'N:/Mi unidad/git/connectscap src/welcome.py')

# system(paste0('C:\\Users\\Admin\\AppData\\Local\\r-miniconda\\envs\\cola\\python.exe ',
#       'N:/"Mi unidad"\\git\\connecting-landscapes\\src\\welcome.py'))


library(reticulate)

(condaLists <- reticulate::conda_list())
(pyCola <- subset(condaLists, name == 'cola')$python)

folñe
system(paste(pyCola, ' welcome.py'))
