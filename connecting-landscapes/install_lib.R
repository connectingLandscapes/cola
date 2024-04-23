### ------------- v2

system('sudo apt -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libudunits2-dev')

install.packages("devtools")

instPack <- function(x){
  sapply(x, FUN = function(y){
    rn=rownames(installed.packages());
    if(!y %in% rownames(installed.packages())) {install.packages(y)}
  })
}

instPack(c(
  # "rgeos", "rgdal", 'raster', 
  'markdown', 'rmarkdown', 'knitr',
  "devtools", 'units',
  "reshape2", 'bit', 'digest', 'dplyr',
  'tidyverse', 'DT', 'ggplot2', 
  'htmlwidgets', 'htmltools', 'magrittr', 'RColorBrewer',  
  
  'sf', 'terra',
  'rlang', "leaflet", "leaflet.extras", 
  'shiny', "shinydashboard",  "shinycssloaders", 
  'shinydashboardPlus', 'shinyjs', 'shinyWidgets', 'dashboardthemes',
  'shiny', 'shinydashboard', 'shinycssloaders', 
  "highcharter", 'plotly' ))


## Install from Github
devtools:::install_github("gearslaboratory/gdalUtils")
# rgeos,. matools, reshape2

#remotes::install_github("r-earthengine/rgeeExtra")
#####

### check
library(bit)
library(digest)
library(dplyr)
library(ggplot2)
library(highcharter)
library(htmlwidgets)
library(htmltools)
library(leaflet)
library(leaflet.extras)
library(knitr)
library(magrittr)
library(raster)
library(RColorBrewer)
library(rgdal)
library(rmarkdown)
library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(plotly)



### ------------- v0 Ignore
if (FALSE){
  
  # sudo apt -y install libfontconfig1-dev
  
  sapply(c('crayon', 'systemfonts'), function(x){ {install.packages(x)} })
  sapply(c('R6', 'httr'), function(x){ {install.packages(x)} })
  sapply(c('rlang', "openssl", "curl", "lifecycle", "pillar"), 
         function(x){ {install.packages(x)} })
  install.packages(c('magrittr', 'cli'))
  install.packages(c('jsonlite'))
  install.packages(c('ps', 'fastmap', 'glue', 'processx', 'cachem'))
  install.packages(c('fs', 'withr'))
  install.packages(c('stringi'))
  install.packages(c('digest', 'stringr'))
  install.packages(c('evaluate','highr', 'xfun', 'knitr'))
  install.packages(c('Rcpp', 'roxygen2'))
  
  install.packages('devtools')
  
  rpi <- rownames(installed.packages())
  sapply(c('rlang', "leaflet", "leaflet.extras", "rgeos", "rgdal", "shinydashboard", "highcharter", 
           "devtools", "shinycssloaders", "reshape2", 'bit', 'digest', 'dplyr', 'ggplot2', 'highcharter', 
           'htmlwidgets', 'htmltools', 'leaflet', 'leaflet.extras', 'knitr', 'magrittr', 'raster', 'RColorBrewer', 'rgdal', 
           'rmarkdown', 'shiny', 'shinydashboard', 'shinycssloaders', 
           'plotly'), function(x){ if (rpi  %in% x) {install.packages(x)} else {0}})
  
  install.packages(c('htmltools', 'promises'))
  install.packages(c( 'httpuv'))
  install.packages(c( 'xtable'))
  install.packages(c( 'shiny', 'shinydashboard'))
  
  
  # https://installati.one/install-r-cran-devtools-ubuntu-22-04/ # sudo apt-get update # sudo apt-get -y install r-cran-devtools
  
  library(devtools)
  #devtools::install_github("rstudio/shinydashboard")
  devtools::install_github("jcheng5/bubbles")
  devtools::install_github("hadley/shinySignals")
  #devtools::install_github("ramnathv/rCharts")
  
  if (!require("remotes")) {
    install.packages("remotes")
  }
  remotes::install_github("jcrodriguez1989/shinyParallel")
  #devtools::install_github("dkahle/ggmap")
  
  ## https://stackoverflow.com/questions/71748279/error-while-downloading-highcharter-package-in-r
  remotes::install_version("rjson", "0.2.20")
  install.packages("highcharter") ## CRAN
  remotes::install_github("jbkunst/highcharter")
  
  # debug insall order: htmltools >> shiny >> shinyWidgets
  
  # install.packages('shinyWidgets')
  install.packages('shinydashboardPlus')
  # install.packages('dashboardthemes')
  
  install.packages(c( 'bit', 'dplyr'))
  install.packages(c( 'ggplot2'))
  install.packages("yaml")
  install.packages("highcharter")
  
  install.packages('rgdal')
  install.packages('raster')
  install.packages('markdown')
  install.packages('rmarkdown')
  install.packages(c( 'leaflet', 'leaflet.extras'))
  install.packages('plotly')
  install.packages('shinycssloaders')
  install.packages(c('jquerylib', 'bslib'))
  install.packages('rgl')
  install.packages('maptools')
  #sudo apt-get install -y libudunits2-dev
  install.packages('units')
  install.packages('sf')
  install.packages('DT')
  #sudo apt-get install -y libharfbuzz-dev libfribidi-dev
  
  install.packages('devtools')
  install.packages('shinyjs')
  install.packages('shinyWidgets')
  install.packages('dashboardthemes')
  install.packages('tidyverse')
  install.packages('reshape2')
  install.packages('markdown')
}