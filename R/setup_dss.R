
instPack <- function(x, devt = FALSE){
  sapply(x, FUN = function(y){

    lib2inst <- ifelse(devt, basename(y), y)
    rnIl <- rownames(installed.packages());
    # if( require( eval(parse(text=paste(y))) ) ){}

    if(!lib2inst %in% rnIl) {
      if (!devt){
        install.packages(lib2inst)
      } else {
        devtools::install_github(y)
      }
    }
  })
}

setup_cola_dss <- function(

   libs2colaDSS = c(

  'markdown', 'rmarkdown',
  'knitr', 'units', 'reticulate',

  "reshape2", 'bit', 'digest', 'dplyr',
  'tidyverse', 'DT', 'ggplot2', 'data.table',

  # debug install order: htmltools >> shiny >> shinyWidgets
  'htmlwidgets', 'htmltools', ## Before shiny
  'magrittr', 'RColorBrewer', 'viridis',

  ## Spatial
  # "rgeos", "rgdal", 'raster',
  'sf', 'terra',
  'rlang', "leaflet", "leaflet.extras",
  'gdalUtilities',

  ## Shiny
  'shiny',  'shinyBS',## Before shiny plugins
  "shinydashboard",  "shinycssloaders",
  'shinydashboardPlus', 'shinyjs',
  'shinyWidgets', 'dashboardthemes',
  "highcharter", 'plotly')){

  if(Sys.info()["sysname"] == 'Linux'){
    cat('   Consider install the next libraries in Linux console before installing R packages: \n    ',
    'sudo apt -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libudunits2-dev')
    Sys.sleep(5)
  }

  if (!require("devtools")) {
    install.packages("devtools")
  }
  library(devtools)

  if(!find_rtools()){
    warning(" No Rtools found. Please install it for your R version")
  }

  if (!require("remotes")) {
    install.packages("remotes")
  }

  cola::instPack(libs2colaDSS )

  # if(!require('gdalUtils') ){
  #   ## Install from Github
  #   devtools:::install_github("gearslaboratory/gdalUtils")
  # }

  if (!require("shinydashboard")) {
    instPack("rstudio/shinydashboard", devt = TRUE)
  }

  if (!require("shinySignals")) {
    instPack("hadley/shinySignals")
    if(!require('shinySignals') ){
      remotes::install_github("hadley/shinySignals")
    }
  }

  # remotes::install_version("rjson", "0.2.20")
  if (!require("shinySignals")) {
    instPack("highcharter") ## CRAN
  }

  if(!require('highcharter') ){
    (remotes::install_github("jbkunst/highcharter"))
  }

  library(bit) #
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(gdalUtilities)
  library(highcharter)
  library(htmlwidgets)
  library(htmltools)
  library(leaflet)
  library(leaflet.extras)
  library(knitr)
  library(magrittr)
  #library(maptools) # deprecated
  #library(raster) # deprecated
  library(RColorBrewer)
  #library(rgdal) # deprecated
  #library(rgeos) # deprecated
  library(rmarkdown)
  library(sf)
  library(shiny)
  library(shinyBS)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(shinyjs)
  library(shinyWidgets)
  library(dashboardthemes)
  library(shinycssloaders)
  library(tidyverse)
  library(shiny)
  library(reshape2)
  library(DT)
  library(tibble)
  library(terra)
  library(viridis)

  cat("\n   === All libraries required for COLA's DSS installed === \n\n")
}
