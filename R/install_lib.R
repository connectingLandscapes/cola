
instPack <- function(x, devt = FALSE){
  sapply(x, FUN = function(y){

    lib2inst <- ifelse(devt, basename(y), y)
    rnIl <- rownames(installed.packages());
    # if( require( eval(parse(text=paste(y))) ) ){}

    if(!lib2inst %in% rnIl) {
      if (!devt){
        tryCatch(install.packages(lib2inst), function(e) e)
      } else {
        tryCatch( devtools::install_github(y), function(e) e)
      }
    }
  })
}

setup_cola_dss <- function(){

  if(Sys.info()["sysname"] == 'Linux'){
    cat('   Consider install the next libraries in Linux console before installing R packages: \n    ',
    'sudo apt -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libudunits2-dev')
    Sys.sleep(10)
  }

  if (!require("devtools")) {
    install.packages("devtools")
  }
  library(devtools)


  if (!require("remotes")) {
    install.packages("remotes")
  }


  instPack(c(
    # "rgeos", "rgdal", 'raster',
    'markdown', 'rmarkdown',
    'knitr', 'units',
    "reshape2", 'bit', 'digest', 'dplyr',
    'tidyverse', 'DT', 'ggplot2',

    # debug install order: htmltools >> shiny >> shinyWidgets

    'htmlwidgets', 'htmltools', ## Before shiny
    'magrittr', 'RColorBrewer',
    'sf', 'terra',
    'rlang', "leaflet", "leaflet.extras",
    'shiny',  ## Before shiny plugins
    "shinydashboard",  "shinycssloaders",
    'shinydashboardPlus', 'shinyjs',
    'shinyWidgets', 'dashboardthemes',
    "highcharter", 'plotly'))

  if( require('gdalUtils') ){
    ## Install from Github
    devtools:::install_github("gearslaboratory/gdalUtils")
  }


  instPack("rstudio/shinydashboard", devt = TRUE)
  instPack("hadley/shinySignals")
  # remotes::install_version("rjson", "0.2.20")
  instPack("highcharter") ## CRAN

  if( require('gdalUtils') ){
    tryCatch(remotes::install_github("jbkunst/highcharter"))
  }

  library(bit) #
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
  #library(maptools) # deprecated
  #library(raster) # deprecated
  library(RColorBrewer)
  #library(rgdal) # deprecated
  #library(rgeos) # deprecated
  library(rmarkdown)
  library(sf)
  library(shiny)
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

  cat("      All libraries required for COLA's DSS installed")
}
