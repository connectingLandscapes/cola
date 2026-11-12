#' @title Install R packages in batch
#' @description Installs R packages, trying CRAN, then github.
#' @param x String. Name of the libraries to install. 
#' @param devt Logical. Default FALSE. Is from a developer (github)?  
#' @return NULL. Prints in console logs regarding different steps
#' @examples
#' if (FALSE){ 
#'   setup_cola_dss( )
#' }
#' @author Ivan Gonzalez <Ivan.Gonzalez@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

instPack <- function(x, devt = FALSE){
  sapply(x, FUN = function(y){

    lib2inst <- ifelse(devt, basename(y), y)
    rnIl <- rownames(utils::installed.packages());
    # if( require( eval(parse(text=paste(y))) ) ){}

    if(!lib2inst %in% rnIl) {
      if (!devt){
        tryCatch(utils::install.packages(lib2inst), error = function(e) print(e))
      } else {
        devtools::install_github(y)
      }
    }
  })
}


#' @title Install  \emph{COLA} graphical user interface
#' @description Installs cola GUI. Requires several libraries to display maps and widgets.
#' @param libs2colaDSS String. Name of the libraries to install. 
#' @return NULL. Prints in console logs regarding different steps
#' @examples
#' if (FALSE){ 
#'   setup_cola_dss( )
#' }
#' @author Ivan Gonzalez <Ivan.Gonzalez@@nau.edu>
#' @author Patrick Jantz <Patrick.Jantz@@gmail.com>
#' @export

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
  'rlang', "leaflet",
  # "leaflet.extras",
  'gdalUtilities',

  ## Shiny
  'shiny',  'shinyBS',## Before shiny plugins
  "shinydashboard",  "shinycssloaders",
  'shinydashboardPlus', 'shinyjs',
  'shinyalert', 'shinyFiles',
  'shinyWidgets', 'dashboardthemes',
  "highcharter", 'plotly')){

  if(Sys.info()["sysname"] == 'Linux'){
    cat('\n   Only for Linux: Install the next libraries in Linux console before installing R DSS packages (not required for command line): \n    ',
    'sudo apt -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libudunits2-dev\n')
    Sys.sleep(5)
  }

  if (!requireNamespace("devtools")) {
    utils::install.packages("devtools")
  }
  requireNamespace("devtools")

  if(Sys.info()["sysname"] == 'Windows'){
    if(!devtools::find_rtools()){
      warning(" No Rtools found. Please install it for your R version")
    }
  }

  if (!requireNamespace("remotes")) {
    utils::install.packages("remotes")
  }

  cola::instPack(libs2colaDSS )

  # if(!require('gdalUtils') ){
  #   ## Install from Github
  #   devtools:::install_github("gearslaboratory/gdalUtils")
  # }

  if (!requireNamespace("shinydashboard")) {
    instPack("rstudio/shinydashboard", devt = TRUE)
  }

  if (!requireNamespace("shinySignals")) {
    instPack("hadley/shinySignals")
    if(!requireNamespace('shinySignals') ){
      remotes::install_github("hadley/shinySignals")
    }
  }

  # remotes::install_version("rjson", "0.2.20")
  if (!requireNamespace("highcharter")) {
    instPack("highcharter") ## CRAN
  } 

  if(!requireNamespace('highcharter') ){
    (remotes::install_github("jbkunst/highcharter"))
  }

  if(!requireNamespace('leaflet.extras') ){
    remotes::install_github("trafficonese/leaflet.extras")
  }

  requireNamespace("bit") #
  requireNamespace("digest")
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("gdalUtilities")
  requireNamespace("highcharter")
  requireNamespace("htmlwidgets")
  requireNamespace("htmltools")
  requireNamespace("leaflet")
  requireNamespace("leaflet.extras")
  requireNamespace("knitr")
  requireNamespace("magrittr")
  #requireNamespace("maptools") # deprecated
  #requireNamespace("raster") # deprecated
  requireNamespace("RColorBrewer")
  #requireNamespace("rgdal") # deprecated
  #requireNamespace("rgeos") # deprecated
  requireNamespace("rmarkdown")
  requireNamespace("sf")
  requireNamespace("shiny")
  requireNamespace("shinyalert")
  requireNamespace("shinyBS")
  requireNamespace("shinydashboard")
  requireNamespace("shinydashboardPlus")
  requireNamespace("shinyjs")
  requireNamespace("shinyWidgets")
  requireNamespace("dashboardthemes")
  requireNamespace("shinycssloaders")
  requireNamespace("tidyverse")
  requireNamespace("shiny")
  requireNamespace("reshape2")
  requireNamespace("DT")
  requireNamespace("tibble")
  requireNamespace("terra")
  requireNamespace("viridis")

  #cat("\n\n   === All libraries required for COLA's DSS installed === \n\n")
  message("   === All libraries required for COLA's DSS installed === ")
}
