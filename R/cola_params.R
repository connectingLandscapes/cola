#' Parameters for installation
#'
#' Default parameters for cola conda environment installation
#' @return The temperature in degrees Celsius
#' @export
cola_params <<- list(
  ## Environment name
  envName = 'cola',
  ## Libraries to install
  libs2Install = c('gdal', 'h5py', # 'osgeo',
                   'numexpr',
                   'rasterio', 'pytables',
                   'pandas',  'cython', 'numba' ,
                   'networkit', 'fiona', 'shapely',
                   'geopandas',
                   'kdepy', # 'KDEpy',
                   'scikit-image'),
  yml = TRUE,
  ## Number steps
  nSteps = 5,
  cola_params <<- list(
    ## Environment name
    envName = 'cola',

    ## Libraries to install in the python conda environment
    libs2Install = c('gdal', 'h5py', # 'osgeo',
                     'numexpr',
                     'rasterio', 'pytables',
                     'pandas',  'cython', 'numba' ,
                     'networkit', 'fiona', 'shapely',
                     'geopandas',
                     'kdepy', # 'KDEpy',
                     'scikit-image'),
    yml = TRUE,
    ## Number steps
    nSteps = 5,

    ## R packages for the DSS
    libs2colaDSS = c(
      'markdown', 'rmarkdown',
      'knitr', 'units',
      "reshape2", 'bit', 'digest', 'dplyr',
      'tidyverse', 'DT', 'ggplot2',
      # debug install order: htmltools >> shiny >> shinyWidgets
      'htmlwidgets', 'htmltools', ## Before shiny
      'magrittr', 'RColorBrewer', 'viridis',
      ## Spatial
      # "rgeos", "rgdal", 'raster', ## old
      'sf', 'terra',
      'rlang', "leaflet", "leaflet.extras",
      'gdalUtilities',
      ## Shiny
      'shiny',  ## Before shiny plugins
      "shinydashboard",  "shinycssloaders",
      'shinydashboardPlus', 'shinyjs',
      'shinyWidgets', 'dashboardthemes',
      "highcharter", 'plotly')
  )
)

# setwd('N:/Mi unidad/git/cola/data/')
# setwd( 'N:/My Drive/git/cola/data/')
# save(cola_params, file = 'cola_params.RData')
