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
  nSteps = 5
)
