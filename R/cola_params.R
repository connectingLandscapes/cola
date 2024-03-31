# Parameters for installation

cola_params <- list(
  ## Environment name
  envName = 'cola',
  ## Libraries to install
  libs2Install = c('gdal', 'h5py', # 'osgeo',
                   'rasterio', 'pytables',
                   'pandas',  'cython', 'numba' ,
                   'networkit', 'fiona', 'shapely',
                   'geopandas',
                   'kdepy', # 'KDEpy',
                   'scikit-image'),
  ## Number steps
  nSteps <- 5
)
