
## Definir lista de librerias a instalar
packages <- c('raster', 'rgdal', 'terra', 'sf', 'gdalUtilities', 'gdalUtils', 'sp', 'foreign')

## Revisar librerías instaladas
installedPackages <- installed.packages()

## Instalar librerias
sapply(packages, function(x){
  if(! x %in% rownames(installedPackages)){ # Instalar solo las que no estan
    print(paste0('Installing ', x))
    tryCatch(install.packages(x), error = function (e) e) # Capturar error si se presenta
  } else {
    print(paste0( x, ' already installed'))
  }
})

## Revisar librerías instaladas
installedPackages <- installed.packages()

## Contrastar si los paquetes requeridos ya están instalados. Deben ser TRUE
packages %in% rownames(installedPackages)


## Opcion para intalar gdalUtils
if (!require(devtools)){
  install.packages("devtools")
}

if (!require(gdalUtils)){
  ## Select option 3: NONE
  devtools:::install_github("gearslaboratory/gdalUtils")## Select option 3: NONE
}

## If there's a problem with gdalUtils
# https://github.com/gearslaboratory/gdalUtils >> CODE (green button) >> download zip
## Unzip the file in C:/path_to_unzziped_folder_pkg or any path
# devtools:::install("C:/path_to_unzziped_folder_pkg")
# Example: devtools:::install("C:/temp/gdalUtils-master") # The unzipped folder

## Cargar librerías
library('raster')
library('rgdal')
library('terra')
library('sf')
library('gdalUtilities')
library('gdalUtils')
library('sp')
library('foreign')


## Eliminar librerias -- en caso que se requiera
# remove.packages( "name_lib_here" )

## Leer funciones
## Cargar funcion de conteo de pixeles
(source("https://raw.githubusercontent.com/gonzalezivan90/SDG15_indicators/main/scripts/tabuleRaster.R"))
## Cargar funcion para encontrar los ejecutables de GDAL en el computador
(source("https://raw.githubusercontent.com/gonzalezivan90/SDG15_indicators/main/scripts/find_gdal.R"))

# Si el anterior genera error por conexión a Github, entonces llamar localmente al achivo:
# El archivo se puede descargar con la opción "Guardar como" desde el navegador, al copiar el link anterior. Es posible que el archivo descarge como R_03_tabuleRaster.R.txt
## source("C:/temp/tabuleRaster.R") ## Cambiar por ruta equivalente, o agregar ".txt" al nombre ## <Dato externo original>
# source("C:/temp/tabuleRaster.R.txt") ## Cambiar por ruta equivalente ## <Dato externo original>
# source("C:/temp/find_gdal.R.txt") ## Cambiar por ruta equivalente ## <Dato externo original>

