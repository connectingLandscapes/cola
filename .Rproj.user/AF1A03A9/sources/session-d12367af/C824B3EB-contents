

# library(gdalUtilities)
# library(gdalUtils)
# library(foreign)
#
# setwd('N:/Mi unidad/connectivity-nasa-USFSIP/01_original-data/peru/')
# rastPath <- 'final_2021.tif'
# origPolPath <- 'ANPNacionalDefinitivas.shp'
# polPath <- 'ap.shp'
# (polLayerName <- tools::file_path_sans_ext(polPath))
#
# ## Copy to a new layer
# #gdalUtilities::ogr2ogr(src_datasource_name = origPolPath, dst_datasource_name = polPath)
#
# gdalUtilities::ogr2ogr(src_datasource_name = origPolPath, dst_datasource_name = polPath,
#                        dialect = 'sqlite', f = "ESRI Shapefile",
#                        sql = paste0('select ST_ExteriorRing(geometry) as geometry from ',
#                                     tools::file_path_sans_ext(origPolPath)))
#
# #ogr2ogr output.shp input.shp -dialect sqlite -f "ESRI Shapefile" -sql "select ST_ExteriorRing(geometry) as geometry from input"
#
# ## Add the table column
# # gdalUtils::gdalinfo(datasetname = polPath, )
# # rgdal::ogrInfo(layer = polPath, sql = paste0('"ALTER TABLE ', polLayerName, ' ADD COLUMN sortid integer"'))
# # ogrInfo(src_datasource_name = origPolPath, dst_datasource_name = polPath)
#
# dbfFile <- gsub('.shp$', '.dbf', x = polPath)
# dbf <- foreign::read.dbf(file = dbfFile, as.is = TRUE)
# dbf$sortid <- 1:nrow(dbf)
# write.dbf(dataframe = dbf, file = dbfFile, factor2char = TRUE)




# ogrinfo train.shp -sql "ALTER TABLE layer ADD COLUMN fid1 integer"
# ogrinfo train.shp -dialect SQLite -sql "UPDATE layer set fid1 = rowid"

# # ogrinfo -sql "ALTER TABLE compiledHUTM18M_v2024 ADD COLUMN km2 NUMERIC(10,15)" compiledHUTM18M_v2024.shp # NO ##### ogrinfo -sql "ALTER TABLE compiledHUTM18M ADD COLUMN km2 VARCHAR(10)" compiledHUTM18M.shp
# ogrinfo -dialect SQLite -sql "UPDATE compiledHUTM18M_v2024 SET km2 = ST_Area(geometry)/1000000" compiledHUTM18M_v2024.shp
# ogr2ogr -sql "SELECT * FROM compiledHUTM18M_v2024 ORDER BY km2 DESC" sortedHUTM18Mv2024.shp compiledHUTM18M_v2024.shp
# NO --- place biggers on top | ogr2ogr -sql "SELECT * FROM compiledHUTM18M_v2024 ORDER BY km2 ASC" sortedHUTM18Mv2024.shp compiledHUTM18M_v2024.shp

# areakm2 <- ogrinfo(datasource_name = paste0(dataPath, '/singleLayers/', metric,'.shp'), dialect = 'SQLite',
#                    sql = paste0("SELECT ", paste0(layerFields, collapse = ', '),", ST_Area(intersection)/1000000 as km2 FROM (SELECT ",
#                                 paste0(layerFields, collapse = ', '), ", ST_Intersection(GEOMETRY, ST_GeomFromText('",
#                                 writeWKT(wkt_pcs, byid = F), "')) AS intersection FROM ", metric, " WHERE ST_Intersects(GEOMETRY, ST_GeomFromText('",
#                                 writeWKT(wkt_pcs, byid = F),"')))") )
# #
