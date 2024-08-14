### COLA web app.
### Ivan Gonzalez - ig299@nau.edu | gonzalezgarzonivan@gmail.com
### Patrick Jantz - Patrick.Jantz@nau.edu | jantzenator@gmail.com

(cat('\n\n >>>> getwd(): ', getwd(), '\n'))

{
  library(cola)

  #library(bit) #
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
  # library(maptools) # deprecated
  # library(raster) # deprecated
  library(RColorBrewer)
  # library(rgdal) # deprecated
  # library(rgeos) # deprecated
  library(rmarkdown)
  library(sf)
  library(shiny)
  library(shinydashboard)
  library(shinydashboardPlus)
  library(shinyjs)
  library(shinyWidgets)
  #library(dashboardthemes)
  library(shinycssloaders)
  library(tidyverse)
  library(shiny)
  library(reshape2)
  library(DT)
  library(tibble)
  library(terra)
  library(viridis)

}


## Init A
{
  ## Initials ---

  os <- Sys.info()[c("sysname")]
  # if ( identical ( unname(Sys.info()[c("sysname", 'nodename')]), c("Windows", 'HP-Z400')) ){
  #   setwd('N:/Mi unidad/IG/server_IG/gedivis')
  #   #setwd('N:/Mi unidad/IG/server_IG/gedivis/')
  # }


  (COLA_DATA_PATH <- Sys.getenv('COLA_DATA_PATH'))
  (dataFolder <<- ifelse( COLA_DATA_PATH != '' & dir.exists(COLA_DATA_PATH),
                          no = paste0(tempdir(), '/'), yes = COLA_DATA_PATH));
  # Keep last slash
  # /data/temp/'; dir.create(dataFolder)

  base::options(scipen = 999)

  (COLA_DSS_UPL_MB <- as.numeric(Sys.getenv('COLA_DSS_UPL_MB')))
  (COLA_DSS_UPL_MB <- as.numeric(
    ifelse( is.numeric(COLA_DSS_UPL_MB) & !is.na(COLA_DSS_UPL_MB),
            yes = COLA_DSS_UPL_MB, no = 250) ))
  base::options('shiny.maxRequestSize' = COLA_DSS_UPL_MB * 1024^2)


  (COLA_VIZ_THREs_PIX <- as.numeric(Sys.getenv('COLA_VIZ_THREs_PIX')))
  (COLA_VIZ_THREs_PIX <- as.numeric(
    ifelse( is.numeric(COLA_VIZ_THREs_PIX) & !is.na(COLA_VIZ_THREs_PIX),
            yes = COLA_VIZ_THREs_PIX, no = 1000000) ))
  base::options('COLA_VIZ_THREs_PIX' = COLA_VIZ_THREs_PIX)


  (COLA_VIZ_RES_NCOL <- as.numeric(Sys.getenv('COLA_VIZ_RES_NCOL')))
  (COLA_VIZ_RES_NCOL <- as.numeric(
    ifelse( is.numeric(COLA_VIZ_RES_NCOL) & !is.na(COLA_VIZ_RES_NCOL),
            yes = COLA_VIZ_RES_NCOL, no = 1000) ))
  base::options('COLA_VIZ_RES_NCOL' = COLA_VIZ_RES_NCOL)


  (COLA_VIZ_RES_NROW <- as.numeric(Sys.getenv('COLA_VIZ_RES_NROW')))
  (COLA_VIZ_RES_NROW <- as.numeric(
    ifelse( is.numeric(COLA_VIZ_RES_NROW) & !is.na(COLA_VIZ_RES_NROW),
            yes = COLA_VIZ_RES_NROW, no = 1000) ))
  base::options('COLA_VIZ_RES_NROW' = COLA_VIZ_RES_NROW)

  #(rootPath <- find.package('cola'))
  (rootPath <- system.file(package = 'cola'))
  path_error <<- '/var/log/shiny-server/'


  source( system.file(package = 'cola', 'app/cola_tools.R') ) # included
  (hs2rs_samp_file <- system.file(package = 'cola', 'sampledata/sampleTif.tif'))
  # file.exists(hs2rs_file)

  # py <- '/home/shiny/anaconda3/envs/cola3/bin/python'
  (py <- Sys.getenv("COLA_PYTHON_PATH"))
  # Sys.getenv("COLA_SCRIPTS_PATH")

  devug <<- TRUE

  (showcasePath <<- base::paste0(rootPath, '/sampledata')); dir.exists(showcasePath)

  uper <- cola::uper ## Unique performance table
  per <- cola::per ## Performance table
  crs_df <- cola::crs_df ## CRS available


  ## Showcase -----
  sh_object <- base::paste0(rootPath, '/docs/showcase/showcase.RData')

  if( base::file.exists(sh_object) ){
    ss <- base::load(sh_object)
  } else {

    sh_hs <- terra::rast(base::paste0(showcasePath, '/HabSui_res375m.tif'))
    sh_sr <- terra::rast(base::paste0(showcasePath, '/Resistance_res375m.tif'))
    sh_pt <- sf::st_read(base::paste0(showcasePath, '/SourcePoints_50.shp'))
    #sh_pt <- spTransform(sh_pt, crs = sf::st_crs('EPSG:4326'))
    sh_pt <- sf::st_transform(sh_pt, crs = sf::st_crs('EPSG:4326'))
    sh_pt[, c('ln', 'lt')] <- st_coordinates(sh_pt)
    sh_crk <- terra::rast(paste0(showcasePath, '/CRK_SP50_DT250k_v1.tif'))
    sh_lcc <- terra::rast(paste0(showcasePath, '/LCC_SP50_DT1mln_CSF5CT5.tif'))


    sh_hs_pal <-leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                      domain= base::range(sh_hs[], na.rm = TRUE) + 0.0,
                                      na.color = "transparent")
    sh_sr_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                      domain = base::range(sh_sr[], na.rm = TRUE)+ 0.0,
                                      na.color = "transparent")
    sh_crk_pal <-leaflet::colorNumeric(palette = "inferno", reverse = TRUE,
                                       domain= base::range(sh_crk[], na.rm = TRUE)+ 0.0,
                                       na.color = "transparent")
    sh_lcc_pal <-leaflet::colorNumeric(palette = "plasma", reverse = TRUE,
                                       domain= base::range(sh_lcc[], na.rm = TRUE)+ 0.0,
                                       na.color = "transparent")

    # addCircleMarkers(lng = cmlng, lat = cmlat, group = "draw")

    ll_sh <- leaflet::leaflet() %>% leaflet::addTiles() %>%
      leaflet::addMeasure( position = "topright",
                           primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
                           activeColor = "#3D535D",completedColor = "#7D4479") %>%
      leaflet::addMiniMap( tiles = leaflet::providers$Esri.WorldStreetMap, toggleDisplay = TRUE) %>%

      addCircleMarkers(lng = sh_pt$ln, lat = sh_pt$lt, group = "Points", radius = 1) %>%

      addRasterImage(sh_hs, colors = sh_hs_pal, opacity = .7,
                     group = "HabitatSuitability", layerId = "HabitatSuitability") %>%
      addRasterImage(sh_sr, colors = sh_sr_pal, opacity = .7,
                     group = "SurfaceResistance", layerId = "SurfaceResistance") %>%
      addRasterImage(sh_lcc, colors = sh_lcc_pal, opacity = .7,
                     group = "Corridors", layerId = "Corridors") %>%
      addRasterImage(sh_crk, colors = sh_crk_pal, opacity = .7,
                     group = "Kernels", layerId = "Kernels") %>%

      # addCircleMarkers(sh_pt, group = "draw") %>%
      # ll_sh %>%
      addLegend(pal =  sh_hs_pal, values = base::range(sh_hs[], na.rm = TRUE),
                group = "HabitatSuitability", layerId = "HabitatSuitability",
                position = 'bottomleft', title = "Hab. suitability")  %>%

      addLegend(pal = sh_crk_pal, values = base::range(sh_crk[], na.rm = TRUE),
                group = "Kernels", layerId = "Kernels",
                position = 'bottomleft', title = "Kernels")  %>%

      addLegend(pal =  sh_sr_pal, values = base::range(sh_sr[], na.rm = TRUE),
                group = "SurfaceResistance", layerId = "SurfaceResistance",
                position = 'bottomleft', title = "Sur. resistance")  %>%

      addLegend(pal = sh_lcc_pal, values = base::range(sh_lcc[], na.rm = TRUE),
                group = "Corridors", layerId = "Corridors",
                position = 'bottomleft', title = "Corridors")  %>%


      leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
      leaflet::addLayersControl(
        baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
        overlayGroups = c('Points',
                          'HabitatSuitability',
                          'SurfaceResistance',
                          'Kernels',
                          'Corridors'
        ),
        options =  leaflet::layersControlOptions(collapsed = FALSE))

    # save(ll_sh, file = sh_object)


    # file.remove(sh_object)
    # rv$hs_sp <- terra::rast(rv$hs)
    # #rng_newtif <- c(newtif@data@min, newtif@data@max)
    # rv$hs_rng <- rng_newtif <- range(rv$hs_sp[], na.rm = TRUE)
    #
    # rv$hs_pal <- hsPal <<- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
    #                                      domain = rng_newtif, na.color = "transparent")
  }
}

## Init B
{
  (sessionID <<- sessionIDgen(folder = TRUE))
  tempFolder <<- paste0(dataFolder, '/', sessionID, '/')
  dir.create(tempFolder)

  (cat(' >>>> COLA_DATA_PATH: ', COLA_DATA_PATH, '\n'))
  (cat(' >>>> tempFolder: ', tempFolder, '\n'))
  (cat(' >>>> R-tempdir(): ', tempdir(), '\n\n'))
}



#  >> SERVER ---------------------------------------------------------------------------
server <- function(input, output, session) {

  getRastVal <- function(clk, grp){
    #clk <- list(lat = 5.9999, lng = 116.89)
    # $lat [1] 5.9999 | $lng     [1] 116.89
    # grp <- c("OpenStreetMap", "draw" )
    # grp <- c("Habitat suitability", "OpenStreetMap", "draw" )

    (avgrp <- length(setdiff(grp, c("OpenStreetMap", "Esri.WorldImagery", "draw")))>0)
    if(!is.null(clk) & !is.null(grp) & avgrp ){

      coords <- cbind(clk$lng, clk$lat)
      pt0 <-  st_sfc(st_point(coords), crs = 4326) # coords = c("x","y")
      #grp <- 'Habitat suitability'
      # rv <- list(hs = '/home/shiny/connecting-landscapes/docs/HS_size5_nd_squared.tif')
      # rv$tiforig <- '/data/temp/XC2024012300322305file1392143e34bb6//out_surface_PF2024012300322605file13921692a78bb.tif'
      # # rv$hs_sp <- terra::rast(rv$h)
      # sapply(X = c(rv$hs_sp, rv$tif_sp, rv$lcc_sp, rv$crk_sp, rv$pritif_sp), FUN = is.null)

      rw2add <- NULL

      if("Habitat suitability" %in%  grp){
        rast_sp <- rv$hs_sp

        rast_prj <- terra::crs(rast_sp, proj = TRUE)
        pt_ <- sf::st_transform(pt0, crs = rast_prj)
        rw2add <- rbind.data.frame(
          rw2add,
          cbind.data.frame(
            lay='Habitat suitability',
            val = as.numeric(terra::extract(rast_sp, st_coordinates(pt_) ) ) )
        )
      }

      if("Surface resistance" %in%  grp){
        rast_sp <- rv$tif_sp

        rast_prj <- terra::crs(rast_sp, proj = TRUE)
        pt_ <- sf::st_transform(pt0, crs = rast_prj)
        rw2add <- rbind.data.frame(
          rw2add,
          cbind.data.frame(lay='Surface resistance',
                           val = as.numeric(terra::extract(rast_sp,
                                                           st_coordinates(pt_) ) ) )
        )
      }

      if("Corridors" %in%  grp){
        rast_sp <- rv$lcc_sp

        rast_prj <- terra::crs(rast_sp, proj = TRUE)
        pt_ <- sf::st_transform(pt0, crs = rast_prj)
        rw2add <- rbind.data.frame(
          rw2add,
          cbind.data.frame(lay='Corridors',
                           val = as.numeric(terra::extract(rast_sp,
                                                           st_coordinates(pt_) ) ) )
        )
      }

      # "Habitat suitability" = rv$hs; "Surface resistance" : rv$tiforig; 'Corridors' rv$lcc; 'Kernels' : rv$crk; 'Prioritization' rv$pritif;
      if("Kernels" %in%  grp){
        rast_sp <- rv$crk_sp

        rast_prj <- terra::crs(rast_sp, proj = TRUE)
        pt_ <- sf::st_transform(pt0, crs = rast_prj)
        rw2add <- rbind.data.frame(
          rw2add,
          cbind.data.frame(lay='Kernels',
                           val = as.numeric(terra::extract(rast_sp,
                                                           st_coordinates(pt_) ) ) )
        )
      }

      if("Prioritization" %in%  grp){
        rast_sp <- rv$pritif_sp

        rast_prj <- terra::crs(rast_sp, proj = TRUE)
        pt_ <- sf::st_transform(pt0, crs = rast_prj)
        rw2add <- rbind.data.frame(
          rw2add,
          cbind.data.frame(lay='Prioritization',
                           val = as.numeric(terra::extract(rast_sp,
                                                           st_coordinates(pt_) ) ) )
        )
      }

      text<-paste0('<strong>', rw2add$lay, ':</strong> ', round(rw2add$val, 3), collapse = '<br>')
    } else {
      text <- 'none'
    }


    return(text)
  }

  observeEvent(input$ll_map_lcc_click,{
    clk <- input$ll_map_lcc_click; grp <- input$ll_map_lcc_groups;
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_lcc") %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_crk_click,{
    clk <- input$ll_map_crk_click; grp <- input$ll_map_crk_groups;
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_crk")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_pri_click,{
    clk <- input$ll_map_pri_click; grp <- input$ll_map_pri_groups;
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_pri")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_edi_click,{
    clk <- input$ll_map_edi_click;
    grp <- input$ll_map_edi_groups;
    #print('--------------------------------------------')
    #print(' > clk '); print(clk)
    #print(' > grp '); print(grp)
    # grp: "OpenStreetMap" "Esri.WorldImagery" "draw" || 'none'
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_edi")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_h2r_click,{
    clk <- input$ll_map_h2r_click; grp <- input$ll_map_h2r_groups; ##
    # print(' ---- clk')
    # print(clk)
    # print(' ---- grp')
    # print(grp)
    text <- getRastVal(clk, grp)
    #print(text)
    if (text != 'none') {
      leafletProxy("ll_map_h2r")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_dist_click,{ ##
    clk <- input$ll_map_dist_click; grp <- input$ll_map_dist_groups; ##
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_dist")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_points_click,{ ##
    clk <- input$ll_map_points_click; grp <- input$ll_map_points_groups; ##
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_points")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  observeEvent(input$ll_map_comp_click,{ ##
    clk <- input$ll_map_comp_click; grp <- input$ll_map_comp_groups; ##
    text <- getRastVal(clk, grp)
    # rastVal <- getRastVal(clk, grp)
    # text<-paste(rastVal$lay, ':', round(rastVal$val, 3), '\n')
    if (text != 'none') {
      leafletProxy("ll_map_points")  %>% clearPopups() %>% addPopups(clk$lng, clk$lat, text)
    }
  })

  # clk <- input$ll_map_pri_click; grp <- input$ll_map_pri_groups;
  # clk <- input$ll_map_crk_click; grp <- input$ll_crk_pri_groups;
  # clk <- input$ll_map_edi_click; grp <- input$ll_edi_pri_groups;
  # clk <- input$ll_map_dist_click; grp <- input$ll_map_dist_groups;
  # clk <- input$ll_map_points_click; grp <- input$ll_map_points_groups;
  # clk <- input$ll_map_h2r_click; grp <- input$ll_map_h2r_groups;

  # observe({
  #   # https://stackoverflow.com/questions/41468538/is-it-possible-to-access-r-leaflet-layer-controls-in-shiny-outside-of-leaflet
  #   # selected_groups <- req(input$my_map_groups)
  #
  #   #https://stackoverflow.com/questions/33575321/how-to-create-an-event-when-pressing-on-a-leaflet-popup-in-r/33575567#33575567
  #   # do whatever ...
  #   #event <- input$myMap_shape_click #Critical Line!!!
  #
  #   # print(list(input$MAPID_click)
  #   # selected_groups <- req(input$my_map_groups)
  #
  #
  #   print(list(input$ll_map_pri_click, input$ll_map_lcc_click,
  #              input$ll_map_crk_click,
  #              input$ll_map_plot_click,
  #              input$ll_map_edi_click, input$ll_map_dist_click,
  #              input$ll_map_points_click, input$ll_map_h2r_click
  #   ))
  #
  #   print(' Groups ' )
  #   print(list(input$ll_map_pri_groups, input$ll_map_lcc_groups,
  #              input$ll_map_crk_groups, input$ll_map_map_groups,
  #              input$ll_map_plot_groups,
  #              input$ll_map_edi_groups, input$ll_map_dist_groups,
  #              input$ll_map_points_groups, input$ll_map_h2r_groups
  #   ))
  #   # selected_groups <- req(input$my_map_groups)
  # })


  output$pdfviewer <- renderText({
    #return(paste('<iframe style="height:600px; width:100%" src="', input$pdfurl, '"></iframe>', sep = ""))
    #pdfManFile <- system.file(package = 'cola', 'docs/pdf.pdf')
    #pdfManFile <- 'docs/pdf.pdf'
    pdfManFile <- 'pdf.pdf' # located in 'www/pdf.pdf'
    return(paste0(
      '<iframe style="height:600px; width:100%" src="',pdfManFile,'"></iframe>'))
    # "/home/shiny/connecting-landscapes/R/pdf_logoA.pdf"
  })

  checkEnv <- function(){
    cat('\n\n\tNames rv: \n')
    cat(names(rv), sep = ' - ')
    cat('\n\n\trv:' )
    print(str(rv))
    # sapply(1:length(rv), function(x){ cat('\n\t', names(rv)[x]); cat(rv[[x]])     })
  }

  resampIfNeeded <- function(rastPath){
    # rastPath <- '/data/temp/ZS2023111311113105file4f6823209882/in_surface_fixed_QK2023111311142905file4f687fa39935.tif'
    # rastPath <- '/data/temp/LI2024011611295205file1a5d191fa71355/in_lcc_CQ2024011611304405file1a5d19779aa294.tif'
    # rastPath <- '/home/shiny/preCanopyClass30mCost.tif'

    r <- terra::rast(rastPath)
    (totpixels <- terra::ncol(r) * terra::nrow(r))
    #names(r)
    #rasname <- r@pnt@.xData$filenames()
    (resamPath <- gsub(x = rastPath,
                       '.tif$', '_resam.tif'))
    if(totpixels > options('COLA_VIZ_THREs_PIX')){
      #if(file.exists()){
      if (!file.exists(resamPath)){

        gdalUtilities::gdalwarp(srcfile = rastPath,
                                ts = c(options('COLA_VIZ_RES_NCOL'), options('COLA_VIZ_RES_NROW')),
                                dstfile = resamPath)
        print(paste(' ---- >>>> Resampling to ',
                    options('COLA_VIZ_RES_NCOL'),' - ', options('COLA_VIZ_RES_NROW')))
      }
      #}
      return(resamPath)
    } else {
      return(rastPath)
    }
  }

  #test <- burnShp(polDraw, burnval, rastPath, rastCRS)
  burnShp <- function(polPath, burnval = 'val2burn',
                      rastPath, rastCRS = NA,  att = FALSE, lineBuffW = 1){
    # rastPath <- '/data/temp/XZ2024041911393405file9c152374e9a2/in_edit_DJ2024041911410705file9c15429f450d.tif'
    # rast <-terra::rast(rastPath)
    # burnval = -10
    # (load(file = '/data/tempR/draw.RData')) # polDraw

    #if( burnval != 0 & is.numeric(burnval) & !is.na(burnval) ){

    #(polPath <- gsub(x = rastPath, '.tif$', '_pol.shp'))
    (rasterizedPath <- gsub(x = rastPath, '.tif$', '_rasterized.tif'))

    file.copy(rastPath, rasterizedPath, overwrite = TRUE)

    #load(file = '/data/temp/2lines.RData') # polDraw #load(file = '/data/temp/4geom.RData') # polDraw
    #str(polDraw) #polDraw$type # FeatureCollection

    rt <- terra::rast(rastPath)
    rastRes <- res(rt)

    if( is.na(rastCRS)){
      ## rastPath <- '/home/shiny/connecting-landscapes/docs/HS_size5_nd_squared.tif'
      #gi <- rgdal::GDALinfo(rastPath)
      #gi2 <- sf::gdal_crs(rastPath)
      #gi <- strsplit(x = gdalUtilities::gdalinfo(rastPath), '\n')
      #prj <- attr(x = gi, "projection")
      prj <- terra::crs(rt, proj = TRUE)
      (rastCRS <- st_crs(rt))
      # if(!is.na(prj)){
      #   pol2save@proj4string@projargs <- prj
      # }
      #ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:NEW_EPSG_NUMBER -s_srs EPSG:OLD_EPSG_NUMBER output.shp input.shp
      #EPSG:4326
    }

    gdalUtilities::gdal_rasterize(
      src_datasource = polPath,
      at = att,
      dst_filename = rasterizedPath,
      add = TRUE,
      a = 'val2burn') #as.numeric(burnval)

    #file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)
    # rasteri <- gdalUtils::gdal_rasterize(src_datasource = polPath, at = T,
    #                                      dst_filename = rastPath,
    #                                      add = TRUE, a = 'val2burn')
    # plot(raster(rastPath))
    # plot(raster(rasterizedPath), main = 'Rasterized')
    # plot(sf::read_sf(polPath), add = TRUE)
    # file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)

    return(rasterizedPath)
    # } else {
    #   return(NA)
    # }
  }


  #test <- replaceRastShp(polDraw, burnval, rastPath, rastCRS)
  replaceRastShp <- function(polPath, burnval = 'val2burn', rastPath,
                             att = FALSE, rastCRS = NA, lineBuffW = 1){

    #if( burnval != 0 & is.numeric(burnval) & !is.na(burnval) ){
    ## Polygon to write
    #(polPath <- gsub(x = rastPath, '.tif$', '_pol.shp'))
    ## Raster with new features
    (rasterizedPath <- gsub(x = rastPath, '.tif$', '_rasterized2replace.tif'))
    ## Raster to create
    (replacedPath <- gsub(x = rastPath, '.tif$', '_replaced.tif'))


    rt <- terra::rast(rastPath)
    (rastRes <- res(rt))
    if( is.na(rastCRS)){
      prj <- terra::crs(rt, proj = TRUE)
      (rastCRS <- st_crs(rt))
    }

    # pol2Rast <- draws2Features(polDraw, distLineBuf = min(rastRes), rastCRS = rastCRS)
    # # pol2Rastx <- st_sf(data.frame(a = 1:length(pol2Rast), pol2Rast))
    # pol2Rastx <- st_as_sf( pol2Rast)
    # # plot(pol2Rastx, add = TRUE, border = 'blue', col = NA)
    #
    # pol2Rastx$val2burn <- as.numeric(burnval)
    # sf::st_write( obj = pol2Rastx, dsn = dirname(polPath),
    #               layer = tools::file_path_sans_ext(basename(polPath)),
    #               driver = 'ESRI Shapefile',
    #               append = FALSE,
    #               overwrite_layer = TRUE)


    ## Rasterize polygons
    rastExtent <- as.character(as.vector(terra::ext(rt))[c('xmin', 'ymin', 'xmax', 'ymax')])
    gdalUtilities::gdal_rasterize(
      te = rastExtent, # -te <xmin> <ymin> <xmax> <ymax>
      # tr = rastRes,
      ts = c(terra::ncol(rt), terra::nrow(rt)), # <width> <height>
      src_datasource = polPath,
      at = att,
      dst_filename = rasterizedPath,
      add = FALSE,
      a = 'val2burn') #as.numeric(burnval)

    #rft <- rast(rasterizedPath); plot(rft)

    (cmdCalc <- paste0('gdal_calc.py --overwrite ',
                       ' -A ', rastPath, # original
                       ' -B ', rasterizedPath,
                       ' --calc "((B == 0 ) * A ) + ((B != 0 ) * B )"',
                       ' --outfile ', replacedPath))

    runCMD <- system(cmdCalc, intern = TRUE)

    # plot(rast(rastPath))
    # plot(rast(rasterizedPath), main = 'Rasterized')
    # plot(sf::read_sf(polPath), add = TRUE, col = NA, border = 'red')
    # plot(rast(replacedPath), main = 'Replaced')
    # file.remove(rasterizedPath); file.copy(rastPath, rasterizedPath, overwrite = TRUE)

    return(replacedPath)
    #} else { return(NA) }
  }

  makeLL <- function(){
    # https://rstudio.github.io/leaflet/morefeatures.html
    ll0 <- leaflet::leaflet()
    grps <- c()
    {
      # pdebug(devug = devug, sep = '\n-', pre = '\n Make LL\n',
      #        'rv$hs', 'file.exists(rv$hs)', 'rv$hsready',
      #        'rv$tif', 'file.exists(rv$tif)', 'rv$tifready'
      #        )

      ## Debug
      # rv <- list(path = '/data/temp/OK2024011614464605file3f679747355/')
      # rv$tif <- paste0(rv$path, 'in_lcc_fixedPQ2024011614471405file3f66bce19da.tif')
      # rv$crk <- paste0(rv$path, 'out_crk_PS2024011614474105file3f64ad67212.tif')
      # rv$prishp <- paste0(rv$path, 'out_pri_VL2024011614475305file3f6175e5eb2.shp')
      # rv$pritif <- paste0(rv$path, 'out_pri_VL2024011614475305file3f6175e5eb2.tif')

      hs_pal_name <-  "viridis"
      sr_pal_name <- "magma"
      tif_pal_name <- "viridis"
      crk_pal_name <- "inferno"
      lcc_pal_name <- "plasma"

      # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
      #leaflet::colorNumeric(palette = "magma", reverse = TRUE,
      #                         domain= base::range(sh_sr[], na.rm = TRUE)+ 0.0,
      #                         na.color = "transparent")

      if( (rv$hsready) ){
        #pdebug(devug=devug,pre='\n MakeLL - HS\n',sep='\n-','rv$hs_pal','rv$hs_rng')

        grps <- c(grps, "Habitat suitability")

        if (rv$hs0 != rv$hs){
          rv$hs2s <- resampIfNeeded(rv$hs)
          rv$hs0 <- rv$hs
        }
        rv$hs2s_sp <- terra::rast(rv$hs2s)
        rv$hs_rng2 <- getMxMn(rv$hs2s)+ 0.0 # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
        rv$hs_pal2 <- leaflet::colorNumeric(palette = hs_pal_name, reverse = TRUE,
                                            domain = rv$hs_rng2 + 0.0,
                                            na.color = "transparent")

        # pdebug(devug = devug, sep = '\n-', pre = '-', '(rv$hs)', 'rv$hs2s')

        ll0 <- ll0 %>% addRasterImage(x = rv$hs2s_sp, colors = rv$hs_pal2,
                                      opacity = .7,
                                      group = "Habitat suitability", layerId = "HabitatSuitability") %>%
          addLegend(pal =  rv$hs_pal2, values = rv$hs_rng2,
                    group = "Habitat suitability", layerId = "Habitat suitability",
                    position = 'bottomleft', title = "Suitability"#, opacity = .3
                    #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          )

        #inPoints <<- input$in_points_ly
        #pdebug(devug = devug, pre = '\n', sep = '\n-',  '(rv$hs)', 'rv$hs2s', 'rv$point_choices', "input$in_points_ly", 'inPoints')
        rv$point_choices <<- unique(c(rv$point_choices, 'HabitatSuitability'))
        updateSelectizeInput(session, inputId = "in_points_ly",
                             selected = 'HabitatSuitability',
                             choices = rv$point_choices,
                             server = TRUE)
        # pdebug(devug = devug, pre = '\n', sep = '\n-', 'rv$point_choices')
      }


      if(rv$tifready){
        grps <- c(grps, "Surface resistance")

        # pdebug(devug=devug,pre='\n  MakeLL - TIF', sep='\n','rv$tiforig','rv$tif0', 'rv$tif0 != rv$tiforig', 'rv$tif2s', 'rv$tif_rng2')

        if (!(rv$tif0 == rv$tiforig)){
          rv$tif2s <- resampIfNeeded(rv$tiforig)
          rv$tif0 <- rv$tiforig
        }
        rv$tif2s_sp <- terra::rast(rv$tif2s)
        rv$tif_rng2 <- getMxMn(rv$tif2s) + 0.0 # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
        #
        rv$tif_rng2 <- c(max(0, rv$tif_rng2[1]), rv$tif_rng2[2])

        rv$tif_pal2 <-leaflet::colorNumeric(palette = tif_pal_name, reverse = TRUE,
                                            domain = rv$tif_rng2 + 0.0,
                                            na.color = "transparent")
        # pdebug(devug=devug,pre='\n  MakeLL - TIF', sep='\n','rv$tiforig','rv$tif0', ' rv$tif0 != rv$tiforig', 'rv$tif2s', 'rv$tif_rng2')

        ll0 <- ll0 %>%
          addRasterImage(rv$tif2s_sp, colors = rv$tif_pal2,
                         opacity = .7,
                         group = "Surface resistance",
                         layerId = "Surfaceresistance") %>%

          addLegend(pal =  rv$tif_pal2, values = rv$tif_rng2,
                    group = "Surface resistance",
                    layerId = "Surfaceresistance",
                    position = 'bottomleft', title = "Resistance"#, opacity = .3
                    #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          )

        #pdebug(devug = devug, pre = '\n', sep = '\n-', 'rv$tif', 'rv$point_choices', "input$in_points_ly")
        rv$point_choices <- unique(c(rv$point_choices, 'SurfaceResistance'))
        updateSelectizeInput(session, "in_points_ly",
                             choices = rv$point_choices,
                             selected = 'SurfaceResistance',
                             server = TRUE)
        # pdebug(devug = devug, pre = '\n', sep = '\n-', 'rv$point_choices')
        #pdebug(devug = devug, pre = '\n', sep = '\n-', 'rv$tif', 'rv$point_choices', "input$in_points_ly")
      }

      if((rv$crkready)){
        grps <- c(grps, 'Kernels')

        if (rv$crk0 != rv$crk){
          rv$crk2s <- resampIfNeeded(rv$crk)
          rv$crk0 <- rv$crk
        }
        rv$crk2s_sp <- terra::rast(rv$crk2s)
        rv$crk_rng2 <- getMxMn(rv$crk2s) + 0.000 # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
        print(rv$crk_rng2)
        rv$crk_pal2 <-leaflet::colorNumeric(palette = crk_pal_name, reverse = TRUE,
                                            domain = rv$crk_rng2 + 0.001,
                                            na.color = "transparent")


        ll0 <- ll0 %>% addRasterImage(rv$crk2s_sp, colors = rv$crk_pal2, opacity = .7,
                                      group = "Kernels", layerId = "Kernels") %>%
          addLegend(pal =  rv$crk_pal2, values = rv$crk_rng2,
                    layerId = "Kernels", group = "Kernels",
                    position = 'bottomleft', title = "Kernels"#, opacity = .3
                    #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          )
      }

      if((rv$lccready)){

        # pdebug(devug=devug,pre='\n MakeLL - LCC ',sep='\n','rv$lcc')
        # rv <- list(lcc = '/data/temp/YM2024011518570905file1a4cf94fe3cee0//out_lcc_XL2024011518571705file1a4cf962d6830f.tif')

        grps <- c(grps, 'Corridors')

        if (rv$lcc0 != rv$lcc){
          rv$lcc2s <- resampIfNeeded(rv$lcc)
          rv$lcc0 <- rv$lcc
        }
        rv$lcc2s_sp <- terra::rast(rv$lcc2s)
        rv$lcc_rng2 <- getMxMn(rv$lcc2s) + 0.001 # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
        rv$lcc_pal2 <-leaflet::colorNumeric(palette = lcc_pal_name, reverse = TRUE,
                                            domain = rv$lcc_rng2 + 0.0,
                                            na.color = "transparent")

        ll0 <- ll0 %>% addRasterImage(rv$lcc2s_sp, colors = rv$lcc_pal2, opacity = .7,
                                      group = "Corridors", layerId = "Corridors") %>%
          addLegend(pal =  rv$lcc_pal2, values = rv$lcc_rng2,
                    layerId = "Corridors", group = "Corridors",
                    position = 'bottomleft', title = "Corridors"#, opacity = .3
                    #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          )
      }


      if((rv$priready)){

        #rv$pritif <- out_crk
        #rv$pritif_sp <- terra::rast(rv$pritif); #plot(newtif)
        #rv$prishp_sp <- sf::read_sf(rv$prishp); #plot(newtif)

        grps <- c(grps, 'Prioritization')

        if (rv$pritif0 != rv$pritif){
          rv$pritif2s <- resampIfNeeded(rv$pritif)
          rv$pritif0 <- rv$pritif
        }
        rv$pritif2s_sp <- terra::rast(rv$pritif2s)
        rv$pritif_rng2 <- getMxMn(rv$pritif2s) + 0.001 # hs = "viridis" | sr "magma" | crk "inferno" | lcc "plasma"
        rv$pritif_pal2 <- leaflet::colorNumeric(palette = 'Blues', reverse = FALSE,
                                                domain = rv$pritif_rng2 + 0.0,
                                                na.color = "transparent")

        rv$prishp_sp$ID <- order(rv$prishp_sp$cp1)
        rv$prishp_wgs <- st_transform(rv$prishp_sp, '+proj=longlat +datum=WGS84')
        patchTif <- (gsub('out_pri_', 'out_pri_patch_', rv$pritif))
        patchTifR <- rast(patchTif)

        patchpal <- colorBin(palette = "Dark2", domain = c(1, max(getMxMn(patchTif))), bins = c( 1, 43),
          reverse = FALSE, na.color = "transparent")

        ll0 <- ll0 %>%
          addRasterImage(rv$pritif_sp, colors = rv$pritif_pal2, opacity = .5,
                                      group = "Prioritization", layerId = "Priorit. corr") %>%
           addRasterImage(patchTifR, colors = patchpal,
                          opacity = .5, group = "Prioritization", layerId = "Priorit. patch") %>%
          addCircleMarkers(data = rv$prishp_wgs, label = ~ID, group = 'Prioritization',  radius = 5, color = 'red')  %>%
          addLegend(pal =  rv$pritif_pal2, values = rv$pritif_rng2,
                    layerId = "Prioritization", group = "Prioritization",
                    position = 'bottomleft', title = "Prioritization"#, opacity = .3
                    #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          )
      }

      if((rv$ptsready)){
        grps <- c(grps, 'Points')
        #rv$pts_sp$ID <- 1:nrow(rv$pts_sp)
        ll0 <- ll0 %>%  addCircleMarkers(data = rv$pts_sp, label = ~sortID,
                                         group = 'Points',  radius = 5)
      }
    }


    grps <<- rev(grps)
    if( length(grps) != 0){

      ll0 <- ll0 %>%leaflet::addTiles() %>% #clearBounds() %>%
        leaflet::addLayersControl(baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
                                  overlayGroups = rev(grps),
                                  options =  leaflet::layersControlOptions(collapsed = FALSE)) %>%
        leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
        leaflet::addMeasure( position = "topright",
                             primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
                             activeColor = "#3D535D",completedColor = "#7D4479")  %>%
        leaflet::addMiniMap( tiles = leaflet::providers$Esri.WorldStreetMap, toggleDisplay = TRUE)


      # grps <- c("Habitat suitability", "Surface resistance", "Corridors", "Kernels", 'Points')
      # grps <- c("Surface resistance", "Corridors", "Kernels", 'Points')
      # grps <- c("Surface resistance", "Kernels", 'Points')

      (grpTIF <<- rev(grps[grps != 'Points'])) # Layers to hide
      if(length(grpTIF) > 1){
        (hideGr <<- grpTIF[2:(length(grpTIF))])
        ll0 <- ll0 %>% hideGroup(hideGr)
      }
      #  pdebug(devug = TRUE, sep = '\n', pre = '\n', 'grps', 'grpTIF', 'hideGr')

    }

    if (FALSE){

      # if(TRUE){
      #   setwd('/data/temp/DE2023090812244005file12665b147993/')
      #   pts_sp = sf::read_sf('out_simpts_TV2023090812282705file126638a90d5c.shp') # path
      #   pts_sp@data$ID <- 1:nrow(pts_sp)
      #   pts_sp <- spTransform(pts_sp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
      #   rv <<- list(
      #     hs_sp = terra::rast('in_surface_TV2023090812282705file126638a90d5c.tif'), # path
      #     tif_sp =terra::rast('out_surface_TV2023090812282705file126638a90d5c.tif'), # path
      #     lcc_sp =terra::rast('out_lcc_LV2023090813140005file12667acee23d.tif'), # spatial object
      #     crk_sp =terra::rast('out_crk_PM2023090813154705file126665ed88ec.tif') # spatial object
      #   )

      #   setwd('/data/temp/LG2023100214170305file84944f27902c/')
      # rv <<- list( edi_sp =terra::rast('in_edit_fixed_KN2023100214171605file849445dd4363c.tif'))

      #   rv$pts_sp <- pts_sp
      #   #rv$hs_sp <- 100 - rv$tif_sp
      #
      #   rv$hs_rng <- range(rv$tif_sp[], na.rm = TRUE)
      #   rv$tif_rng <- range(rv$tif_sp[], na.rm = TRUE)
      #   rv$lcc_rng <- range(rv$lcc_sp[], na.rm = TRUE)
      #   rv$crk_rng <- range(rv$crk_sp[], na.rm = TRUE)
      #
      #   ## Color pal
      #   rv$hs_pal <-leaflet::colorNumeric(palette = "Accent", reverse = TRUE,
      #                             domain = rv$hs_rng, na.color = "transparent")
      #   rv$tif_pal <-leaflet::colorNumeric(palette = "RdYlBu", reverse = TRUE,
      #                              domain = rv$tif_rng, na.color = "transparent")
      #   rv$lcc_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
      #                              domain = rv$lcc_rng+0.0, na.color = "transparent")
      #   rv$crk_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
      #                              domain = rv$crk_rng+0.0, na.color = "transparent")
      #   rv$hs <- rv$tif <- rv$lcc <- rv$crk <- rv$pts <- 1
      #   makeLL(rv)
      # }
    }

    rv$ll <- ll0
    updateLL(ll0)
    return(ll0)
  } # llz <- makeLL()

  updateLL <- function(ll){
    output$ll_map_cdp <- output$ll_map_pri <- output$ll_map_lcc <-
      output$ll_map_crk <- output$ll_map_map <- output$ll_map_plot <-
      output$ll_map_edi <- output$ll_map_dist <-
      output$ll_map_points <- output$ll_map_h2r <- leaflet::renderLeaflet({
        ll
      })

    ## add draw tools to scenarios and comparisson
    output$ll_map_com <- output$ll_map_edi <- leaflet::renderLeaflet({ ll %>%
        leaflet.extras::addDrawToolbar(singleFeature = FALSE,
                                       targetGroup='draw', polylineOptions = T,
                                       rectangleOptions = T, circleOptions = T,
                                       markerOptions = T, circleMarkerOptions = F,
                                       editOptions = leaflet.extras::editToolbarOptions())
    })
  }

  output$ll_map_show <- leaflet::renderLeaflet({ ll_sh })
  output$ll_map_showPriv <- leaflet::renderLeaflet({  llmap })

  updateParamsTEXT <- function(params_txt = params_txt,
                               start = FALSE, pprint = FALSE,
                               hs = FALSE, sr = FALSE,
                               pts = FALSE, dst = FALSE,
                               lcc = FALSE, crk = FALSE){
    # hs = FALSE; sr = FALSE; pts = FALSE; lcc = FALSE; crk = FALSE

    color_darkgreen_hex <- '#006400'
    color_darkred_hex <- '#FF0000'

    readyMsg <- paste0("<font color=\"", color_darkgreen_hex,"\"><b> Ready</b></font>")
    notReadyMsg <- paste0("<font color=\"", color_darkred_hex,"\"><b> Not ready</b></font>")


    # if( pprint){
    #   print('======== Params Before')
    #   print(rbind(do.call(cbind, params_txt),do.call(cbind, params_txt_loc)))
    # }

    if (start){
      output$out_par_surA <- renderText({ paste0("Input: Habitat suitability", notReadyMsg)})
      output$out_par_ptsA <- renderText({paste0("Input: Raster ", notReadyMsg)})
      output$out_par_prioA <<- output$out_par_ediA <<- output$out_par_lccB <<- output$out_par_crkB <<-
        renderText({paste0("Input: Surface resistance", notReadyMsg)})
      output$out_par_lccA <- output$out_par_crkA <- renderText({ paste0("Input: Points", notReadyMsg)})
      output$out_par_prioB <- renderText({ paste0("Input: Corridors", notReadyMsg)})
      output$out_par_prioC <- renderText({ paste0("Input: Kernels", notReadyMsg)})


      output$out_par_distA <- renderText({ paste0("Input: Surface resistance", notReadyMsg)})
      output$out_par_distB <- renderText({ paste0("Input: Points", notReadyMsg)})
      output$out_par_distC <- renderText({ paste0("Output: Distance", notReadyMsg)})

      output$out_par_cdpoA <- renderText({ paste0("Input: Points", notReadyMsg)})
      output$out_par_cdpoB <- renderText({ paste0("Input: Distance matrix", notReadyMsg)})
    }

    params_txt_loc <- params_txt
    hs <<- params_txt_loc$hs <- (hs | params_txt_loc$hs)
    sr <<- params_txt_loc$sr <- (sr | params_txt_loc$sr)
    params_txt <<- params_txt_loc

    ## Turn labels to 'ready'
    if( hs  ){
      output$out_par_surA <- renderText({ paste0("Input: Habitat suitability", readyMsg)})
      output$out_par_ptsA <- renderText({paste0("Input: Raster ", readyMsg)})
    }

    if( sr ){
      output$out_par_distA <<- output$out_par_prioA <<- output$out_par_ediA <<- output$out_par_lccB <<- output$out_par_crkB <<-
        renderText({paste0("Input: Surface resistance", readyMsg)})
      output$out_par_ptsA <- renderText({paste0("Input: Raster ", readyMsg)})
    }

    if( pts ){
      output$out_par_cdpoA <- output$out_par_distB <- output$out_par_lccA <- output$out_par_crkA <- renderText({ paste0("Input: Points", readyMsg)})
    }

    if( dst ){
      output$out_par_cdpoB <- renderText({ paste0("Input: Distance", readyMsg)})
      output$out_par_distC <- renderText({ paste0("Output: Distance", readyMsg)})
    }

    if( lcc ){
      output$out_par_prioB <- renderText({ paste0("Input: Corridors", readyMsg)})

    }

    if( crk ){
      output$out_par_prioC <- renderText({ paste0("Input: Kernels", readyMsg)})
    }

    return(params_txt_loc)
    # paste0("Input: Habita
  }

  output$vout_com <- renderText({isolate('Waiting for inputs')})

  updateVTEXT <- function(txt, devug = FALSE){
    if(devug){print(txt)}
    output$vout_cdp <-  output$vout_h2r <- output$vout_pri <- output$vout_points <-
      output$vout_dist <- output$vout_crk <- output$vout_cdpop <-
      output$vout_lcc <- output$vout_edi <- renderText({isolate(txt)})
  }

  output$distPlot <- renderPlot({
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    hist(x, breaks = bins, col = "darkgray", border = "white", main = input$title)
  })

  DEFAULT_VALS_STARTER <- function(x){x}
  rv <<- reactiveValues(
    layersList = NULL,
    sessionID = sessionID,  tempFolder = tempFolder,

    data = NULL, orig = NULL,
    cdpopRun = NULL, out_cdpop_files = c(''),

    log = paste0('Your session ID: ', sessionID, '\nWaiting for inputs ... '),

    llmap0 = leaflet::leaflet() %>%leaflet::addTiles(),

    surfmap = NULL, pointsmap = NULL,
    distmap = NULL, distrast = NULL, distshp = NULL,

    point_choices = NULL,

    ## Check readyness ...
    hsready = FALSE,
    editready = FALSE,
    tifready = FALSE,
    ptsready = FALSE,
    cdmready = FALSE,
    lccready = FALSE,
    crkready = FALSE,
    priready = FALSE,
    sceready = FALSE,
    comready = FALSE,

    ## Spatial file paths
    hs = NULL, # path
    tiforig = NULL, # path
    tif = NULL, # path
    edi = NULL, # path
    pts = NULL, # path
    shp = NULL, # spatial object
    lcc = NULL, # spatial object
    crk = NULL, # spatial object
    pritif = NULL, # spatial object
    prishp = NULL, # spatial object
    cdm = NULL, # csv
    com = NULL, # csv

    ## Spatial file 2 show (2s) paths
    hs2s = NULL, # path
    tiforig2s = NULL, # path
    tif2s = NULL, # path
    edi2s = NULL, # path
    pts2s = NULL, # path
    shp2s = NULL, # path
    lcc2s = NULL, # path
    crk2s = NULL, # path
    pritif2s = NULL, # path
    cdm2s = NULL, # csv

    ## Original spatial file 2 show (2s) paths -- for avoid overt
    hs0 = "", # path
    tiforig0 = "", # path
    tif0 = "", # path
    edi0 = "", # path
    pts0 = "", # path
    shp0 = "", # path
    lcc0 = "", # path
    crk0 = "", # path
    pritif0 = "", # path
    com0 = "", # path
    cdm0 = "", # csv

    ## Spatial objects
    hs_sp = NULL, # spatial object
    edi_sp = NULL, # spatial object
    tif_sp = NULL, # spatial object
    pts_sp = NULL, # spatial object
    lcc_sp = NULL, # spatial object
    crk_sp = NULL, # spatial object
    pritif_sp = NULL, # spatial object
    cdm_sp = NULL, # csv

    ## Spatial objects to show
    hs2s_sp = NULL, # spatial object
    edi2s_sp = NULL, # spatial object
    tif2s_sp = NULL, # spatial object
    pts2s_sp = NULL, # spatial object
    lcc2s_sp = NULL, # spatial object
    crk2s_sp = NULL, # spatial object
    pritif2s_sp = NULL, # spatial object
    cdm2s_sp = NULL, # csv

    ## Color pal
    hs_pal = NULL, # path
    tif_pal = NULL, # path
    shp_pal = NULL, # spatial object
    lcc_pal = NULL, # spatial object
    pritif_pal = NULL, # spatial object
    crk_pal = NULL, # spatial object

    # SessionIDs
    inLccSessID = NULL,
    incrkSessID = NULL,
    inDistSessID = NULL,
    inPointsSessID = NULL,
    inSurSessID = NULL,
    inProjSessID = NULL,

    newtifPath = NULL,
    newtifPath_pts = NULL,
    newtifPath_dist = NULL,
    newshpPath_dist = NULL,

    tifpathlcc = NULL,
    tifpathcrk = NULL,
    tifpathlccfix = NULL,
    tifpathcrkfix = NULL,

    # Projection
    # inProjSessID
    crs_db = NULL,
    crs_pts_temp = "",
    crs_pts = "",
    crs_tif_temp = "",
    crs_tif = "",

    last = NULL)

  # rv <- list()

  rv$sessionID <- sessionID
  rv$tempFolder <- tempFolder

  llmap <<- leaflet::leaflet() %>%leaflet::addTiles() %>%
    leaflet::addLayersControl(baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
                              options =  leaflet::layersControlOptions(collapsed = FALSE)) %>%
    leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
    leaflet::addMeasure( position = "topright",
                         primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
                         activeColor = "#3D535D",completedColor = "#7D4479") %>%
    # leaflet.extras::addDrawToolbar(singleFeature = FALSE,
    #                                targetGroup='draw', polylineOptions = TRUE,
    #                                rectangleOptions = TRUE, circleOptions = TRUE,
    #                                markerOptions = TRUE, circleMarkerOptions = F,
    #                                editOptions = leaflet.extras::editToolbarOptions()) %>%
    leaflet::addMiniMap( tiles = leaflet::providers$Esri.WorldStreetMap, toggleDisplay = TRUE)


  rv$llmap0 <- rv$llmap <- llmap

  # rv$llmap rv$hsready rv$tifready rv$ptsready
  # ll_map_corr lcc vout_corr in_lcc_3 4 5
  # ll_map_crk crk vout_crk in_crk_3 4
  # ll_map_plot ll_map_map
  # rv$hsready = FALSE
  # rv$tifready = FALSE
  # rv$ptsready = FALSE
  output$ll_coord <- leaflet::renderLeaflet({ rv$llmap })
  updateLL(rv$llmap)
  updateVTEXT(rv$log) #'Waiting for inputs')
  params_txt <- list(hs = FALSE, sr = FALSE, pts = FALSE, lcc = FALSE, crk = FALSE)
  params_txt <- updateParamsTEXT(params_txt = params_txt, start = TRUE)

  output$hccdpop1 <- output$hccdpop2 <- output$hccdpop3 <- highcharter::renderHighchart({ highchart() })
  output$hccomp1 <- output$hccomp2 <- highcharter::renderHighchart({ highchart() })


  # output$ll_map_lcc <- output$ll_map_crk <- output$ll_map_map <- output$ll_map_plot <-
  #   output$ll_map_dist <- output$ll_map_points <- output$ll_map_h2r <- leaflet::renderLeaflet({
  #   rv$llmap
  # })

  # vtext <- "Waiting for the habitat sutiability TIF"
  # output$vout <- renderText({isolate(vtext)})
  #
  # vout_points <- "Waiting for the surface resistance TIF"
  # output$vout_points <- renderText({isolate(vout_points)})
  #
  # vout_dist <- "Waiting for the surface resistance TIF and points"
  # output$vout_dist <- renderText({isolate(vout_dist)})
  #
  # vout_crk <- "Waiting for the surface resistance TIF and points"
  # output$vout_crk <- renderText({isolate(vout_crk)})
  #
  # vout_lcc <- "Waiting for the surface resistance TIF and points"
  # output$vout_lcc <- renderText({isolate(vout_lcc)})


  # SRV FIX PROJ  ------------------
  # in_uncrs_pts in_uncrs_tif coo_pts coo_pts ll_coord | pts_uncrs
  # selectizeInput("sel_crs", "Select", choices = NULL), #

  updateSelectizeInput(session, inputId = 'sel_crs',
                       choices = c("", crs_df$label),
                       selected = NA, server = TRUE)

  updateSelectizeInput(session, inputId = 'sel_crs2',
                       choices = c("", crs_df$label),
                       selected = NA, server = TRUE)

  observeEvent(input$in_uncrs_pts, {
    invisible(suppressWarnings(
      tryCatch(file.remove(c(rv$crs_pts_orig, rv$crs_pts)),
               error = function(e) NULL)))
    tempID <- sessionIDgen()
    rv$crs_pts_orig <- input$in_uncrs_pts$datapath
    rv$crs_pts <- file.path(tempFolder, paste0('/in_uncrs_',
                                               tempID, '_', basename(rv$crs_pts_orig)))
    pdebug(devug=devug,pre='\n\t Load uncrs PTS\n', sep='\n','rv$crs_pts_orig', 'rv$crs_pts')
    file.copy(rv$crs_pts_orig, rv$crs_pts);
    #rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______

    #rv <- list(crs_pts = '~/small_test_10pts.xy')
    rv$pts_uncrs <- read.csv(rv$crs_pts)
    rv$pts_uncrs$lng <- rv$pts_uncrs[, grep('[xX]',colnames(rv$pts_uncrs))[1]]
    rv$pts_uncrs$lat <- rv$pts_uncrs[, grep('[yY]',colnames(rv$pts_uncrs))[1]]
    rv$pts_uncrs <- st_as_sf(rv$pts_uncrs, coords=c("lng","lat"))
    #coordinates(rv$pts_uncrs) =~ lng + lat
    rv$pts_uncrs_extent <- as.polygons(terra::ext(rv$pts_uncrs))
    #plot(rv$pts_uncrs)
    #plot(pts_uncrs_extent, add = TRUE)
  })

  observeEvent(input$in_uncrs_tif, {
    invisible(suppressWarnings(
      tryCatch(file.remove(c(rv$crs_tif_orig, rv$crs_tif)),
               error = function(e) NULL)))
    tempID <- sessionIDgen()
    rv$crs_tif_orig <- input$in_uncrs_tif$datapath
    rv$crs_tif <- file.path(tempFolder, paste0('/in_uncrs_',
                                               tempID, '_', basename(rv$crs_tif_orig)))
    pdebug(devug=devug,sep='\n\t Load uncrs TIF', pre='','rv$crs_tif_orig', 'rv$crs_tif')
    file.copy(rv$crs_tif_orig, rv$crs_tif);

    #rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______
    rv$tif_uncrs <- terra::rast(rv$crs_tif)
    rv$tif_uncrs_extent <- as.polygons(terra::ext(rv$tif_uncrs))
    #plot(rv$pts_uncrs)
  })

  # rv <- list(crs_pts = '/data/temp/TN2023090900452605file23e6f76efff8b///in_uncrs_IP2023090900454205file23e6f70719c16_0.xy',
  #           crs_tif = '/data/temp/TN2023090900452605file23e6f76efff8b///in_uncrs_OG2023090900453505file23e6f47071d2d_0.rsg')

  # # crs_tif_temp = "", # # crs_tif = "", # input$sel_crs # "" # input$in_uncrs_pts # input$in_uncrs_tif

  observe({
    ## Create an ID
    if(input$sel_crs != '' | input$sel_crs2 != '') {
      if(is.null(rv$inProjSessID)){
        rv$inProjSessID <- sessionIDgen()
      }
      change <<- 0

      if(input$sel_crs != '') {
        crs_selected <<- crs_df$crs_code[crs_df$label %in% input$sel_crs]
      }
      if(input$sel_crs2 != '') {
        crs_selected2 <<- crs_df$crs_code[crs_df$label %in% input$sel_crs2]
      }
      # (crs_selected <- crs_df$crs_code[grep('Asia_South_Albers_Equal_Area_Conic', crs_df$label)])
      # Asia_South_Albers_Equal_Area_Conic
      # print(crs_selected)
      pdebug(devug=devug,sep='\n',pre='\n---- uncrs_pts \n','crs_selected') # _____________

      output$ll_coord <- leaflet::renderLeaflet({

        llcrs <- llmap

        if(!is.null(rv$tif_uncrs_extent)){
          # tif_uncrs_extent_p <- SpatialPolygonsDataFrame(
          #   rv$tif_uncrs_extent, data.frame(ID = 1))
          tif_uncrs_extent_p <- sf::st_as_sf(rv$tif_uncrs_extent)
          # tif_uncrs_extent_p@proj4string <- CRS(crs_selected)
          st_crs(tif_uncrs_extent_p) <- crs_selected
          tif_uncrs_extent_p <- sf::st_transform(tif_uncrs_extent_p,
                                                 crs = sf::st_crs("+proj=longlat +ellps=GRS80"))

          llcrs <- llcrs %>% leaflet::addPolygons(data = tif_uncrs_extent_p) %>%
            addLegend("bottomright", colors = c("blue"), labels = c("Raster"), opacity = 1)
        }


        if(!is.null(rv$pts_uncrs_extent)){
          pts_uncrs_extent_p <- sf::st_as_sf(rv$pts_uncrs_extent)
          st_crs(pts_uncrs_extent_p) <- (crs_selected2)
          pts_uncrs_extent_p <- sf::st_transform(pts_uncrs_extent_p,
                                                 crs = sf::st_crs("+proj=longlat +ellps=GRS80"))

          pts_uncrs_p <- rv$pts_uncrs
          #pts_uncrs_p@proj4string <- CRS(crs_selected)
          st_crs(pts_uncrs_p) <- (crs_selected2)
          pts_uncrs_p <- sf::st_transform(pts_uncrs_p, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          #pts_uncrs_p@data[, c('x0', 'y0')] <- st_coordinates(pts_uncrs_p)
          pts_uncrs_p$x0<- st_coordinates(pts_uncrs_p)[, 1]
          pts_uncrs_p$y0<- st_coordinates(pts_uncrs_p)[, 2]

          llcrs <- llcrs %>%  addCircleMarkers(data = pts_uncrs_p,
                                               lng =  ~x0, lat = ~y0, color = 'red') %>%
            addPolygons(data = pts_uncrs_extent_p,color = 'red', fillColor = 'red') %>%
            addLegend("bottomright", colors = c("red"), labels = c("Points"), opacity = 1)
        }

        llcrs
      })

      # crs_db = NULL, # crs_pts_temp = "", # crs_pts = "",
      # crs_tif_temp = "", # crs_tif = "",

      # input$sel_crs # "" # input$in_uncrs_pts # input$in_uncrs_tif
      # input$coo_pts # input$coo_pts
    }
  })

  observeEvent(input$coo_pts, {
    if(!is.null(rv$pts_uncrs_extent)){

    }
  })

  observeEvent(input$coo_tif, {
    if(!is.null(rv$tif_uncrs_extent)){

    }
  })



  ####### SRV PERFORMANCE  ------------------

  output$perftable <- DT::renderDataTable(
    dat <- datatable(per, editable = F,
                     options = list(
                       paging =TRUE # , pageLength =  nrow(rv$data)
                     )
    )
  )

  output$scetable <- DT::renderDataTable(
    dat <- datatable(
      uper , options = list(
        paging =TRUE # , pageLength =  nrow(rv$data)
      )
    )
  )

  observe({

    logx <- input$logx
    logy <- input$logy
    factx <- input$factx
    xaxis <- c('npix', 'spix', 'size')[c('Total pixels', 'Side-pixels', 'Size order') %in% input$xaxis]
    soft <- c('lcc', 'crk', 'mat')[c('Least cost path',
                                     'Kernel density',
                                     'Distance matrix') %in% input$soft]
    # input <- list(xaxis = 'size')

    # factx logx logy xaxis: npix spix size | 'Total pixels', 'Side-pixels', 'Size order'

    # soft <- c('lcc', 'crk', 'mat')[1] #
    # xaxis <- c('npix', 'spix')[1]
    xlab <- input$xaxis# c('Number of total pixels', 'Number of side pixels')[1]
    ylabs <- c(ram = 'max RAM (GB)', cpu = 'CPU time (minutes)')

    # ram <- subset(per, variable == 'ramgb' & func == soft)
    # cpu <- subset(per, variable == 'mins' && func == soft)
    # which(per$variable == 'ramgb' & per$func == soft)

    per$xaxis <- per[, c(xaxis)]
    per$yaxis <- per$value

    # factx logx logy npix spix  size

    if (logx){
      per$xaxis <- log(per$xaxis)
      xlab <- paste0('Log(', xlab, ')')
    }

    if (logy){
      per$yaxis <-log(per$yaxis)
      ylabs <- paste0('Log(', ylabs, ')')
    }

    if (factx & (xaxis != 'size') ){
      per$xaxis <- as.factor(per$xaxis)
    }


    ram <- per[which(per$variable == 'ramgb' & per$func == soft), ]
    cpu <- per[which(per$variable == 'mins' & per$func == soft), ]

    ram <- ram[order(ram$xaxis), ]
    cpu <- cpu[order(cpu$xaxis), ]

    # head(per)
    # tail(per)
    # tail(cpu)

    # https://rpubs.com/rsaidi/676158

    output$hcout1 <- highcharter::renderHighchart({ # cpu
      CPU_chart <<- highchart() %>% hc_exporting(enabled = TRUE) %>%
        hc_add_series(data = subset(cpu, soft == 'COLA'),
                      type = "line", dashStyle = "DashDot",
                      hcaes(x = xaxis,y = yaxis, #color = npts,
                            group = npts)) %>%
        hc_add_series(data = subset(cpu, soft != 'COLA'),
                      type = "line",
                      hcaes(x = xaxis,y = yaxis, #color = npts,
                            group = npts))  %>%
        hc_yAxis(title = list(text = as.character(ylabs['cpu']))) %>%
        hc_xAxis(title = list(text = xlab)) %>%
        hc_title(text = paste0('Performance in CPU time (mins) of ', input$soft) ) %>%
        hc_plotOptions(series = list(marker = list(symbol = "circle"))) %>%
        hc_plotOptions(series = list(marker = list(symbol = "circle"))) %>%

        hc_legend(align = "right", verticalAlign = "top") %>%
        hc_tooltip(shared = F, borderColor = "black",
                   pointFormat = paste0("Software: {point.soft}<br>",
                                        "Scenario: {point.scen}<br>",
                                        "N. pixels: {point.npix}<br>",
                                        "Side pixels: {point.spix}<br>",
                                        "N. points: {point.npts}<br>",
                                        "RAM: {point.value:.2f}<br>") # #pointFormat = paste0("Scenario: {point.scen}<br>")
        )%>% hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))
      CPU_chart
    }) #

    output$hcout2 <- highcharter::renderHighchart({ # ram

      RAM_chart <<- highchart() %>% hc_exporting(enabled = TRUE) %>%
        hc_add_series(data = subset(ram, soft == 'COLA'),
                      type = "line", dashStyle = "DashDot",
                      hcaes(x = xaxis,y = yaxis, #color = npts,
                            group = npts)) %>%
        hc_add_series(data = subset(ram, soft != 'COLA'),
                      type = "line",
                      hcaes(x = xaxis,y = yaxis, #color = npts,
                            group = npts))  %>%
        hc_yAxis(title = list(text = as.character(ylabs['ram']))) %>%
        hc_xAxis(title = list(text = xlab)) %>%
        hc_title(text = paste0('Performance in RAM (GB) of ', input$soft) ) %>%
        hc_plotOptions(series = list(marker = list(symbol = "circle"))) %>%
        hc_legend(align = "right", verticalAlign = "top") %>%
        hc_tooltip(shared = F, borderColor = "black",
                   pointFormat = paste0("Software: {point.soft}<br>",
                                        "Scenario: {point.scen}<br>",
                                        "N. pixels: {point.npix}<br>",
                                        "Side pixels: {point.spix}<br>",
                                        "N. points: {point.npts}<br>",
                                        "RAM: {point.value:.2f}<br>") # #pointFormat = paste0("Scenario: {point.scen}<br>")
        ) %>% hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))

      RAM_chart
    }) #

  })

  ####### SRV CDPOP  ------------------
  # in_cdpop_cd
  # in_cdpop_age
  # in_cdpop_xy

  #   name size type datapath
  # input <- list(in_cdpop_par$datapath = "/data/temp/2023082820202105_file3c454758300e/",
  #               in_cdpop_age$datapath = "/tmp/Rtmp083gNo/9489e1ac4d4f142bb1718552/0.csv",
  #               in_cdpop_cd$datapath = "",
  #               in_cdpop_par$datapath = "")

  if(FALSE){
    #tempFolder <- '/data/temp/2023082821124805_file5762732c645'
    #setwd(tempFolder)

    cdpop_path_invars <- paste0(tempFolder, '/invars.csv')
    cdpop_path_xy <- paste0(tempFolder, '/xy.csv')
    cdpop_path_age <- paste0(tempFolder, '/age.csv')
    cdpop_path_cdmat <- paste0(tempFolder, '/cdmat.csv')

    writeParamCDPOP <- read.csv(file = cdpop_path_invars)
    #invarsorig <- writeParamCDPOP
    xyfile <- read.csv(file = cdpop_path_xy)
    age <- read.csv(file = cdpop_path_age)
    cdmat <- read.csv(file = cdpop_path_cdmat, header = FALSE)
  }

  observeEvent(input$crv, {
    checkEnv()
  })

  observeEvent(input$in_cdpop_par, {
    if(devug){ print('cdpop_par: '); print(input$in_cdpop_par)}
    file <- input$in_cdpop_par
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    invarsorig <- (read.csv(file$datapath, header = input$header))
    rv$orig <- t(invarsorig)
    colnames(rv$orig) <- paste0('Scen', 1:ncol(rv$orig))
    rv$data <- (rv$orig)
    cdpop_path_invars <<- paste0(tempFolder, '/invars.csv')
    #cdpop_path_invars0 <<- paste0(tempFolder, '/invars0.csv')
    write.csv(invarsorig, file = cdpop_path_invars, row.names = FALSE, quote = FALSE)
    #write.csv(invarsorig, file = cdpop_path_invars0, row.names = FALSE)
  })


  observeEvent(input$in_cdpop_xy, {
    if(devug){ print('in_cdpop_xy: '); print(input$in_cdpop_xy)}
    xyfile <- input$in_cdpop_xy
    ext <- tools::file_ext(xyfile$datapath)
    req(xyfile)
    validate(need(ext == "csv", "Please upload a csv file"))
    if(devug){ print(' reading in_cdpop_xy: '); print(xyfile$datapath)}
    xyfile <<- xyfile <- (read.csv(xyfile$datapath, header = input$header))
    print(dim(xyfile))
    #colnames(rv$orig) <- paste0('Scen', 1:ncol(rv$orig))
    #rv$data <- (rv$orig)
    cdpop_path_xy <<- paste0(tempFolder, '/xy.csv')
    file.copy(xyfile$datapath, cdpop_path_xy);
    #try(file.remove(xyfile$datapath) )
    #write.csv(xyfile, file = cdpop_path_xy, row.names = FALSE, quote = FALSE)
    if(!file.exists(cdpop_path_xy)){
      print('Error copy cdpop xy')
      write.csv(xyfile, file = cdpop_path_xy, row.names = FALSE, quote = FALSE)
      if(!file.exists(cdpop_path_xy)){
        print('Error copy cdpop xy 2')
        file.copy(xyfile$datapath, cdpop_path_xy);
        #try(file.remove(xyfile$datapath))
      }
    }

  })

  observeEvent(input$in_cdpop_age, {
    if(devug){ print('in_cdpop_age: '); print(input$in_cdpop_age)}

    agefile <- input$in_cdpop_age
    ext <- tools::file_ext(agefile$datapath)
    req(agefile)
    validate(need(ext == "csv", "Please upload a csv file"))
    age <<- (read.csv(agefile$datapath, header = input$header))
    #colnames(rv$orig) <- paste0('Scen', 1:ncol(rv$orig))
    cdpop_path_age <<- paste0(tempFolder, '/age.csv')
    write.csv(age, file = cdpop_path_age, row.names = FALSE, quote = FALSE)
  })


  #input <- list(in_cdpop_cd$)
  #cdmatfile <- list(datapath = '/tmp/RtmpxR2z6w/4f50c5b3fc19eb4be395c7f1/0.csv')
  observeEvent(input$in_cdpop_cd, {
    if(devug){ print('in_cdpop_cd: '); print(input$in_cdpop_cd)}
    #print(input$in_cdpop_cd)

    cdmatfile <- input$in_cdpop_cd
    ext <- tools::file_ext(cdmatfile$datapath)
    req(cdmatfile)
    validate(need(ext == "csv", "Please upload a csv file"))
    cdmat <<- (read.csv(cdmatfile$datapath, header = FALSE))
    #colnames(rv$orig) <- paste0('Scen', 1:ncol(rv$orig))
    cdpop_path_cdmat <<- paste0(tempFolder, '/cdmat.csv')
    write.table(cdmat, file = cdpop_path_cdmat, row.names = FALSE, col.names = F, sep = ',')
    file.copy(cdmatfile$datapath, cdpop_path_cdmat);
    #try(file.remove(cdmatfile$datapath))
  })

  # output$selectUI<-renderUI({
  #   req(rv$data)
  #   selectInput(inputId='selectcolumn', label='select column', choices = names(rv$data))
  # })

  #splitcolumn
  observeEvent(input$Splitcolumn, {
    rv$data <- splitColumn(rv$data, input$selectcolumn)
  })

  #delterows
  observeEvent(input$deleteRows,{
    if (!is.null(input$table1_rows_selected)) {
      rv$data <- rv$data[-as.numeric(input$table1_rows_selected),]
    }
  })


  ####### CDPOP  ------------------

  # output$table1 <- renderDT({
  #   datatable(rv$data, editable = TRUE)
  # })

  output$table1 <- DT::renderDataTable(
    dat <- datatable(rv$data, editable = TRUE,
                     options = list(
                       paging =TRUE,
                       pageLength =  nrow(rv$data)
                     )))

  observeEvent(input$table1_cell_edit, {
    roww  <- input$table1_cell_edit$row
    clmn <- input$table1_cell_edit$col
    rv$data[roww, clmn] <- input$table1_cell_edit$value
  })

  observeEvent(input$replacevalues, {
    rv$data <- fillvalues(rv$data, input$textbox, input$selectcolumn)
  })
  observeEvent(input$removecolumn, {
    #rv$data <- removecolumn(rv$data,input$selectcolumn)
    rv$data <- removecolumn(rv$data)
  })
  observeEvent(input$addcolumn, {
    rv$data <- addcolumn(rv$data)
  })

  observeEvent(input$Undo, {
    rv$data <- rv$orig
  })

  observeEvent(input$cdpop_run, {
    #xyfilename no requires .csv
    #agefilename requires .csv
    rv$data <- rv$orig
  })

  output$cdpop_box1 <- output$dist_box1 <- shinydashboard::renderValueBox({
    valueBox( "Not ready", '', #"Ready", icon = icon("thumbs-up", lib = "glyphicon"),
              color = "red"
    )
  })

  observeEvent(input$cdpop_check1, {
    writeParamCDPOP <- t(rv$data)[1, ]
    #writeParamCDPOP <- t(writeParamCDPOP)
    writeParamCDPOP$xyfilename <- 'xy'
    writeParamCDPOP$agefilename <- 'age.csv'
    writeParamCDPOP$matecdmat <- writeParamCDPOP$dispcdmat <- 'cdmat'
    write.csv(writeParamCDPOP, file = cdpop_path_invars, row.names = FALSE, quote = FALSE)
    cond <- TRUE
    if( cond){
      output$cdpop_box1 <- shinydashboard::renderValueBox({
        valueBox(
          "YES", "Ready", icon = icon("thumbs-up", lib = "glyphicon"),
          color = "green"
        )
      })
    }
  })

  output$cdpop_box2 <- shinydashboard::renderValueBox({
    valueBox(
      "Not yet", "Ready to run CDPOP", icon = icon("thumbs-up", lib = "glyphicon"),
      color = "red"
    )
  })

  observeEvent(input$cdpop_check2, {
    cdmat2 <- cdmat[, apply(!is.na(cdmat), 2, sum) != 0]
    (cond <- (ncol(cdmat2) == nrow(cdmat2)) &
        ncol(cdmat2) > 2 &
        (nrow(cdmat2) == nrow(xyfile)))

    if( cond){
      output$cdpop_box2 <- shinydashboard::renderValueBox({
        valueBox(
          "YES", "Ready", icon = icon("thumbs-up", lib = "glyphicon"),
          color = "green"
        )
      })

      cdpopRun <- runCDPOP(py, datapath = tempFolder)
      rv$cdpopRun <- cdpopRun
      output$out_cdpop_files<-renderUI({
        req(rv$cdpopRun)
        #cdpop_outfiles <- gsub(pattern = cdpopRun$cdpopPath, replacement = '', x = cdpopRun$newFiles)
        selectInput(inputId='out_cdpop_files', label='select column', choices = basename(cdpopRun$newFiles))
      })
    }
  })

  isolate(observeEvent(input$run_cdpop, {

    # tempFolder = '/tmp/RtmpYiPPnn/colaGBF2024072213145005'; setwd(tempFolder)
    # rv <- list(pts = 'out_simpts_MFU2024072213183205.shp')

    # tempFolder = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305'; setwd(tempFolder)
    # rv <- list(pts = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305/out_simpts_ZFE2024080118415505.shp',
    #            cdm = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305/out_cdmatrix_XUV2024080118421705.csv')
    # pref <- 'SLW'

    newxy <- gsub('.shp', '.csv', x = rv$pts)
    shp2xy(shapefile = rv$pts, outxy = newxy, tempDir = tempFolder)
    rv$ptsxy <- newxy
    # rv$cdm <- 'cdmat.csv'
    (pref <- sessionIDgen(only3 = TRUE))

    output$vout_cdp <- isolate(renderText({

      tStartCDP <- Sys.time()

      cdpop_ans <<- cdpop_py(inputvars = NULL, agevars = NULL,
                            cdmat = rv$cdm, xy = rv$ptsxy,
                            tempFolder = tempFolder, prefix = pref)
      #save(cdpop_out, file = 'cdpop_out.RData'); load('cdpop_out.RData')
      rv$cdpop_ans <<- cdpop_ans
      cdpop_grids <<- grep(pattern = 'grid.+.csv$', x = cdpop_ans$newFiles, value = TRUE)


      tStopCDP <- Sys.time() - tStartCDP
      textElapCDP <- paste(round(as.numeric(tStopCDP), 2), attr(tStopCDP, 'units'))
      mssg2Display <- paste0(' CDPOP simulation finished on ', textElapCDP, '. Generated ', length(cdpop_grids), ' years.')
    }))
      # inputvars = NULL; agevars = NULL; cdmat = rv$cdm; xy = rv$ptsxy; tempFolder = tempFolder; prefix = pref

    output$ll_map_cdp <- leaflet::renderLeaflet({

      # rv <- list()
     # cdpop_out <- '/mnt/c/tempRLinux/colaZYT2024080514114105/GJL__1722885170/batchrun0mcrun0/output.csv'
      cdpop_out <- grep(pattern = 'output.csv$', x = cdpop_ans$newFiles, value = TRUE)
      rv$cdpop_out <- read.csv(cdpop_out)
      rv$cdpFolder <- dirname(dirname(cdpop_out))

      cdpop2Plot <- sapply(rv$cdpop_out[, c('Year', 'Population_Age1.', 'Alleles', 'He', 'Ho')],
                           function(x){
                             #x = rv$cdpop_out$Population_Age1.
                             as.numeric(gsub('\\|.+|\\|', '', x))
                             }
                           )

      cdpop_grids_Num <- as.numeric(gsub('grid|\\.csv', '', basename(cdpop_grids) ) )
      cdpop_grids <- cdpop_grids[order( cdpop_grids_Num )]
      cdpop_grids_Num <- sort(cdpop_grids_Num)

      output$hccdpop1 <- highcharter::renderHighchart({ # cpu
        highchart() %>% hc_exporting(enabled = TRUE) %>%
          hc_add_series(data = cdpop2Plot[, c('Year', 'Population_Age1.')],
                        type = "line", dashStyle = "DashDot", name = 'Population'
                        #, hcaes(x = 'Year', y = 'Population')
          ) %>%
          hc_yAxis(title = list(text = 'Population')) %>%
          hc_xAxis(title = list(text = 'Year')) %>%
          hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))

      })

      output$hccdpop2 <- highcharter::renderHighchart({ # cpu
        highchart() %>% hc_exporting(enabled = TRUE) %>%
          hc_add_series(data = cdpop2Plot[, c('Year', 'Alleles')],
                        type = "line", dashStyle = "DashDot", name = 'Alleles'
                        #, hcaes(x = 'Year', y = 'Population')
          ) %>%
          hc_yAxis(title = list(text = 'Alleles')) %>%
          hc_xAxis(title = list(text = 'Year')) %>%
          hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))
      })

      output$hccdpop3 <- highcharter::renderHighchart({ # cpu
        highchart() %>% hc_exporting(enabled = TRUE) %>%
          hc_add_series(data = cdpop2Plot[, c('Year', 'He')],
                        type = "line", dashStyle = "DashDot", name = 'He'
                        #, hcaes(x = 'Year', y = 'Population')
          ) %>% hc_add_series(data = cdpop2Plot[, c('Year', 'Ho')],
                              type = "line", dashStyle = "DashDot", name = 'Ho'
                              #, hcaes(x = 'Year', y = 'Population')
          ) %>%
          hc_yAxis(title = list(text = 'Heterozygosity')) %>%
          hc_xAxis(title = list(text = 'Year')) %>%
          hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))
      }) #

      updateSelectizeInput(session, inputId = 'cdpop_ans_yy',
                           choices = cdpop_grids_Num,
                           selected = tail(cdpop_grids_Num, 1), server = TRUE)


      # rv$tif <- "/mnt/c/tempRLinux/RtmpNGsi0b/colaADC2024073113022305/out_surface_ZRY2024073113024005.tif"

      densMap <- cdpop_mapdensity(grids = cdpop_grids[1], template = rv$tif,
                                  method = 'average',
                                  bandwidths = 'None', type = 'count', crs = 'None')
      # grids = cdpop_grids[1]; template = rv$tif; method = 'thin_plate_spline'; neighbors = 'all'; crs = 'None'

      struMap <- cdpop_mapstruct(grids = cdpop_grids[1], template = rv$tif,
                                 method = 'thin_plate_spline',
                                 neighbors = 'all', crs = 'None')

      rv$struRA <- struRA <- rast(struMap$newFiles[1])
      rv$struRB <- struRB <- rast(struMap$newFiles[2])
      rv$densR <- densR <- rast(densMap$file)

      rng_strA <- getMxMn(struRA); palA <- leaflet::colorNumeric(palette = "viridis", reverse = TRUE, domain = rng_strA+0.0, na.color = "transparent")
      rng_strB <- getMxMn(struRB); palB <- leaflet::colorNumeric(palette = "viridis", reverse = TRUE, domain = rng_strB+0.0, na.color = "transparent")
      rng_dens <- getMxMn(densR); palC <- leaflet::colorNumeric(palette = "viridis", reverse = TRUE, domain = rng_dens+0.0, na.color = "transparent")

      llcdp <<- leaflet() %>% addTiles() %>%
        addRasterImage(struRA, colors = palA, opacity = .7, group = "Alleles", layerId = 'Alleles') %>%
        addLegend(pal=palA, values=rng_strA, group = 'Alleles', position = 'topleft', title="Alleles") %>%
        addRasterImage(struRB, colors = palB, opacity = .7, group = "Heterozygosity", layerId = 'Heterozygosity') %>%
        addLegend(pal=palB, values=rng_strB, group = 'Heterozygosity', position = 'topleft', title="Heterozygosity") %>%
        addRasterImage(densR, colors= palC, opacity = .7, group = "Density", layerId = 'Density') %>%
        addLegend(pal=palC, values=rng_dens, group = 'Density', position = 'topleft', title="Density") %>%
        leaflet::addLayersControl(
          baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
          overlayGroups = c('Alleles', "Heterozygosity", "Density"),
          options =  leaflet::layersControlOptions(collapsed = FALSE)) %>%
        leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery")
    })

  }))

  # observeEvent(mapcdpop, {
  #   input$cdpop_ans_yy
  #   if(input$cdpop_ans_yy != ''){
  #     rv$cdpop_ans
  #   }
  # })


  observeEvent(input$cdpop_check3, {
    if(input$out_cdpop_files != ''){
      selfile <- cdpopRun$newFiles[basename(cdpopRun$newFiles) %in% input$out_cdpop_files]
      out_cdpop_file_X <- as.data.frame(readLines(
        paste0(tempFolder,'/', selfile)))

      if(devug){
        print(' ---- input$out_cdpop_files: ')
        print(input$out_cdpop_files);
        print(' ---- rv$cdpopRun: ')
        print(rv$cdpopRun)
        print(dim(out_cdpop_file_X))
        print(head(out_cdpop_file_X))
      }

      output$out_cdpop_filestable <- DT::renderDataTable(
        dat <- datatable(out_cdpop_file_X,
                         options = list( paging =TRUE,
                                         pageLength =  nrow(out_cdpop_file_X)
                         )
        ))
    }
  })

  ####### LOAD MAPS  ------------------

  ####### > SURFACE  ------------------

  observeEvent(input$in_sur_tif, {

    output$ll_map_h2r <- leaflet::renderLeaflet({

      # if(is.na(newtifPath)){
      #   rv$log <- paste0(rv$log, '\n -- Error uploading the "Habitat suitability" TIF file')
      #   updateVTEXT(rv$log)
      # } else {
      #   rv$newtifPath <- newtifPath
      #   rv$edi <- newtifPath
      #   rv$tif <- newtifPath
      #   rv$ediready <- TRUE
      #   rv$tifready <- TRUE
      #
      #   rv$tif_sp <- terra::rast(rv$tif)
      #   rv$tif_rng <- rng_newtif <- range(rv$tif_sp[], na.rm = TRUE)
      #   rv$tif_pal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
      #                                domain = rng_newtif+0.0, na.color = "transparent")
      #
      #   rv$edi_sp <- terra::rast(rv$edi)
      #   rv$edi_rng <- rng_newtif <- range(rv$edi_sp[], na.rm = TRUE)
      #   rv$edi_pal <- ediPal <<- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
      #                                          domain = rng_newtif, na.color = "transparent")
      #   rv$log <- paste0(rv$log, '--- DONE'); updateVTEXT(rv$log)
      #   makeLL()

      # (inEdiSessID <<- sessionIDgen())
      # rv$inEdiSessID <<- inEdiSessID
      # tifpath <<- paste0(tempFolder, '/in_edit_', inEdiSessID, '.tif')
      # tifpathfixed <- paste0(tempFolder, '/in_edit_fixed_', inEdiSessID, '.tif')
      #
      #
      # file.copy(input$in_edi_tif$datapath, tifpath);
      # rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems')
      # updateVTEXT(rv$log)
      #
      # newtifPath <- fitRaster2cola(inrasterpath = tifpath, outrasterpath = tifpathfixed)
      # newtifPath <- ifelse(is.na(newtifPath), yes = tifpath, no = newtifPath)


      (inSurSessID <<- sessionIDgen())
      rv$inSurSessID <<- inSurSessID
      rv$orighs <- tifpath <- paste0(tempFolder, '/in_surface_', rv$inSurSessID, '.tif')
      tifpathfixed <- paste0(tempFolder, '/in_surface_fixed_', rv$inSurSessID, '.tif')


      file.copy(input$in_sur_tif$datapath, tifpath);
      if(file.exists(tifpath)){file.remove(input$in_sur_tif$datapath)}

      # try(file.remove(input$in_sur_tif$datapath))

      pdebug(devug=F,sep='\n',pre='- A',"tempFolder","inSurSessID",
             "rv$inSurSessID", "input$in_sur_tif$datapath",
             "rv$hsready", "tifPath", "tifpathfixed",
             "paste0(tempFolder, '/in_surface_', inSurSessID, '.tif')"
      )

      rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems')
      updateVTEXT(rv$log)

      newtifPath <- fitRaster2cola(inrasterpath = tifpath, outrasterpath = tifpathfixed)
      newtifPath <<- ifelse(is.na(newtifPath), yes = tifpath, no = newtifPath)


      if(is.na(newtifPath)){
        rv$log <- paste0(rv$log, '\n -- Error uploading the "Habitat suitability" TIF file')
        updateVTEXT(rv$log)
        makeLL()
      } else {
        params_txt <- updateParamsTEXT(params_txt = params_txt, hs = TRUE)
        rv$newtifPath <- newtifPath
        rv$hs <- newtifPath
        rv$hsready <- TRUE

        pdebug(devug=devug,sep='\n',pre='-',"tempFolder","inSurSessID",
               "rv$inSurSessID", "newtifPath")

        rv$log <- paste0(rv$log, '--- DONE')
        updateVTEXT(rv$log)

        # newtifPath <- "/data/temp/GA2023090812182205file1266634e12b//in_surface_MN2023090812183705file12666abd01f1.tif"
        #pdebug(devug=devug,sep='\n',pre='-',"tifpath", "newtifPath", "rv$newtifPath", "rv$hs", "rv$hsready")


        rv$hs_sp <- terra::rast(rv$hs)
        #rng_newtif <- c(newtif@data@min, newtif@data@max)
        rv$hs_rng <- rng_newtif <- range(rv$hs_sp[], na.rm = TRUE)

        updateTextInput(session, inputId = "in_sur_3", value = rv$hs_rng[1])
        updateTextInput(session, inputId = "in_sur_4", value = rv$hs_rng[2])


        rv$hs_pal <- hsPal <<- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                                     domain = rng_newtif, na.color = "transparent")

        makeLL()
        }
    })
  })

  ## Run h23
  observeEvent(input$h2r, {

    if(rv$hsready){
      # rv <- list(newtifPath = '/data/temp//E-2023082911285005_file3112135d2b4c//in_surface_V-2023082911285705_file3112303ea820.tif',
      #            inSurSessID = 'V-2023082911285705_file3112303ea820')
      # input <- list(in_sur_3 = 0, in_sur_4 =100, in_sur_5 = 100, in_sur_6 = 1, in_sur_7 = -9999)
      output$ll_map_h2r <- leaflet::renderLeaflet({

        (inSurSessID <- sessionIDgen())
        #(inSurSessID2 <<- sessionIDgen())
        outs2r <- paste0(tempFolder, '/out_surface_', inSurSessID, '.tif')
        #pdebug(devug=devug,sep='\n',pre='---H2S\n',"inSurSessID", 'inSurSessID2', 'outs2r') # = = = = = = =  = = =  = = =  = = =  = = =

        rv$log <- paste0(rv$log,  # _______
                         '\nCreating resistance surface');updateVTEXT(rv$log) # _______

        in_sur_7 <- ifelse(input$in_sur_7 == '', yes = -9999, no = input$in_sur_7)



        hs2rs_file <- s2res_py(py = py,
                               intif = rv$hs,
                               outtif = outs2r,
                               param3 = as.numeric(input$in_sur_3),
                               param4 = as.numeric(input$in_sur_4),
                               param5 = as.numeric(input$in_sur_5),
                               param6 = as.numeric(input$in_sur_6),
                               param7 = in_sur_7,
                               param8 = 'None')

        if(!is.na(hs2rs_file$file)){

          params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)

          rv$log <- paste0(rv$log,  # _______
                           ' ... DONE');updateVTEXT(rv$log) #
          rv$tifready <- TRUE
          rv$tif <- hs2rs_file$file
          rv$tiforig <- hs2rs_file$file

          rv$tif_sp <- hs2rs_tif <- terra::rast(hs2rs_file$file)
          rv$tif_rng <- rng_rstif <- range(hs2rs_tif[], na.rm = TRUE)
          #rv$tif_rng <- rng_rstif <- minmax(rv$tif_sp)[1:2]
          rv$tif_pal <- rsPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                        domain = rng_rstif, na.color = "transparent")


          # rv$llmap rv$hsready rv$tifready rv$ptsready #  rv$llmap
          #rv$llmap <<- rv$llmap %>%
          #leafsurface <<- leaflet::leaflet() %>%leaflet::addTiles() %>%

          # pdebug(devug=devug,sep='\n',pre='---H2S\n'," hs2rs_tif[]") # = = = = = = =  = = =  = = =  = = =  = = =
          makeLL( )

        } else {
          rv$log <- paste0(rv$log, '\n -- Error creating the "Surface resistance" TIF file')
          updateVTEXT(rv$log)
        }


      })
    }
  })



  observeEvent(input$h2rsample, {

    # rv <- list(newtifPath = '/data/temp//E-2023082911285005_file3112135d2b4c//in_surface_V-2023082911285705_file3112303ea820.tif',
    #            inSurSessID = 'V-2023082911285705_file3112303ea820')
    # input <- list(in_sur_3 = 0, in_sur_4 =100, in_sur_5 = 100, in_sur_6 = 1, in_sur_7 = -9999)
    # (hs2rs_samp_file <- system.file(package = 'cola', 'sampledata/sampleTif.tif'))
    # rv <<- list(hs = hs2rs_samp_file)
    output$ll_map_h2r <- leaflet::renderLeaflet({

      (inSurSessID <<- sessionIDgen())
      rv$inSurSessID <- inSurSessID

      params_txt <- updateParamsTEXT(params_txt = params_txt, hs = TRUE)

      rv$hsready <- TRUE
      rv$hs <- hs2rs_samp_file

      rv$hs_sp <- terra::rast(hs2rs_samp_file)
      #rv$hs_rng <- rng_rstif <- range(hs2rs_tif[], na.rm = TRUE)
      rv$hs_rng <- rng_rstif <- getMxMn(rv$hs_sp)[1:2]

      updateTextInput(session, inputId = "in_sur_3", value = rv$hs_rng[1])
      updateTextInput(session, inputId = "in_sur_4", value = rv$hs_rng[2])

      rv$hs_pal <- rsPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                   domain = rng_rstif, na.color = "transparent")
      # print('XXXX')
      # print(hs2rs_samp_file)
      # save(rv, file = 'rv.RData')
      # rv$llmap rv$hsready rv$tifready rv$ptsready #  rv$llmap
      # rv$llmap <<- rv$llmap %>%
      # leafsurface <<- leaflet::leaflet() %>%leaflet::addTiles() %>%

      # pdebug(devug=devug,sep='\n',pre='---H2S\n'," hs2rs_tif[]") # = = = = = = =  = = =  = = =  = = =  = = =
      makeLL( )

    })
  })



  ####### > EDIT  ------------------
  ####### > Draw notes ------------------
  # observeEvent(input$ll_map_edi_draw_new_feature, {
  # https://rdrr.io/cran/leaflet.extras/src/inst/examples/shiny/draw-events/app.R
  #   polDraw <- input$ll_map_edi_draw_new_feature # LEAFLETWIDGET_draw_new_feature
  #   save(polDraw, file = '/data/temp/draw.RData')
  # })
  # observeEvent(input$ll_map_edi_draw_start, {
  #   print("Start of drawing")
  #   print(input$ll_map_edi_draw_start)
  # })
  # # Stop of Drawing
  # observeEvent(input$ll_map_edi_draw_stop, {
  #   print("Stopped drawing")
  #   print(input$ll_map_edi_draw_stop)
  # })
  # # New Feature
  # observeEvent(input$ll_map_edi_draw_new_feature, {
  #   print("New Feature")
  #   print(input$ll_map_edi_draw_new_feature)
  # })
  # # Edited Features
  # observeEvent(input$ll_map_edi_draw_edited_features, {
  #   print("Edited Features")
  #   print(input$ll_map_edi_draw_edited_features)
  # })
  # # Deleted features
  # observeEvent(input$ll_map_edi_draw_deleted_features, {
  #   print("Deleted Features")
  #   print(input$ll_map_edi_draw_deleted_features)
  # })
  # # We also listen for draw_all_features which is called anytime
  # # features are created/edited/deleted from the map
  # observeEvent(input$ll_map_edi_draw_all_features, {
  #   print("VVVV----- All Features ------VVVV")
  #   print(input$ll_map_edi_draw_all_features)
  #   print("^^^^^^----- All Features ------^^^^")
  # })


  ## Load tif
  observeEvent(input$in_edi_tif, {

    #try(file.remove(c(tifpath, newtifPath)))
    invisible(suppressWarnings(tryCatch(file.remove(c(tifpath, newtifPath, tifpathfixed)),
                                        error = function(e) NULL)))

    output$ll_map_edi <- leaflet::renderLeaflet({
      # tempFolder <- '/data/temp//T2023082911164705_file3112795957d7/'
      #tifpath <- '/data/temp//T2023082911164705_file3112795957d7//in_surface_C2023082911165605_file31126374e76.tif'
      #inSurSessID <- 'C2023082911165605_file31126374e76'
      (inEdiSessID <<- sessionIDgen())
      rv$inEdiSessID <<- inEdiSessID
      rv$tiforig <- tifpath <<- paste0(tempFolder, '/in_edit_', inEdiSessID, '.tif')
      tifpathfixed <- paste0(tempFolder, '/in_edit_fixed_', inEdiSessID, '.tif')

      file.copy(input$in_edi_tif$datapath, tifpath);
      if (file.exists(tifpath)){ file.remove(input$in_edi_tif$datapath)}
      #try(file.remove(input$in_sur_tif$datapath))
      #if(devug){ print(' ----- input$in_sur_tif'); print(input$in_sur_tif); print(tifpath); file.exists(tifpath)}

      rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems')
      updateVTEXT(rv$log)

      newtifPath <- fitRaster2cola(inrasterpath = tifpath, outrasterpath = tifpathfixed)
      newtifPath <- ifelse(is.na(newtifPath), yes = tifpath, no = newtifPath)

      if(is.na(newtifPath)){
        rv$log <- paste0(rv$log, '\n -- Error uploading the "Habitat suitability" TIF file')
        updateVTEXT(rv$log)
      } else {
        params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)

        rv$newtifPath <- newtifPath
        rv$edi <- newtifPath
        rv$tif <- newtifPath
        rv$ediready <- TRUE
        rv$tifready <- TRUE


        rv$tif_sp <- terra::rast(rv$tif)
        rv$tif_rng <- rng_newtif <- range(rv$tif_sp[], na.rm = TRUE)
        #rv$tif_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)


        rv$tif_pal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                             domain = rng_newtif+0.0, na.color = "transparent")

        rv$edi_sp <- terra::rast(rv$edi)
        rv$edi_rng <- rng_newtif <- range(rv$edi_sp[], na.rm = TRUE)
        #rv$tif_rng <- rng_newtif <- range(minmax(rv$edi_sp)[1:2], na.rm = TRUE)

        rv$edi_pal <- ediPal <<- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                                       domain = rng_newtif+0.0, na.color = "transparent")

        rv$log <- paste0(rv$log, '--- DONE'); updateVTEXT(rv$log)

        # newtifPath <- "/data/temp/GA2023090812182205file1266634e12b//in_surface_MN2023090812183705file12666abd01f1.tif"
        #pdebug(devug=devug,sep='\n',pre='-',"tifpath", "newtifPath", "rv$newtifPath", "rv$hs", "rv$hsready")

        makeLL( )
      }
    })
  })

  ## Load shp
  observeEvent(input$in_edi_shp, {

    # invisible(suppressWarnings(tryCatch(file.remove(c(in_lcc_shp, newin_lcc_shp)), error = function(e) NULL)))
    pdebug(devug=devug,sep='\n',pre='\n---- LCC - SHP\n','rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________

    rv$log <- paste0(rv$log, '\nLoading scenario shapefile');updateVTEXT(rv$log) # _______

    ## Create session IF if started from this tab
    if(is.null(rv$inEdiSessID)){
      (inEdiSessID <<- sessionIDgen())
      rv$inEdiSessID <- inEdiSessID
    }

    pdebug(devug=devug,sep='\n',pre='--', 'inEdiSessID', 'rv$inEdiSessID') # _____________

    inFiles <- input$in_edi_shp #

    inFiles$newFile <- paste0(tempFolder, '/', basename(inFiles$name))
    pdebug(devug=devug,sep='\n',pre='\n--','print(inFiles)', 'inFiles$newFile', 'tempFolder') # _____________

    file.copy(inFiles$datapath, inFiles$newFile);
    # try(file.remove(inFiles$datapath))
    #if(devug){save(inFiles, file = paste0(tempFolder, '/shpfiles.RData'))}

    inShp <<- loadShp(inFiles, tempFolder, rv$inEdiSessID)
    print( ' =================== print( inShp$shp ')
    print( inShp$shp )

    save(inShp, inFiles, file = paste0(tempFolder, '/debug_edi_shp.RData'))

    if ( any(class(inShp$shp) %in% 'sf') ){
      # if(class(inShp$shp) == 'SpatialPointsDataFrame'){

      rv$sceready <- TRUE
      rv$sce <- inShp$layer
      rv$shp <- inShp$shp
      rv$sce_sp <- inShp$shp
      rv$log <- paste0(rv$log, ' -- Scenario shapefile loaded');
      updateVTEXT(rv$log) # _______

      pdebug(devug=devug,sep='\n',pre='---- LOAD SCE SHP\n','rv$sceready', 'rv$sce', 'rv$inEdiSessID') # _____________

      sce <- st_transform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
      sce$sortID <- 1:nrow(sce)
      rv$edi_sp <- sce

      print( ' =================== print(rv$edi_sp) ')
      print(rv$edi_sp)

      bounds <- sce %>% st_bbox() %>% as.character()
      proxy <- leafletProxy("ll_map_edi")
      #proxy <- leaflet() %>% addTiles()
      #print (getwd()); save(sce, file = 'sce_debug_edi_shp.RData'); print (sce)
      proxy %>%  #remove(layerId = 'SurfaceResistance')  %>%
        removeShape('layerId') %>%
        leaflet::addPolygons(data = sce, color = 'darkblue', fillColor = 'darkblue', group = "layerId") %>%
        fitBounds(bounds[1], bounds[2], bounds[3], bounds[4])

      # output$ll_map_edi <- leaflet::renderLeaflet({
      #   makeLL( )
      # })
    }
  })


  ## Run Edi sum ---
  isolate(observeEvent(input$edi, {

    #polDraw <- input$ll_map_edi_draw_new_feature # LEAFLETWIDGET_draw_new_feature
    polDraw <- input$ll_map_edi_draw_all_features # LEAFLETWIDGET_draw_new_feature
    save(polDraw, file = paste0(tempFolder, '/draw.RData'))

    if(input$in_edi_val != 0 & input$in_edi_val != "0" & rv$tifready & (!is.null(polDraw) | isTRUE(rv$sceready)) ){
      # rv <- list(tif = '/data/temp/XS2023100319220605file859285936e77a/in_edit_TG2023100319221605file8592817dc90f8.tif')
      # rv$tif_sp <- terra::rast(rv$tif)
      # input <- list(in_sur_3 = 0, in_sur_4 =100, in_sur_5 = 100, in_sur_6 = 1, in_sur_7 = -9999)

      # print(input$in_edi_val)
      # print(polDraw)
      # print(num2Burn)
      # (load('/data/tempR/draw4pol.RData'))
      # (load('/home/shiny/cola/inst/app/sce_debug_edi_shp.RData')); print (sce)


      ## Not required as string or numeric is passed to the
      # if (in_edi_val %in% colnames(sce)){
      #     num2Burn <- in_edi_val
      #   }

      output$ll_map_edi <- leaflet::renderLeaflet({
        (inEdiSessID2 <<- sessionIDgen())
        rv$inEdiSessID2 <<- inEdiSessID2
        editRastpath <<- paste0(tempFolder, '/in_editrast_', inEdiSessID2, '.tif')
        editShppath <<- paste0(tempFolder, '/in_editrast_', inEdiSessID2, '.shp')

        rv$log <- paste0(rv$log,  # _______
                         '\n -- Creating scenario');updateVTEXT(rv$log) #

        num2Burn <<- input$in_edi_val
        lineBuffW <<- input$in_edi_wid

        # num2Burn <- 'try'
        # in_edi_val <- input$in_edi_val
        # #in_edi_val <- 'ID'
        pdebug(devug=devug,sep='\n',pre='-',"num2Burn", "input$in_edi_val", "rv$tif")

        ## Define if use shape or draw
        if( isTRUE(rv$sceready) ){

          ## Uses shapefile
          pol2Rastx <- rv$sce_sp # polPath
          # pol2use <- sce # polPath

          ## Assign column
          if(num2Burn %in% colnames(pol2Rastx)){ # num2Burn <- 'ID'
            pol2Rastx$val2burn <- as.data.frame(pol2Rastx)[, num2Burn] # polPath
          } else {
            pol2Rastx$val2burn <- as.numeric(num2Burn)
          }

          polPath <- rv$sce
          print(paste0(' ----- SAVE: ', tempFolder, '/polPath_pol2Rastx.RData'))
          save(pol2Rastx, polPath, num2Burn, file = paste0(tempFolder, '/polPath_pol2Rastx.RData'))
          polPath <- paste0(tools::file_path_sans_ext(basename(polPath)), '_scepol.shp')

        } else {

          ## Uses draw polygon
          rt <- terra::rast(rv$tif)
          rastRes <- res(rt)
          (rastCRS <- st_crs(rt))

          pol2Rast <- draws2Features(polDraw, rastCRS = rastCRS, distLineBuf = min(rastRes) * lineBuffW )
          # pol2Rastx <- st_sf(data.frame(a = 1:length(pol2Rast), pol2Rast))
          pol2Rastx <- st_as_sf( pol2Rast)
          pol2Rastx$val2burn <- as.numeric(num2Burn)

          polPath <- gsub('.tif$', '_scepol.shp', rv$tif)
        }

        print(paste0('polPath:', polPath ))

        ## Write polygon

        sf::st_write( obj = pol2Rastx,
                      dsn = dirname(polPath),
                      layer = tools::file_path_sans_ext(basename(polPath)),
                      driver = 'ESRI Shapefile',
                      append = FALSE,
                      overwrite_layer = TRUE)


        ## Burn the value of the polygon into the rast
        burned <<- burnShp(polPath = polPath,
                           burnval = 'val2burn',
                           rastPath = rv$tif,
                           lineBuffW = as.numeric(input$in_edi_wid),
                           att = input$in_edi_che,
                           rastCRS = NA)

        pdebug(devug=devug,sep='\n',pre='-',"burned")

        if(!is.na(burned)){
          rv$log <- paste0(rv$log,  # _______
                           ' ... DONE');updateVTEXT(rv$log) #

          rv$tifready <- TRUE
          rv$tif <- burned
          rv$tiforig <- burned
          rv$tif_sp <- terra::rast(rv$tif)
          rv$tif_rng <- getMxMn(rv$tif_sp)

          # rv$tif2s <- resampIfNeeded(burned)
          # rv$tif0 <- burned
          # rv$tif_sp <- terra::rast(burned)
          # rv$tif2s_sp <- terra::rast(rv$tif2s)

          rv$tif_rng <- rng_rstif <- range(rv$tif_sp[], na.rm = TRUE)
          #rv$tif_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)

          rv$tif_pal <- rsPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                        domain = rng_rstif, na.color = "transparent")

          makeLL( )

        } else {
          rv$log <- paste0(rv$log, '\n -- Error creating the "Surface resistance" TIF file')
          updateVTEXT(rv$log)
        }
      })
    }
  })
  )

  ## Run Edi replace ---
  isolate(observeEvent(input$rpl, {
    polDraw <- input$ll_map_edi_draw_all_features # LEAFLETWIDGET_draw_new_feature

    if(input$in_edi_val != 0 & input$in_edi_val != "0" & rv$tifready & (!is.null(polDraw) | isTRUE(rv$sceready)) ){

      output$ll_map_edi <- leaflet::renderLeaflet({

        (inEdiSessID2 <<- sessionIDgen())
        rv$inEdiSessID2 <<- inEdiSessID2
        editRastpath <<- paste0(tempFolder, '/in_editrast_', inEdiSessID2, '.tif')
        editShppath <<- paste0(tempFolder, '/in_editrast_', inEdiSessID2, '.shp')

        rv$log <- paste0(rv$log,  # _______
                         '\n -- Creating scenario');updateVTEXT(rv$log) #

        num2Burn <<- input$in_edi_val
        lineBuffW <<- input$in_edi_wid

        # num2Burn <- 'try'
        # in_edi_val <- input$in_edi_val
        # #in_edi_val <- 'ID'
        pdebug(devug=devug,sep='\n',pre='-',"num2Burn", "input$in_edi_val", "rv$tif")

        ## Define if use shape or draw
        if( isTRUE(rv$sceready) ){

          ## Uses shapefile
          pol2Rastx <- rv$sce_sp # polPath
          # pol2use <- sce # polPath

          ## Assign column
          if(num2Burn %in% colnames(pol2Rastx)){ # num2Burn <- 'ID'
            pol2Rastx$val2burn <- as.data.frame(pol2Rastx)[, num2Burn] # polPath
          } else {
            pol2Rastx$val2burn <- num2Burn
          }

          polPath <- rv$sce
          pdebug(devug=devug,sep='\n',pre='-', "'Sce shape'", "polPath", polPath)



        } else {

          ## Uses draw polygon
          rt <- terra::rast(rv$tif)
          rastRes <- res(rt)
          (rastCRS <- st_crs(rt))

          pol2Rast <- draws2Features(polDraw, rastCRS = rastCRS, distLineBuf = min(rastRes) * lineBuffW )
          # pol2Rastx <- st_sf(data.frame(a = 1:length(pol2Rast), pol2Rast))
          pol2Rastx <- st_as_sf( pol2Rast)
          pol2Rastx$val2burn <- num2Burn

          polPath <- gsub('.tif$', '_scepol.shp', rv$tif)
          pdebug(devug=devug,sep='\n',pre='-', "'Sce draw'", "polPath", polPath)

        }



        ## Write polygon
        sf::st_write( obj = pol2Rastx,
                      dsn = dirname(polPath),
                      layer = tools::file_path_sans_ext(basename(polPath)),
                      driver = 'ESRI Shapefile',
                      append=FALSE,
                      overwrite_layer = TRUE)


        ## Burn the value of the polygon into the rast
        burned <<- replaceRastShp(polPath = polPath,
                                  lineBuffW = as.numeric(input$in_edi_wid),
                                  att = input$in_edi_che,
                                  burnval = 'val2burn',
                                  rastPath = rv$tif,
                                  rastCRS = NA)

        pdebug(devug=devug,sep='\n',pre='-',"burned")

        if(!is.na(burned)){
          rv$log <- paste0(rv$log,  # _______
                           ' ... DONE');updateVTEXT(rv$log) #

          rv$tifready <- TRUE
          rv$tif <- burned
          rv$tiforig <- burned
          rv$tif_sp <- terra::rast(burned)

          rv$tif_rng <- rng_rstif <- range(rv$tif_sp[], na.rm = TRUE)
          #rv$tif_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)

          rv$tif_pal <- rsPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                        domain = rng_rstif, na.color = "transparent")

          makeLL( )

        } else {
          rv$log <- paste0(rv$log, '\n -- Error creating the "Surface resistance" TIF file')
          updateVTEXT(rv$log)
        }
      })
    }
  })
  )


  ####### > POINTS  ------------------

  ## Load Res tif in pts tab
  observeEvent(input$in_points_tif, {
    # invisible(suppressWarnings(tryCatch(file.remove(c(rv$tifpathpts, rv$tifpathptsfix)),
    #                                     error = function(e) NULL)))

    if(is.null(rv$inSurSessID)){
      #pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
      pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
    }

    if(is.null(rv$inPointsSessID)){
      (inPointsSessID <- sessionIDgen())
      rv$inPointsSessID <- inPointsSessID
    }

    rv$tiforig <- paste0(tempFolder, '/in_points_', rv$inPointsSessID, '.tif')
    file.copy(input$in_points_tif$datapath, rv$tiforig);
    # try(file.remove(input$in_points_tif$datapath))
    # pdebug(devug=devug,sep='\n',pre='---H2S\n'," hs2rs_tif[]") # = = = = = = =  = = =  = = =  = = =  = = =

    rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______

    rv$tifpathptsfix <- paste0(tempFolder, '/in_points_fixed_', rv$inPointsSessID, '.tif')
    pdebug(devug=devug,sep='\n',pre='-','rv$tiforig')
    newtifPath_pts <- fitRaster2cola(inrasterpath = rv$tiforig, outrasterpath =  rv$tifpathptsfix)
    newtifPath_pts <- ifelse(is.na(newtifPath_pts), yes = rv$tiforig, no = newtifPath_pts)
    rv$newtifpathpts <- newtifPath_pts

    if (file.exists(rv$newtifpathpts)){
      params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)

      rv$log <- paste0(rv$log, ' --- DONE');updateVTEXT(rv$log) # _______
      rv$tifready <- TRUE
      rv$tif <- newtifPath_pts

      output$ll_map_points <- leaflet::renderLeaflet({

        rv$tif_sp <- newtif_pts <- terra::rast(rv$newtifpathpts)
        #rv$tif_rng <- rng_newtif_pts <- range(newtif_pts[], na.rm = TRUE)
        rv$tif_rng <- rng_newtif_pts <- getMxMn(rv$tif_sp)


        updateTextInput(session, inputId = "in_points_3", value = rv$tif_rng[1])
        updateTextInput(session, inputId = "in_points_4", value = rv$tif_rng[2])

        rv$tif_pal <- ptsPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                       domain = rng_newtif_pts, na.color = "transparent")

        makeLL( )
        #
        # llmap <<- rv$llmap %>% removeImage(layerId = 'SurfaceResistance') %>%
        #   removeControl('legendSurface') %>%
        #   addRasterImage(newtif_pts, colors = ptsPal, opacity = .7,
        #                  group = "Surface resistance", layerId = 'SurfaceResistance') %>%
        #   addLegend(pal =  ptsPal, values = newtif_pts[], layerId = 'legendSurface',
        #             position = 'topleft',
        #             title= "Resistance"#, opacity = .3
        #             #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
        #   ) %>%  leaflet::addLayersControl(
        #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
        #     overlayGroups = c("Habitat suitability", "Surface resistance"),
        #     options =  leaflet::layersControlOptions(collapsed = FALSE)
        #   ) %>% clearBounds() %>%  leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" )
        #
        #
        # rv$llmap <<- llmap
        # updateLL(llmap)
        # # leafsurface
        # #llmap
        # rv$llmap
      })
    }
  })

  ## Load HS in pts point
  observeEvent(input$in_points_hs, {
    # invisible(suppressWarnings(tryCatch(file.remove(c(rv$tifpathpts, rv$tifpathptsfix)),
    #                                     error = function(e) NULL)))

    if(is.null(rv$inSurSessID)){
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
      pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
    }

    if(is.null(rv$inPointsSessID)){
      (inPointsSessID <- sessionIDgen())
      rv$inPointsSessID <- inPointsSessID
    }

    rv$tifpathpts <- paste0(tempFolder, '/in_pointshs_', rv$inPointsSessID, '.tif')
    file.copy(input$in_points_hs$datapath, rv$tifpathpts);
    rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______

    rv$tifpathptsfix <- paste0(tempFolder, '/in_pointshs_fixed_', rv$inPointsSessID, '.tif')
    newtifPath_pts <- fitRaster2cola(inrasterpath = rv$tifpathpts, outrasterpath =  rv$tifpathptsfix)
    newtifPath_pts <<- ifelse(is.na(newtifPath_pts), yes = rv$tifpathpts, no = newtifPath_pts)

    pdebug(devug=devug,sep='\n',pre='---PTS\n',"newtifPath_pts") # = = = = = = =  = = =  = = =  = = =  = = =


    if (file.exists(newtifPath_pts)){
      params_txt <- updateParamsTEXT(params_txt = params_txt, hs = TRUE)
      rv$log <- paste0(rv$log, ' --- DONE');updateVTEXT(rv$log) # _______
      rv$hsready <- TRUE
      rv$hs <- newtifPath_pts
      # pdebug(devug=devug,sep='\n',pre='---PTS\n',"rv$hs") # = = = = = = =  = = =  = = =  = = =  = = =

      output$ll_map_points <- leaflet::renderLeaflet({

        # rv <- list(hs = '/data/temp/TO2024011617101505filebfa367c3d30//in_pointshs_QP2024011617103005filebfa42265b6d.tif')
        rv$hs_sp <- newtif_pts <- terra::rast(rv$hs)
        #rv$hs_rng <- rng_newtif_pts <- r ange(newtif_pts[], na.rm = TRUE)
        rv$hs_rng <- rng_newtif_pts <- getMxMn(rv$hs)[1:2]
        pdebug(devug=devug,sep='\n',pre='-','rv$hs_rng')


        updateTextInput(session, inputId = "in_points_3", value = rv$hs_rng[1])
        updateTextInput(session, inputId = "in_points_4", value = rv$hs_rng[2])

        rv$hs_pal <- ptsPal <<- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                                      domain = rng_newtif_pts, na.color = "transparent")
        makeLL( )
      })
    }
  })


  observeEvent(input$points_py, {
    if(! (rv$tifready | rv$hsready)){
      rv$log <- paste0(rv$log, ' \n Creating points -- No raster yet!');updateVTEXT(rv$log) # _______
    } else {
      rv$log <- paste0(rv$log, ' \nCreating points');updateVTEXT(rv$log) # _______

      if(is.null(rv$inSurSessID)){
        pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
        (inSurSessID <- sessionIDgen())
        rv$inSurSessID <- inSurSessID
        pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
      }

      if(is.null(rv$inPointsSessID)){
        (inPointsSessID <- sessionIDgen())
        rv$inPointsSessID <- inPointsSessID
      }

      out_pts <- paste0(tempFolder, '/out_simpts_', rv$inSurSessID, '.shp')

      in_points_ly <<- input$in_points_ly

      pdebug(devug=devug,sep='\n',pre='---PTS\n',
             "inPts", 'rv$in_points_ly','in_points_ly',
             'rv$hs', 'rv$tif') # = = = = = = =  = = =  = = =  = = =  = = = http://18.190.026.82:8787/p/f2dab63b/#shiny-tab-tab_surface

      if(in_points_ly == 'SurfaceResistance'){
        inPts <<- rv$hs
        pdebug(devug=devug,sep='\n',pre='---PTS\n',
               "'SR'", 'rv$in_points_ly','in_points_ly',
               'rv$hs', 'rv$tif')

        points_file <- points_py(py = py,
                                 intif = as.character(rv$tif),
                                 outshp = out_pts,
                                 param3 = as.numeric(input$in_points_3),
                                 param4 = as.numeric(input$in_points_4),
                                 param5 = as.numeric(input$in_points_5))
      }

      if(in_points_ly == 'HabitatSuitability'){
        inPts <<- rv$tif
        pdebug(devug=devug,sep='\n',pre='---PTS\n',
               "'HS'", 'rv$in_points_ly','in_points_ly',
               'rv$hs', 'rv$tif')
        points_file <- points_py(py = py,
                                 intif = as.character(rv$hs),
                                 outshp = out_pts,
                                 param3 = as.numeric(input$in_points_3),
                                 param4 = as.numeric(input$in_points_4),
                                 param5 = as.numeric(input$in_points_5))
      }

      # inPts <<- switch (in_points_ly,
      #                  SurfaceResistance = rv$hs,
      #                  HabitatSuitability = rv$tif)
      # print(points_file)

      # rv$log <- paste0(rv$log, ' \nCreating points');updateVTEXT(rv$log) # _______

      if (!file.exists(points_file$file)){
        rv$log <- paste0(rv$log, '  --- Error creating points');updateVTEXT(rv$log) # _______
      } else {
        params_txt <- updateParamsTEXT(params_txt = params_txt, pts = TRUE)

        rv$pts <- points_file$file
        rv$ptsready <- TRUE

        output$ll_map_points <- leaflet::renderLeaflet({

          #points_file <- "/data/temp/L2023090100204905file18e703e3d6298/out_simpts_J2023090100210305file18e7061e66c55.shp"
          points_shpO <- sf::read_sf(points_file$file)
          points_shp <- sf::st_transform(points_shpO, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          points_shp$sortID <- 1:nrow(points_shp)
          #points_shp@data[, c('lng', 'lat')] <- points_shp@coords
          rv$pts_sp <- points_shp

          rv$log <- paste0(rv$log, '  --- DONE');updateVTEXT(rv$log) # _______

          #temLL <- rv$llmap
          #save(temLL, file = '/data/tempR/ll.RData')
          #load('/data/tempR/ll.RData') # rv <- list(llmap = temLL); llmap = temLL

          makeLL( )

        })
      }
    }
  })

  ## update values
  observeEvent(input$in_points_ly, {
    in_points_ly <<- input$in_points_ly

    if(in_points_ly == 'SurfaceResistance'){
      updateTextInput(session, inputId = "in_points_3", value = rv$tif_rng[1])
      updateTextInput(session, inputId = "in_points_4", value = rv$tif_rng[2])
    }

    if(in_points_ly == 'HabitatSuitability'){
      updateTextInput(session, inputId = "in_points_3", value = rv$hs_rng[1])
      updateTextInput(session, inputId = "in_points_4", value = rv$hs_rng[2])
    }
  })


  ####### > DISTANCE  ------------------

  ##> vout_dist; ll_map_dist; dist_py; in_distance_3,
  ##> in_distance_shp in_dist_tif, inDistSessID distmap newtifPath_dist newshpPath_dist
  # distmap = NULL, distrast = NULL, distshp = NULL,

  observeEvent(input$in_dist_tif, {
    #invisible(suppressWarnings(tryCatch(file.remove(c(rv$tifpathdist, rv$newtifPath_dist)), error = function(e) NULL)))

    if(is.null(rv$inSurSessID)){
      #pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
      #pdebug(devug=devug,sep='\n',pre='-','rv$inSurSessID')
    }

    if(is.null(rv$inDistSessID)){
      (inDistSessID <- sessionIDgen())
      rv$inDistSessID <- inDistSessID
    }


    rv$tiforig <- paste0(tempFolder, '/in_dist_', rv$inDistSessID, '.tif')
    file.copy(input$in_dist_tif$datapath, rv$tiforig);
    # try(file.remove(input$in_dist_tif$datapath))

    # pdebug(devug=devug,sep='\n',pre='---H2S\n'," hs2rs_tif[]") # = = = = = = =  = = =  = = =  = = =  = = =

    rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______

    rv$tifpathdistfix <- paste0(tempFolder, '/in_dist_fixed_', rv$inDistSessID, '.tif')
    newtifPath_dist <- fitRaster2cola(inrasterpath = rv$tiforig, outrasterpath = rv$tifpathdistfix)
    newtifPath_dist <- ifelse(is.na(newtifPath_dist), yes = rv$tifpathpts, no = newtifPath_dist)
    rv$newtifPath_dist <- newtifPath_dist


    if (file.exists(rv$newtifPath_dist)){
      rv$log <- paste0(rv$log, ' --- DONE');updateVTEXT(rv$log) # _______

      rv$tifready <- TRUE
      rv$tif <- newtifPath_dist
      params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)


      output$ll_map_dist <- leaflet::renderLeaflet({
        rv$tif_sp <- newtif <- terra::rast(newtifPath_dist)
        rv$tif_rng <- rng_newtif <- range(newtif[], na.rm = TRUE)
        #rv$tif_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)

        rv$tif_pal <- tifPal <<-leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                      domain = rng_newtif, na.color = "transparent")

        makeLL( )
        # llmap <<- rv$llmap %>% removeImage('SurfaceResistance')  %>% removeControl('legendSurface') %>%
        #   addRasterImage(newtif, colors = tifPal, opacity = .7,
        #                  group = "Surface resistance", layerId = 'SurfaceResistance') %>%
        #   addLegend(pal =  tifPal, values = newtif[], layerId = "legendSurface",
        #             position = 'topleft',
        #             title= "Resistance"#, opacity = .3
        #             #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
        #   ) %>%  leaflet::addLayersControl(
        #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
        #     overlayGroups = c("Habitat suitability", "Surface resistance"),
        #     options =  leaflet::layersControlOptions(collapsed = FALSE)
        #   ) %>% clearBounds() %>%  leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" )
        #
        #
        # rv$llmap <<- llmap
        # updateLL(llmap)
        # # leafsurface
        # #llmap
        # rv$llmap
      })
    }
  })


  ##> vout_dist; ll_map_dist; distance_py; in_distance_3,
  ##> in_distance_shp in_dist_tif, inDistSessID distmap newtifPath_dist newshpPath_dist
  observeEvent(input$in_dist_shp, {

    pdebug(devug=devug,sep='\n',pre='--','names(input)', 'str(input$in_dist_shp)') # _____________

    # invisible(suppressWarnings(tryCatch(file.remove(c(in_distance_shp, newin_distance_shp)), error = function(e) NULL)))

    rv$log <- paste0(rv$log, '\nLoading shapefile');updateVTEXT(rv$log) # _______


    ## Create session IF if started from this tab
    if(is.null(rv$inSurSessID)){
      pdebug(devug=devug,sep='\n',pre='--','rv$inSurSessID') # _____________
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
      pdebug(devug=devug,sep='\n',pre='--','rv$inSurSessID') # _____________
    }


    if(is.null(rv$inDistSessID)){
      pdebug(devug=devug,sep='\n',pre='--','rv$inDistSessID') # _____________
      (inDistSessID <<- sessionIDgen()) # rv <- list()
      rv$inDistSessID <- inDistSessID
      pdebug(devug=devug,sep='\n',pre='--','rv$inDistSessID') # _____________
    }


    pdebug(devug=devug,sep='\n',pre='--','is.null(rv$distrast)') # _____________



    if(!(rv$tifready)){
      rv$log <- paste0(rv$log, '\nSTOP: load a valid surface raster first');updateVTEXT(rv$log) # _______
    } else {

      inFiles <- input$in_dist_shp #
      inFiles$newFile <- paste0(tempFolder, '/', basename(inFiles$name))
      pdebug(devug=devug,sep='\n',pre='--','(inFiles)', 'print(inFiles)') # _____________

      file.copy(inFiles$datapath, inFiles$newFile);
      # try(file.remove(inFiles$datapath))

      inShp <<- loadShp(inFiles, tempFolder, rv$inDistSessID)
      pdebug(devug=devug,sep='\n',pre='--','is.null(rv$newtifPath_dist)', 'rv$newtifPath_dist') # _____________
      rv$ptsready <- TRUE
      rv$pts <- inShp$layer
      rv$shp <- inShp$shp

      rv$log <- paste0(rv$log, '\nShapefile loaded');updateVTEXT(rv$log) # _______

      if (any(class(inShp$shp) %in% 'sf')){

        output$ll_map_dist <- leaflet::renderLeaflet({

          params_txt <- updateParamsTEXT(params_txt = params_txt, pts = TRUE)

          shp <- st_transform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          shp$sortID <- 1:nrow(shp)
          rv$pts_sp <- shp

          make(LL)
          # llmap <<- rv$llmap  %>% removeMarker(layerId = 'Points') %>%
          #   addMarkers(data = rv$pts_sp, label = ~ID, group = 'Points')
          # rv$llmap <<- llmap
          # updateLL(llmap)
          # # leafsurface
          # #llmap
          # rv$llmap
        })
      }
    }
  })

  observeEvent(input$dist_py, {
    pdebug(devug=devug,sep='\n', pre ='\nRUN CDMat\n', 'rv$tif','rv$pts','rv$tifready','rv$ptsready') # _____________

    if( ! all(rv$ptsready, rv$tifready)){ # not ready
      pdebug(devug=devug,sep='\n', pre ='\n|||', 'rv$distPy_sessID', 'outcdmat') # _____________
    } else{ # Working
      rv$distPy_sessID <- distPy_sessID <- sessionIDgen()

      outcdmat <<- paste0(tempFolder, '/out_cdmatrix_', rv$distPy_sessID, '.csv')
      # outcdmat <- '/data/temp/G2023083001195205file35162c836424/out_cdmatrix_L2023083001200105file3516441548c2.csv'

      pdebug(devug=devug,sep='\n ', pre ='\n', 'rv$distPy_sessID', 'outcdmat') # _____________

      rv$log <- paste0(rv$log, '\nGenerating matrix - ');updateVTEXT(rv$log) # _______
      Sys.sleep(1)
      tStartMat <- Sys.time()

      pdebug(devug=devug,sep='\n ', pre ='\n', 'rv$pts', 'rv$tif', 'outcdmat') # _____________
      cdmat_file <- cdmat_py (py = py, inshp = rv$pts, intif = rv$tif,
                              outcsv = outcdmat, param4 = as.numeric(input$in_dist_3))
      rv$cdm <- cdmat_file$file
      tElapMat <- Sys.time() - tStartMat
      textElapMat <- paste(round(as.numeric(tElapMat), 2), attr(tElapMat, 'units'))


      if(file.exists(outcdmat)){
        rv$cdm_sp <- headMat <- data.table::fread(outcdmat, header = F)
        rv$cdm_nvalid <- validCels <- sum(headMat > 1, na.rm = T)/2
        params_txt <- updateParamsTEXT(params_txt = params_txt, dst = TRUE)

        rv$log <- paste0(rv$log,
                         paste0(' --- DONE\nMatrix generated, dim:',
                                ncol(headMat), ' cols, ', nrow(headMat), ' rows.',
                                ' Time elapsed: ', textElapMat));updateVTEXT(rv$log) # _______
        rv$log <- paste0(rv$log,
                         paste0('\n  Pairs of points connected: ', validCels,
                                 '\n Check before sunning CDPOP. If number is low, increase the distance'));updateVTEXT(rv$log) # _______

        output$dist_box1 <- shinydashboard::renderValueBox({
          valueBox("YES", "Matrix: Done",
                   #, icon = icon("thumbs-up", lib = "glyphicon"),
                   color = "green" )
        })
      }
    }
  })


  ####### > LCC  ------------------


  observeEvent(input$in_lcc_tif, {
    pdebug(devug=devug,
           sep='\n',pre='\n---- LCC - TIF\n',
           'rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________

    # invisible(suppressWarnings(tryCatch(file.remove(c(rv$tifpathdist, rv$newtifPath_dist)), error = function(e) NULL)))

    if(is.null(rv$inLccSessID)){
      (inLccSessID <- sessionIDgen())
      rv$inLccSessID <- inLccSessID
    }

    rv$tiforig <- paste0(tempFolder, '/in_lcc_', rv$inLccSessID, '.tif')
    file.copy(input$in_lcc_tif$datapath, rv$tiforig);

    # try(file.remove(input$in_lcc_tif$datapath))

    rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______
    tifFixed <- paste0(tempFolder, '/in_lcc_fixed', rv$inLccSessID, '.tif')
    pdebug(devug=devug,sep='\n',pre='\n','input$in_lcc_tif$datapath', 'rv$tiforig', 'tifFixed') # _____________

    # rv <- list(tempFolder = '/data/temp/VK2024011517312305file1a4cf91dbd4cbd',
    #            tiforig = '/data/tempR/Rtmp4uis8l/ZAG2024051319133705//in_lcc_SQP2024051319134805.tif')
    # tifFixed <- '/data/tempR/Rtmp4uis8l/IQY2024051319031005//in_lcc_PWZ2024051319032405_fixed.tif'

    newtifPath_lcc <<- fitRaster2cola(inrasterpath = rv$tiforig,
                                      outrasterpath = tifFixed)

    rv$tif <- ifelse(is.na(newtifPath_lcc),
                     yes =  rv$tiforig,
                     no = newtifPath_lcc)

    #rv$tif <- paste0(tempFolder, '/in_lcc_fixed', rv$inLccSessID, '.tif')
    #rv$tif <- paste0(tempFolder, '/in_lcc_', rv$inLccSessID, '.tif')


    pdebug(devug=devug,sep='\n',pre='\n','newtifPath_lcc', 'rv$tif') # _____________


    if (file.exists( rv$tif)){
      rv$log <- paste0(rv$log, ' --- DONE');updateVTEXT(rv$log) # _______
      rv$tifready <- TRUE


      pdebug(devug=devug,sep='\n',pre='---- LOAD TIF LCC\n','rv$tifready', 'rv$tif', 'rv$inLccSessID') # _____________


      output$ll_map_lcc <- leaflet::renderLeaflet({

        rv$tif_sp <- terra::rast(rv$tif)
        params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)

        rv$tif_rng <- rng_newtif <- range(rv$tif_sp[], na.rm = TRUE)
        #rv$tif_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)
        rv$tif_pal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                             domain = rng_newtif+0.0, na.color = "transparent")

        makeLL( )
        # llmap <<- rv$llmap %>% removeImage('SurfaceResistance')  %>% removeControl('legendSurface') %>%
        #   addRasterImage(rv$tif_sp, colors = rv$tif_pal, opacity = .7,
        #                  group = "Surface resistance", layerId = 'SurfaceResistance') %>%
        #   addLegend(pal = rv$tif_pal, values = rv$tif_sp[], layerId = "legendSurface",
        #             position = 'topleft',
        #             title= "Resistance"#, opacity = .3
        #             #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
        #   ) %>%  leaflet::addLayersControl(
        #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
        #     overlayGroups = c("Habitat suitability", "Surface resistance"),
        #     options =  leaflet::layersControlOptions(collapsed = FALSE)
        #   ) %>% clearBounds() %>%  leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" )
        #
        #
        # rv$llmap <<- llmap
        # updateLL(llmap)
        # # leafsurface
        # #llmap
        # rv$llmap
      })
    }
  })


  observeEvent(input$in_lcc_shp, {

    # invisible(suppressWarnings(tryCatch(file.remove(c(in_lcc_shp, newin_lcc_shp)), error = function(e) NULL)))
    pdebug(devug=devug,sep='\n',pre='\n---- LCC - SHP\n','rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________

    rv$log <- paste0(rv$log, '\nLoading shapefile');updateVTEXT(rv$log) # _______

    ## Create session IF if started from this tab
    if(is.null(rv$inLccSessID)){
      (inLccSessID <<- sessionIDgen())
      rv$inLccSessID <- inLccSessID
    }

    pdebug(devug=devug,sep='\n',pre='--',
           'rv$tifready', 'inLccSessID', 'rv$inLccSessID') # _____________

    if(!(rv$tifready)){
      rv$log <- paste0(rv$log, ' -- STOP: load a valid surface raster first');updateVTEXT(rv$log) # _______
    } else {

      inFiles <- input$in_lcc_shp #

      inFiles$newFile <- paste0(tempFolder, '/', basename(inFiles$name))
      pdebug(devug=devug,sep='\n',pre='\n--','print(inFiles)', 'inFiles$newFile', 'tempFolder') # _____________

      file.copy(inFiles$datapath, inFiles$newFile);
      # try(file.remove(inFiles$datapath))
      #if(devug){save(inFiles, file = paste0(tempFolder, '/shpfiles.RData'))}

      inShp <<- loadShp(inFiles, tempFolder, rv$inlccSessID)


      if (any(class(inShp$shp) %in% 'sf')){
        # if(class(inShp$shp) == 'SpatialPointsDataFrame'){

        rv$ptsready <- TRUE
        rv$pts <- inShp$layer
        rv$shp <- inShp$shp
        rv$log <- paste0(rv$log, ' -- Shapefile loaded');updateVTEXT(rv$log) # _______
        params_txt <- updateParamsTEXT(params_txt = params_txt, pts = TRUE)

        pdebug(devug=devug,sep='\n',pre='---- LOAD SHP LCC\n','rv$ptsready', 'rv$pts', 'rv$inLccSessID') # _____________


        output$ll_map_lcc <- leaflet::renderLeaflet({

          # shp<-spTransform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          shp <- st_transform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          shp$ID <- 1:nrow( shp )
          rv$pts_sp <- shp

          makeLL( )

          # llmap <<- rv$llmap  %>% clearGroup('Points') %>%  #removeMarker(layerId = 'Points') %>%
          #   addMarkers(data = rv$pts_sp, label = ~ID, group = 'Points')
          # rv$llmap <<- llmap
          # updateLL(llmap)
          # # leafsurface
          # #llmap
          # rv$llmap
          #
          # ## try
          # llx <- makeLL()
          # rv$ll

        })
      }
    }
  })

  observeEvent(input$lcc, {
    pdebug(devug=devug,sep='\n',pre='\n---- RUN LCC\n','rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________
    condDist <<- 0
    #pdebug(devug=devug, ... = 'condDist') # _____________
    if(rv$ptsready & rv$tifready){
      condDist <<- 1
    }
    #pdebug(devug=devug, ... = 'condDist') # _____________

    #if(is.null(rv$inLccSessID)){
    (inLccSessID <<- sessionIDgen())
    rv$inLccSessID <- inLccSessID
    #}


    if( condDist == 1){
      #input <- c(in_dist_3 = 25000)
      rv$log <- paste0(rv$log, '\n Generating corridors');updateVTEXT(rv$log) # _______

      output$ll_map_lcc <- leaflet::renderLeaflet({

        out_lcc <- paste0(tempFolder, '/out_lcc_', inLccSessID, '.tif')
        tStartLcc <- Sys.time()
        #pdebug(devug=devug,sep='\n',pre='\n \t lcc.py\n', 'rv$pts', 'rv$tif', 'out_lcc', 'condDist') # _____________
        out_lcc <- lcc_py(py = py, inshp = rv$pts, intif = rv$tif, outtif = out_lcc,
                          param4 = as.numeric(input$in_lcc_4),
                          param5 = as.numeric(input$in_lcc_5),
                          param6 = as.numeric(input$in_lcc_6))

        # out_lcc <- '/data/temp/QU2024011518271005file1a4cf934de5d47/out_lcc_MQ2024011518271905file1a4cf965d2605a.tif'

        tElapLcc <- Sys.time() - tStartLcc
        textElapLcc <- paste(round(as.numeric(tElapLcc), 2), attr(tElapLcc, 'units'))

        #rv$log <- paste0(rv$log, ' - Time elapsed: ', tElapLcc);updateVTEXT(rv$log) # _______

        # pdebug(devug=devug,sep='\n',pre='\n \t LCC ','rv$out_lcc','file.exists(rv$out_lcc)', 'out_lcc') # _____________

        if(!file.exists(out_lcc$file)){
          rv$log <- paste0(rv$log, ' --- ERROR');updateVTEXT(rv$log) # _______
          rv$llmap
        } else {
          rv$log <- paste0(rv$log, ' --- DONE: ', textElapLcc);updateVTEXT(rv$log) # _______
          params_txt <- updateParamsTEXT(params_txt = params_txt, lcc = TRUE)

          rv$lcc <- out_lcc$file
          rv$lccready <- TRUE
          rv$lcc_sp <- terra::rast(out_lcc$file)

          rv$lcc_rng <- rng_newtif <- range(rv$lcc_sp[], na.rm = TRUE)
          #rv$lcc_rng <- rng_newtif <- range(minmax(rv$lcc_sp)[1:2], na.rm = TRUE)


          rv$lcc_pal <- tifPal <<- leaflet::colorNumeric(c("red3", "gold", "navyblue"),
                                                         reverse = TRUE,
                                                         domain = rv$lcc_rng+0.0,
                                                         na.color = "transparent")
          makeLL( )
        }
      })

    }
  })

  observeEvent(input$lcc2, {
    pdebug(devug=devug,sep='\n',pre='\n---- RUN LCC\n','rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________
    condDist <<- 0
    #pdebug(devug=devug, ... = 'condDist') # _____________
    if(rv$ptsready & rv$tifready){
      condDist <<- 1
    }
    #pdebug(devug=devug, ... = 'condDist') # _____________

    if(is.null(rv$inLccSessID)){
      (inLccSessID <<- sessionIDgen())
      rv$inLccSessID <- inLccSessID
    }


    if( condDist == 1){
      #input <- c(in_dist_3 = 25000)
      rv$log <- paste0(rv$log, '\n Generating corridors');updateVTEXT(rv$log) # _______

      output$ll_map_lcc <- leaflet::renderLeaflet({

        out_lcc <- paste0(tempFolder, '/out_lcc_', rv$inLccSessID, '.tif')
        tStartLcc <- Sys.time()
        #pdebug(devug=devug,sep='\n',pre='\n \t lcc.py\n', 'rv$pts', 'rv$tif', 'out_lcc', 'condDist') # _____________
        out_lcc <- lccHeavy_py(py = py, tempFolder = tempFolder,
                               inshp = rv$pts,
                               intif = rv$tif,
                               outtif = out_lcc,
                               param4 = as.numeric(input$in_lcc_4),
                               param5 = as.numeric(input$in_lcc_5),
                               param6 = as.numeric(input$in_lcc_6))

        tElapLcc <- Sys.time() - tStartLcc
        textElapLcc <- paste(round(as.numeric(tElapLcc), 2), attr(tElapLcc, 'units'))

        rv$log <- paste0(rv$log, ' - Time elapsed: ', tElapLcc);updateVTEXT(rv$log) # _______

        pdebug(devug=devug,sep='\n',pre='\n \t |||| ','rv$out_lcc','file.exists(rv$out_lcc)', 'out_lcc') # _____________

        if(!file.exists(out_lcc$file)){
          rv$log <- paste0(rv$log, ' --- ERROR');updateVTEXT(rv$log) # _______
          rv$llmap
        } else {
          rv$log <- paste0(rv$log, ' --- DONE: ', textElapLcc);updateVTEXT(rv$log) # _______

          rv$lcc <- out_lcc$file
          rv$lccready <- TRUE
          rv$lcc_sp <- terra::rast(out_lcc$file)
          params_txt <- updateParamsTEXT(params_txt = params_txt, lcc = TRUE)

          #rv$lcc_rng <- rng_newtif <- range(out_lcc[], na.rm = TRUE)
          #rv$lcc_rng <- rng_newtif <- range(minmax(rv$lcc_sp)[1:2], na.rm = TRUE)
          rv$lcc_rng <- rng_newtif <- getMxMn(rv$lcc)


          rv$lcc_pal <- tifPal <<- leaflet::colorNumeric(c("red3", "gold", "navyblue"),
                                                         reverse = TRUE,
                                                         domain = rng_newtif+0.0,
                                                         na.color = "transparent")
          makeLL( )
        }
      })

    }
  })


  ####### > CRK  ------------------


  observeEvent(input$in_crk_tif, {
    #invisible(suppressWarnings(tryCatch(file.remove(c(rv$tifpathcrk, rv$newtifPath_crk)), error = function(e) NULL)))

    pdebug(devug=devug,sep='\n',pre='---\nloadTIFCRK\n','rv$inSurSessID', 'rv$incrkSessID')
    if(is.null(rv$inSurSessID)){
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
    }

    (incrkSessID <- sessionIDgen()) # rv <- list()
    #(incrkSessID <<- sessionIDgen()) # rv <- list()
    if(is.null(rv$incrkSessID)){
      rv$incrkSessID <- incrkSessID
    }
    pdebug(devug=devug,sep='\n',pre='','rv$inSurSessID', 'rv$incrkSessID')


    rv$tiforig <- paste0(tempFolder, '/in_crk_',  rv$incrkSessID, '.tif')
    file.copy(input$in_crk_tif$datapath,rv$tiforig);
    # try(file.remove(input$in_crk_tif$datapath))

    rv$log <- paste0(rv$log, '\nUpdating raster: making pixels squared,-9999 as no data and checking coordinates systems');updateVTEXT(rv$log) # _______

    rv$tifpathcrkfix <- paste0(tempFolder, '/in_crk_fixed_', rv$incrkSessID, '.tif')
    newtifPath_crk <- fitRaster2cola(inrasterpath = rv$tiforig, outrasterpath = rv$tifpathcrkfix)
    newtifPath_crk <- ifelse(is.na(newtifPath_crk), yes = rv$tifpathpts, no = newtifPath_crk)
    rv$newtifPath_crk <- newtifPath_crk


    if (file.exists(rv$newtifPath_crk)){
      rv$log <- paste0(rv$log, ' --- DONE');updateVTEXT(rv$log) # _______
      rv$tifready <- TRUE
      rv$tif <- newtifPath_crk
      params_txt <- updateParamsTEXT(params_txt = params_txt, sr = TRUE)

      output$ll_map_crk <- leaflet::renderLeaflet({

        rv$tif_sp <- newtif <- terra::rast(newtifPath_crk)
        rv$tif_rng <- rng_newtif <- range(newtif[], na.rm = TRUE)
        #rv$lcc_rng <- rng_newtif <- range(minmax(rv$tif_sp)[1:2], na.rm = TRUE)

        rv$tif_pal <- tifPal <<- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                                       domain = rng_newtif, na.color = "transparent")

        makeLL()

      })
    }
  })


  ##> vout_crk; ll_map_crk; crkance_py; in_crkance_3,
  ##> in_crkance_shp in_crk_tif, incrkSessID crkmap newtifPath_crk newshpPath_crk
  observeEvent(input$in_crk_shp, {

    #pdebug(devug=devug,sep='\n',pre='--','names(input)', 'str(input$in_crk_shp)') # _____________

    invisible(suppressWarnings(tryCatch(file.remove(c(in_crk_shp, newin_crk_shp)), error = function(e) NULL)))

    rv$log <- paste0(rv$log, '\nLoading shapefile');updateVTEXT(rv$log) # _______

    pdebug(devug=devug,sep='\n',pre='---\nloadSHPCRK\n','rv$inSurSessID', 'rv$incrkSessID')

    ## Create session IF if started from this tab
    if(is.null(rv$inSurSessID)){
      (inSurSessID <- sessionIDgen())
      rv$inSurSessID <- inSurSessID
    }


    if(is.null(rv$incrkSessID)){
      (incrkSessID <<- sessionIDgen()) # rv <- list()
      rv$incrkSessID <- incrkSessID
    }
    pdebug(devug=devug,sep='\n',pre='\n','rv$inSurSessID', 'rv$incrkSessID', 'incrkSessID')



    if(!(rv$tifready)){
      rv$log <- paste0(rv$log, '\nSTOP: load a valid surface raster first');updateVTEXT(rv$log) # _______
    } else {

      inFiles <- input$in_crk_shp #
      inFiles$newFile <- paste0(tempFolder, '/', basename(inFiles$name))
      pdebug(devug=devug,sep='\n',pre='--','(inFiles)', 'print(inFiles)') # _____________

      file.copy(inFiles$datapath, inFiles$newFile);
      # try(file.remove(inFiles$datapath))

      inShp <<- loadShp(inFiles, tempFolder, rv$incrkSessID)
      pdebug(devug=devug,sep='\n',pre='--','is.null(rv$newtifPath_crk)', 'rv$newtifPath_crk') # _____________
      rv$ptsready <- TRUE
      rv$pts <- inShp$layer
      rv$shp <- inShp$shp

      rv$log <- paste0(rv$log, '\nShapefile loaded');updateVTEXT(rv$log) # _______

      #crksurface0 <<- crksurface ## Create bkp if new load shp

      if (any(class(inShp$shp) %in% 'sf')){

        params_txt <- updateParamsTEXT(params_txt = params_txt, pts = TRUE)

        output$ll_map_crk <- leaflet::renderLeaflet({
          shp <- st_transform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
          shp$ID <- 1:nrow( shp )
          rv$pts_sp <- shp

          makeLL()

          # llmap <<- rv$llmap  %>% clearGroup('Points') %>%
          #   addMarkers(data = rv$pts_sp, label = ~ID, group = 'Points')
          # rv$llmap <<- llmap
          # updateLL(llmap)
          # # leafsurface
          # #llmap
          # rv$llmap
        })
      }
    }
  })

  observeEvent(input$crk, {
    pdebug(devug=devug,' rv$distshp','rv$distshp', 'rv$distrast', 'inShp$files') # _____________
    condDist <- 0
    if(rv$ptsready & rv$tifready){
      condDist <- 1
    }

    #if(is.null(rv$incrkSessID)){
    (incrkSessID <<- sessionIDgen()) # rv <- list()
    rv$incrkSessID <- incrkSessID
    #}

    if( condDist == 1){
      #input <- c(in_dist_3 = 25000)
      rv$log <- paste0(rv$log, '\n Generating kernels');updateVTEXT(rv$log) # _______
      out_crk <- paste0(tempFolder, '/out_crk_', incrkSessID, '.tif')

      output$ll_map_crk <- leaflet::renderLeaflet({

        tStartCrk <- Sys.time()
        out_crk <- crk_py(py = py, inshp = rv$pts, intif = rv$tif, outtif = out_crk,
                          param4 = as.numeric(input$in_crk_4),
                          param5 = (input$in_crk_5),
                          param6 = as.numeric(input$in_crk_6))
        #out_crk_no_data <- gdal_nodata

        tElapCrk <- Sys.time() - tStartCrk
        textElapCrk <- paste(round(as.numeric(tElapCrk), 2), attr(tElapCrk, 'units'))

        rv$crk <- out_crk$file

        if(!file.exists(out_crk$file)){
          rv$log <- paste0(rv$log, ' --- ERROR');updateVTEXT(rv$log) # _______
          rv$llmap

        } else {
          rv$log <- paste0(rv$log, ' --- DONE: ', textElapCrk);updateVTEXT(rv$log) # _______

          # rv$lcc <- out_lcc
          # rv$lcc_sp <- out_lcc <- terra::rast(out_lcc)
          # out_crk <- '/data/temp//Z2023090113392605file84467aef57c/out_crk_W2023090113393905file8444afbe785.tif'
          params_txt <- updateParamsTEXT(params_txt = params_txt, crk = TRUE)

          rv$crkready <- TRUE
          rv$crk <- out_crk$file
          rv$crk_sp <- terra::rast(out_crk$file);
          rv$crk_rng <- rng_newtif <- range(rv$crk_sp[], na.rm = TRUE)
          #rv$crk_rng <- rng_newtif <- range(minmax(rv$crk_sp)[1:2], na.rm = TRUE)

          # newtif[newtif[] == 0] <-

          # newtif <- (newtif- min(rng_newtif))/(max(rng_newtif)- min(rng_newtif))
          # # plot(newtif)

          rv$crk_pal <- tifPal <<- leaflet::colorNumeric(palette = "plasma", reverse = TRUE,
                                                         domain = rng_newtif+0.01, na.color = "transparent")
          # "viridis", "magma", "inferno", or "plasma".

          makeLL()

          # llmap <<- rv$llmap %>% removeImage('Kernel')  %>% removeControl('legendKernel') %>%
          #   addRasterImage(out_crk, colors = tifPal, opacity = .7,
          #                  group = "Surface resistance", layerId = 'Kernel') %>%
          #   addLegend(pal =  tifPal, values = out_crk[], layerId = "legendKernel",
          #             position = 'topleft',
          #             title= "Dispersal Kernel"#, opacity = .3
          #             #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          #   ) %>%  leaflet::addLayersControl(
          #     overlayGroups = c('Points', "Habitat suitability", "Surface resistance", 'Corridor'),
          #     options =  leaflet::layersControlOptions(collapsed = FALSE)
          #   ) %>% clearBounds() %>%  leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" )
          #
          # rv$llmap <<- llmap
          # updateLL(llmap)
          # rv$llmap
          #
          # # ## try
          # llx <- makeLL()
          # rv$ll
        }
      })
    }
  })

  ####### > PRIORI  ------------------

  observeEvent(input$pri, {
    condDist <- 0
    if(rv$crkready & rv$lccready){
      condDist <- 1
    }

    #if(is.null(rv$inpriSessID)){
    (inpriSessID <<- sessionIDgen()) # rv <- list()
    rv$inpriSessID <- inpriSessID
    #}

    if( condDist == 1){
      #input <- c(in_dist_3 = 25000)
      rv$log <- paste0(rv$log, '\n Generating prioritization');updateVTEXT(rv$log) # _______
      out_pri_tif <- paste0(tempFolder, '/out_pri_', rv$inpriSessID, '.tif')
      out_pri_tif_patch <- paste0(tempFolder, '/out_pri_patch_', rv$inpriSessID, '.tif')
      out_pri_shp <- paste0(tempFolder, '/out_pri_', rv$inpriSessID, '.shp')
      out_pri_shp_pol <- paste0(tempFolder, '/out_pri_pol_', rv$inpriSessID, '.shp')
      out_pri_shp_patch <- paste0(tempFolder, '/out_pri_patch_', rv$inpriSessID, '.shp')


      # tempFolder <- '/data/tempR/colaZGI2024051609530805//'; list.files(path = tempFolder)
      # rv <- list(crk = '/data/tempR/colaZGI2024051609530805/out_crk_CFR2024051610192005.tif',
      #            tif = '/data/tempR/colaZGI2024051609530805/in_crk_IIF2024051609552305.tif',
      #            lcc = '/data/tempR/colaZGI2024051609530805/out_lcc_FRE2024051610211205.tif')
      #rv$inpriSessID <- 'sessA'

      #
      # out_pri_tif <- '/data/tempR/colaZGI2024051609530805/out_pri_IF2024011520212105file176c0d2f0a3994.tif'
      # out_pri_shp <- '/data/tempR/colaZGI2024051609530805/out_pri_IF2024011520212105file176c0d2f0a3994.shp'
      # input <- list(in_pri_5 = 0.5, in_lcc_6 = 50000)

      output$ll_map_pri <- leaflet::renderLeaflet({

        tStartPri <- Sys.time()
        out_pri <- tryCatch(pri_py(tif = rv$tif,
                                   incrk = rv$crk ,
                                   inlcc = rv$lcc,
                                   maskedcsname = paste0(tempFolder, '/out_pri_temp_', rv$inpriSessID, '.tif'),
                                   outshppoint = out_pri_shp,
                                   outshppol = out_pri_shp_pol,
                                   outshppatch = out_pri_shp_patch,
                                   outtif = out_pri_tif,
                                   outtifpatch = out_pri_tif_patch,
                                   param7 = as.numeric(input$in_pri_5), # 0.5
                                   param8 = as.numeric(input$in_lcc_6)), error = function(e) e)
        ## missing param7 and 8 by user


        #out_crk_no_data <- gdal_nodata

        tElapPri <- Sys.time() - tStartPri
        textElapPri <- paste(round(as.numeric(tElapPri), 2), attr(tElapPri, 'units'))

        rv$pritif <- out_pri$tif
        rv$prishp <- out_pri$shp

        if(is.na(out_pri$shp)){
          rv$log <- paste0(rv$log, ' --- ERROR ');updateVTEXT(rv$log) # _______
          return(rv$llmap)

        } else {
          rv$log <- paste0(rv$log, ' --- DONE: ', textElapPri);updateVTEXT(rv$log) # _______

          # rv$lcc <- out_lcc
          # rv$lcc_sp <- out_lcc <- terra::rast(out_lcc)

          # out_crk <- '/data/temp//Z2023090113392605file84467aef57c/out_crk_W2023090113393905file8444afbe785.tif'
          rv$priready <- TRUE
          rv$pritif_sp <- terra::rast(rv$pritif); #plot(newtif)
          rv$prishp_sp <- sf::read_sf(rv$prishp); #plot(newtif)

          rv$pri_rng <- rng_newtif <- range(rv$crk_sp[], na.rm = TRUE)
          #rv$crk_rng <- rng_newtif <- range(minmax(rv$crk_sp)[1:2], na.rm = TRUE)
          rv$pri_pal <- tifPal <<- leaflet::colorNumeric(palette = "plasma", reverse = TRUE,
                                                         domain = rng_newtif+0.01, na.color = "transparent")
          # "viridis", "magma", "inferno", or "plasma".

          makeLL()

          # llmap <<- rv$llmap %>% removeImage('Kernel')  %>% removeControl('legendKernel') %>%
          #   addRasterImage(out_crk, colors = tifPal, opacity = .7,
          #                  group = "Surface resistance", layerId = 'Kernel') %>%
          #   addLegend(pal =  tifPal, values = out_crk[], layerId = "legendKernel",
          #             position = 'topleft',
          #             title= "Dispersal Kernel"#, opacity = .3
          #             #, labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))
          #   ) %>%  leaflet::addLayersControl(
          #     overlayGroups = c('Points', "Habitat suitability", "Surface resistance", 'Corridor'),
          #     options =  leaflet::layersControlOptions(collapsed = FALSE)
          #   ) %>% clearBounds() %>%  leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" )
          #
          # rv$llmap <<- llmap
          # updateLL(llmap)
          # rv$llmap
          #
          # # ## try
          # llx <- makeLL()
          # rv$ll
        }
      })
    }
  })

  ####### > COMPARE  ------------------

  ## Load shp
  observeEvent(input$in_com_shp, {

    # invisible(suppressWarnings(tryCatch(file.remove(c(in_lcc_shp, newin_lcc_shp)), error = function(e) NULL)))
    pdebug(devug=devug,sep='\n',pre='\n---- LCC - SHP\n','rv$ptsready', 'rv$pts', 'rv$ptsready', 'rv$pts','rv$inLccSessID') # _____________

    rv$log <- paste0(rv$log, '\nLoading scenario shapefile');updateVTEXT(rv$log) # _______

    ## Create session IF if started from this tab
    if(is.null(rv$inComSessID)){
      (inComSessID <<- sessionIDgen())
      rv$inComSessID <- inComSessID
    }

    pdebug(devug=devug,sep='\n',pre='--', 'inComSessID', 'rv$inComSessID') # _____________

    inFiles <- input$in_com_shp #

    inFiles$newFile <- paste0(tempFolder, '/', basename(inFiles$name))
    pdebug(devug=devug,sep='\n',pre='\n--','print(inFiles)', 'inFiles$newFile', 'tempFolder') # _____________

    file.copy(inFiles$datapath, inFiles$newFile);
    # try(file.remove(inFiles$datapath))
    #if(devug){save(inFiles, file = paste0(tempFolder, '/shpfiles.RData'))}

    inShp <<- loadShp(inFiles, tempFolder, rv$inComSessID)
    #save(inShp, inFiles, file = 'debug_edi_shp.RData')

    if (any(class(inShp$shp) %in% 'sf')){
      # if(class(inShp$shp) == 'SpatialPointsDataFrame'){

      rv$comready <- TRUE
      rv$com <- inShp$layer
      rv$comshp <- inShp$shp
      rv$log <- paste0(rv$log, ' -- Comparisson shapefile loaded');
      updateVTEXT(rv$log) # _______

      pdebug(devug=devug,sep='\n',pre='---- LOAD SCE SHP\n','rv$comready', 'rv$com', 'rv$inComSessID') # _____________

      com <- st_transform(inShp$shp, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
      rv$com_sp <- com


      proxy <- leafletProxy("ll_map_com")
      # print(getwd())
      # save(com, file = 'com_debug_shp.RData'); print (com)
      # (load('com_debug_shp.RData'))

      # proxy <- leaflet() %>% addTiles()
      bounds <- com %>% st_bbox() %>% as.character()

      proxy %>%  #remove(layerId = 'SurfaceResistance')  %>%
        removeShape('layerId') %>%
        leaflet::addPolygons(data = com, color = 'darkblue', fillColor = 'darkblue', group = "layerId") %>%
        fitBounds(bounds[1], bounds[2], bounds[3], bounds[4])

      # output$ll_map_edi <- leaflet::renderLeaflet({
      #   makeLL( )
      # })
    }
  })


  ## Run comparisson
  observeEvent(input$com_py, {
    rv$log <- paste0(rv$log, ' \n Comparing scenarios --- ');updateVTEXT(rv$log) # _______

    in_com_ly <- input$in_com_ly
    # in_com_ly <- 'Dispersal kernels'
    # tempFolder <- '/data/temp/scenario_folder'
    layer_type_compare <- switch(in_com_ly,
                                 #'Surface resistance' = 'out_surface_.+.tif$',
                                 'Dispersal kernels' = 'out_crk_.+.tif$',
                                 'Least cost path corridos' = 'out_lcc_.+.tif$')

    # tempFolder <- '/tmp/RtmplWdZFP/colaBJM2024073022213505/';
    # layer_type_compare <- 'out_crk_.+.tif$'
    (avail_layers <- list.files(path = tempFolder, pattern = layer_type_compare,
                                full.names = TRUE))
    (avail_layers <- grep('resam.tif$', avail_layers, value = TRUE, invert = TRUE))
    # avail_layers <- rev(avail_layers)
    avail_layers <- avail_layers[order( gsub('[[:punct:]]|[a-zA-Z]', '', basename(avail_layers)) )]

    # mssg2Display <- paste0(length(avail_layers), ' layer found for ', in_com_ly, ': ', paste0(basename(avail_layers), collapse = ' '))

    (compID <- sessionIDgen(short = TRUE))
    if(in_com_ly == 'Least cost path corridos'){
      prefcomp <- '_lcc_'
    } else if(in_com_ly == 'Dispersal kernels'){
      prefcomp <- '_crk_'
    }

    (outComFolder <- paste0(tempFolder, '/comp', prefcomp, compID))
    rv$comFolder <- outComFolder
    (outComPngAbs <- paste0(outComFolder, '/compAbs.png'))
    (outComPngRel <- paste0(outComFolder, '/compRel.png'))

    (outComCsvAbs <- paste0(outComFolder, '/compAbs.csv'))
    (outComCsvRel <- paste0(outComFolder, '/compRel.csv'))

    dir.create(outComFolder, recursive = TRUE)
    # "C:/size7_crk.tif,C:/size7_s1_crk.tif,C:/size7_s2_crk.tif"

    pdebug(devug=devug,sep='\n',pre='---COMP\n', "rv$comready", 'rv$com','in_com_ly') # = = = = = = =  = = =  = = =  = = =  = = =
    inCompShp <-  'None'
    inCompShpField <-  'None'
    if (isTRUE(rv$comready) & is.null(rv$com)){
      if (file.exists(rv$com)){
        inCompShp <- rv$com
      }
    }


    ( cond <- all(file.exists(avail_layers)) & (length(avail_layers) >=2  ) )

    if (cond ){

      if(in_com_ly == 'Least cost path corridos'){

        comp_out <- tryCatch(
          lcc_compare_py(intif = avail_layers[1],
                         intifs = paste0('"', paste0(avail_layers, collapse = ','), '"'),
                         outcsvabs = outComCsvAbs,
                         outcsvrel = outComCsvRel,
                         outpngabs = outComPngAbs,
                         outpngrel = outComPngRel,
                         outfolder = outComFolder,
                         inshp = inCompShp,
                         shpfield = inCompShpField), error = function(e) NULL)

      } else if(in_com_ly == 'Dispersal kernels'){

        comp_out <- tryCatch(
          crk_compare_py(intif = avail_layers[1],
                         intifs = paste0('"', paste0(avail_layers, collapse = ','), '"'),
                         outcsvabs = outComCsvAbs,
                         outcsvrel = outComCsvRel,
                         outpngabs = outComPngAbs,
                         outpngrel = outComPngRel,
                         outfolder = outComFolder,
                         inshp = inCompShp,
                         shpfield = inCompShpField), error = function(e) NULL)

      }

      print(" --- comp_out")
      print(comp_out)
      if (!is.null(comp_out) & file.exists(comp_out$file)){

        output$ll_map_com <- leaflet::renderLeaflet({
          csvAbs <- read.csv(outComCsvAbs)
          csvAbs$val <- csvAbs[, 2]
          csvRel <- read.csv(outComCsvRel)
          csvRel$val <- csvRel[, 2]

          output$png1 <- renderImage({
            ## Following three lines CREATE a NEW image. You do not need them
            #outfile <- tempfile(fileext = '.png')
            #png(outfile, width = 400, height = 300) # Generate the PNG
            #dev.off()
            list(src = outComPngAbs, contentType = 'image/png', width = 400, height = 300,
                 alt = "Absolute values")
          }, deleteFile = TRUE)

          output$png2 <- renderImage({
            list(src = outComPngRel, contentType = 'image/png', width = 400, height = 300,
                 alt = "This is alternate text")
          }, deleteFile = TRUE)


          output$hccomp1 <- highcharter::renderHighchart({
            hcchart1 <<- highchart() %>% hc_exporting(enabled = TRUE) %>%
              hc_add_series(data = csvAbs,
                            type = "column", hcaes(x = 'Scenario', y = 'val')) %>%
              hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))
          })

          output$hccomp2 <- highcharter::renderHighchart({
            hcchart2 <<- highchart() %>% hc_exporting(enabled = TRUE) %>%
              hc_add_series(data = csvRel, name = 'Relative difference',
                            type = "column", hcaes(x = 'Scenario', y = 'val')) %>%
              hc_add_theme(hc_theme(chart = list(backgroundColor = 'white')))
          })


          #outComFolder <- '/tmp/RtmplWdZFP/colaBJM2024073022213505/compJJV2024073022250605'
          com_tifs <- list.files(outComFolder, pattern = '.tif$', full.names = TRUE)
          com_rast <- lapply(as.list(com_tifs), terra::rast)
          com_stack <- do.call(c, com_rast)
          #terra::NAflag(com_stack) <- 0
          com_rng <- terra::global(com_stack, fun="range")
          com_rng2 <- range(com_rng)

          com_pal <- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                           domain = com_rng2+0.001, na.color = "transparent")
          # "viridis", "magma", "inferno", or "plasma".

          ori_rast <- lapply(as.list(avail_layers), terra::rast)
          ori_stack <- do.call(c, ori_rast)
          ori_rng <- terra::global(ori_stack, fun="range")
          ori_rng2 <- range(ori_rng)

          ori_pal <- leaflet::colorNumeric(palette = "viridis", reverse = TRUE,
                                           domain = ori_rng2+0.001, na.color = "transparent")


          # ## OPT1:  sync 3 layers
          # lls <- leaflet::leaflet() %>% leaflet::addTiles() %>% #clearBounds() %>%
          #
          #   leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
          #   leaflet::addMeasure( position = "topright",
          #                        primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
          #                        activeColor = "#3D535D",completedColor = "#7D4479")
          #
          # llsA <- lls %>%
          #   addRasterImage(x = ori_stack[[1]], colors = ori_pal,
          #                  opacity = .7,
          #                  group = names(ori_stack)[1],
          #                  layerId = names(ori_stack)[1]) %>%
          #   addLegend(pal = ori_pal, values = ori_rng2,
          #             position = 'bottomleft', title = names(ori_stack)[1]) %>%
          #   leaflet::addLayersControl(
          #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
          #     overlayGroups = c(names(ori_stack)[1]),
          #     options =  leaflet::layersControlOptions(collapsed = FALSE))
          #
          # llsB <- lls %>%
          #   addRasterImage(x = ori_stack[[2]], colors = ori_pal,
          #                  opacity = .7,
          #                  group = names(ori_stack)[2],
          #                  layerId = names(ori_stack)[2]) %>%
          #   addLegend(pal = ori_pal, values = ori_rng2,
          #             position = 'bottomleft', title = names(ori_stack)[2]) %>%
          #   leaflet::addLayersControl(
          #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
          #     overlayGroups = c(names(ori_stack)[2]),
          #     options =  leaflet::layersControlOptions(collapsed = FALSE))
          #
          # llsC <- lls %>%
          #   addRasterImage(x = com_stack[[1]], colors = com_pal,
          #                  opacity = .7,
          #                  group = names(com_stack)[1],
          #                  layerId = names(com_stack)[1]) %>%
          #   addLegend(pal = com_pal, values = com_rng2,
          #             position = 'bottomleft', title = names(com_stack)[1]) %>%
          #   leaflet::addLayersControl(
          #     baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
          #     overlayGroups = c(names(com_stack)[1]),
          #     options =  leaflet::layersControlOptions(collapsed = FALSE))
          #
          #
          # leafsync::sync(llsA, llsB, llsC, no.initial.sync = TRUE)



          # OPT2:  all layers
          #if (FALSE) {
          llc <- leaflet::leaflet() %>% leaflet::addTiles() %>% #clearBounds() %>%

            leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
            leaflet::addMeasure( position = "topright",
                                 primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
                                 activeColor = "#3D535D",completedColor = "#7D4479")  %>%
            leaflet::addMiniMap( tiles = leaflet::providers$Esri.WorldStreetMap, toggleDisplay = TRUE)


          ## Orig layers
          for (x1 in 1:length(names(ori_stack)) ){ # x1 <- 1
            llc <- llc %>% addRasterImage(x = ori_stack[[x1]], colors = ori_pal,
                                          opacity = .7,
                                          group = names(ori_stack)[x1],
                                          layerId = names(ori_stack)[x1])
          }

          ## Compare layers
          # terra::NAflag(com_stack) <- 0
          com_pal <- leaflet::colorNumeric(palette = "magma", reverse = TRUE,
                                           domain = com_rng2+0.001, na.color = "transparent")

          for (x2 in 1:length(names(com_stack)) ){ # x2 <- 1
            #xr2 <- com_stack[[x2]]
            # terra::NAflag(xr2) <- 0
            llc <- llc %>% addRasterImage(x = com_stack[[x2]], colors = com_pal,
                                          opacity = .7,
                                          group = names(com_stack)[x2],
                                          layerId = names(com_stack)[x2])
          }

          # in_com_ly <- 'kernels'
          llc <- llc %>%
            addLegend(pal = ori_pal, values = ori_rng2,
                      group = in_com_ly, layerId = in_com_ly,
                      position = 'bottomleft', title = in_com_ly) %>%

            addLegend(pal = com_pal, values = com_rng2,
                      group = in_com_ly, layerId = 'Comparison',
                      position = 'bottomleft', title = 'Comparisson') %>%
            leaflet::addLayersControl(
              baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
              overlayGroups = c(names(com_stack), names(ori_stack)),
              options =  leaflet::layersControlOptions(collapsed = FALSE))
          # }
        })

      }

      #     if (!file.exists(points_file$file)){
      #       rv$log <- paste0(rv$log, '  --- Error creating points');updateVTEXT(rv$log) # _______
      #     } else {
      #       params_txt <- updateParamsTEXT(params_txt = params_txt, pts = TRUE)
      #
      #       rv$pts <- points_file$file
      #       rv$ptsready <- TRUE
      #
      #       output$ll_map_points <- leaflet::renderLeaflet({
      #
      #         #points_file <- "/data/temp/L2023090100204905file18e703e3d6298/out_simpts_J2023090100210305file18e7061e66c55.shp"
      #         points_shpO <- sf::read_sf(points_file$file)
      #         points_shp <- sf::st_transform(points_shpO, crs = sf::st_crs("+proj=longlat +ellps=GRS80"))
      #         points_shp$ID <- 1:nrow(points_shp)
      #         #points_shp@data[, c('lng', 'lat')] <- points_shp@coords
      #         rv$pts_sp <- points_shp
      #
      #         rv$log <- paste0(rv$log, '  --- DONE');updateVTEXT(rv$log) # _______
      #
      #         #temLL <- rv$llmap
      #         #save(temLL, file = '/data/tempR/ll.RData')
      #         #load('/data/tempR/ll.RData') # rv <- list(llmap = temLL); llmap = temLL
      #
      #         makeLL( )
      #
      #       })
      #     }
    }
  })


  observeEvent( input$in_com_ly, {
    # in_com_ly <- 'Dispersal kernels'
    # choices <- c('Surface resistance', 'Dispersal kernels', 'Least cost path corridos')
    #tempFolder <- '/data/temp/scenario_folder'

    in_com_ly <- input$in_com_ly
    layer_type_compare <- switch(in_com_ly,
                                 #'Surface resistance' = 'out_surface_.+.tif$',
                                 'Dispersal kernels' = 'out_crk_.+.tif$',
                                 'Least cost path corridos' = 'out_lcc_.+.tif$')

    # tempFolder <- '/data/temp/scenario_folder'
    avail_layers <- list.files(path = tempFolder, pattern = layer_type_compare,
                               full.names = TRUE)

    mssg2Display <- paste0(length(avail_layers), ' layer(s) found for ', in_com_ly, ': ',
                           paste0(basename(avail_layers), collapse = ' '))
    output$vout_com <- renderText({isolate( mssg2Display )})
  })

  ####### > Change session folder  ------------------
  observeEvent(input$sessPath, {

    sessInput <- input$sessInput
    newTempFolder <- paste0(dataFolder, sessInput)
    print(paste('---New session: ', newTempFolder))
    if(dir.exists(newTempFolder)){
      tempFolder <<- paste0(sessInput, '/')
      output$sessLog <- renderText({isolate( 'Session folder found' )})
    } else {
      output$sessLog <- renderText({isolate( 'No session folder found' )})
    }

  })

  ####### > Errors read  ------------------

  updateSelectizeInput(session, inputId = 'sel_error',
                       choices = validLogs,
                       selected = NA, server = TRUE)

  observeEvent(input$read_error, {
    output$outerror  <- renderText({
      txtPath <- file.path(path_error, input$sel_error)
      if(file.exists(txtPath)){
        txt2show <- read.delim(txtPath,header = FALSE)[, 1]
      }
      else {
        txt2show <- 'No files'
      }
      paste0(txt2show, collapse = '\n')
    })
    #  read.delim('/var/log/shiny-server/cola-shiny-20231011-230516-33145.log', header = FALSE)[, 1]
  })




  ####### > Priv showcase  ------------------

  observeEvent(input$in_priv_rdata, {

    output$ll_map_showPriv <- leaflet::renderLeaflet({


      llgrp <- NULL; ll_priv <- leaflet::leaflet() %>% leaflet::addTiles() %>%
        leaflet::addMeasure( position = "topright",
                             primaryLengthUnit = "kilometers", primaryAreaUnit = "sqkilometers",
                             activeColor = "#3D535D",completedColor = "#7D4479") %>%
        leaflet::addMiniMap( tiles = leaflet::providers$Esri.WorldStreetMap, toggleDisplay = TRUE)

      colPts <- colorFactor(palette = 'RdYlGn', 1:10)
      rdata2load <<- input$in_priv_rdata$datapath
      #rdata2load <- '/home/shiny/scenarios.RData'
      files2load <<- load(rdata2load)
      file.remove(input$in_priv_rdata$datapath)

      inputsPoints <- grep('poinfts', files2load, value = TRUE)
      inputsRaster <- grep('points', value = TRUE, invert = TRUE,
                           grep('userinput__', files2load, value = TRUE))
      outputsRaster <- grep('points', value = TRUE, invert = TRUE,
                            grep('useroutput__', files2load, value = TRUE))

      if(length(inputsPoints) > 0){
        for(i in 1:length(inputsPoints)){ # i = 1
          assign("pt2ll", eval(parse(text = inputsPoints[i])))

          if(class(pt2ll) == 'SpatialPointsDataFrame'){
            pt2ll <- sf::st_transform(pt2ll, crs = sf::st_crs('EPSG:4326'))
            pt2ll[, c('ln', 'lt')] <- st_coordinates(pt2ll)
            ptname <- gsub('userinput__|useroutput__', "", inputsPoints[i])

            ll_priv <- ll_priv %>%
              addCircleMarkers(lng = pt2ll$ln,
                               lat = pt2ll$lt, color = colPts(sample(1:10, 1)),
                               group = ptname, radius = 1)
            llgrp <- c(llgrp, ptname)
          }
        }
      }

      if(length(inputsRaster) > 0){
        for(i in 1:length(inputsRaster)){ # i = 1
          assign("rast2ll", eval(parse(text = inputsRaster[i])))

          if(class(rast2ll) == 'RasterLayer'){
            ptname <- gsub('userinput__|useroutput__', "", inputsRaster[i])

            r_pal <-leaflet::colorNumeric(palette = sample(c('viridis', 'magma',  'inferno',  'plasma'), 1),
                                          reverse = TRUE,
                                          domain= base::range(rast2ll[], na.rm = TRUE) + 0.0,
                                          na.color = "transparent")

            ll_priv <- ll_priv %>%
              addRasterImage(rast2ll, colors = r_pal, opacity = .7,
                             group = ptname,
                             layerId = ptname) %>%
              addLegend(pal =  r_pal, values= base::range(rast2ll[], na.rm = TRUE),
                        group = ptname, layerId = ptname,
                        position = 'bottomleft', title = ptname)

            llgrp <- c(llgrp, ptname)
          }
        }
      }

      if(length(outputsRaster) > 0){
        for(i in 1:length(outputsRaster)){ # i = 1
          assign("rast2ll", eval(parse(text = outputsRaster[i])))

          if(class(rast2ll) == 'RasterLayer'){
            ptname <- gsub('userinput__|useroutput__', "", outputsRaster[i])

            r_pal <- leaflet::colorNumeric(palette = sample(c('viridis', '
                                                     magma',
                                                              'inferno',  'plasma'), 1),
                                           reverse = TRUE,
                                           domain= base::range(rast2ll[], na.rm = TRUE) + 0.0,
                                           na.color = "transparent")

            ll_priv <- ll_priv %>%
              addRasterImage(rast2ll, colors = r_pal, opacity = .7,
                             group = ptname,
                             layerId = ptname) %>%
              addLegend(pal =  r_pal, values= base::range(rast2ll[], na.rm = TRUE),
                        group = ptname, layerId = ptname,
                        position = 'bottomleft', title = ptname)

            llgrp <- c(llgrp, ptname)
          }
        }
      }

      ll_priv <- ll_priv %>%
        leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
        leaflet::addLayersControl(baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
                                  overlayGroups = llgrp,
                                  options =  leaflet::layersControlOptions(collapsed = FALSE))

    })

    # pdebug(devug=devug,pre='\n\t Load RDATA\n', sep='\n','files2load', 'rdata2load')



  })

  #
  ####### > Download buttons  ------------------
  # out_crk$file points_file$file out_pri_tif
  {
    #rv$comFolder <- outComFolder

    ## Download cdpop
    output$cdpDwn <- downloadHandler(
      filename = paste0('cdpop',
                        #"_",
                        #ifelse(!is.null(rv$inPtsSessID), rv$inPtsSessID, rv$sessionID),
                        '.zip'),
      content = function(filename) {
        if(!is.null( rv$cdpFolder) ){
          # rv <- list(tempFolder = '/data/temp/O2023090713414105file522721b3f66/', sessionID = 'O2023090713414105file522721b3f66')
          #filename <- paste0('points_', rv$inPointsSessID , '.zip')
          zip_file <- paste0(tempfile(), '_tempCDPOP.zip')
          #'',
          #file.path(rv$tempFolder ,
          #filename)
          #)
          #zip_file <- file.path(tempdir(), filename)
          print(paste0('Making temp zip file: ', zip_file))

          com_files <- list.files(path = rv$cdpFolder,
                                  recursive = TRUE,
                                  all.files = TRUE, include.dirs = TRUE,
                                  #pattern = "out_simpts_",
                                  full.names = TRUE)
          print(paste0('incluiding: ', paste0(zip_file, collapse = ' ')))


          if (os == 'Windows'){
            zip(zipfile = zip_file, files = paste0(com_files), flags = '-r9X')
          } else {

            # the following zip method works for me in linux but substitute with whatever method working in your OS
            zip_command <<- paste("zip -j",
                                  zip_file,
                                  paste(com_files, collapse = " "))

            pdebug(devug=devug,sep='\n',pre='\n---- Write COMP\n',
                   'filename', 'zip_file') # _____________  , 'zip_command'

            system(zip_command)
          }
          # copy the zip file to the file argument
          file.copy(zip_file, filename)
          # remove all the files created
          try(file.remove(zip_file))
        }
      })




    ## Download comparision

    output$comDwn <- downloadHandler(
      filename = paste0('comparison',
                        #"_",
                        #ifelse(!is.null(rv$inPtsSessID), rv$inPtsSessID, rv$sessionID),
                        '.zip'),
      content = function(filename) {
        if(!is.null( rv$comFolder) ){
          # rv <- list(tempFolder = '/data/temp/O2023090713414105file522721b3f66/', sessionID = 'O2023090713414105file522721b3f66')
          # rv <- list(comFolder = "/mnt/c/tempRLinux/RtmpaEKrMb/colaCQT2024080119231805//comp_crk_XNM2024080119481905")
          #filename <- paste0('points_', rv$inPointsSessID , '.zip')
          zip_file <- paste0(tempfile(), '_tempComp.zip')
          #'',
          #file.path(rv$tempFolder ,
          #filename)
          #)
          #zip_file <- file.path(tempdir(), filename)
          print(paste0('Making temp zip file: ', zip_file))

          com_files <- list.files(path = rv$comFolder,
                                  all.files = TRUE, include.dirs = TRUE,
                                  #pattern = "out_simpts_",
                                  full.names = TRUE)
          otherFiles <- list.files(full.names = TRUE, include.dirs = FALSE,
                                   path = dirname(rv$comFolder),
                                   pattern = gsub('comp_|_.+', '', basename(rv$comFolder)))

          # print(paste0('incluiding: ', paste0(zip_file, collapse = ' ')))


          if (os == 'Windows'){
            zip(zipfile = zip_file, files = c(com_files, otherFiles), flags = '-r9X')
          } else {

            # the following zip method works for me in linux but substitute with whatever method working in your OS
            zip_command <<- paste("zip -j",
                                  zip_file,
                                  paste(com_files, collapse = " "))

            pdebug(devug=devug,sep='\n',pre='\n---- Write COMP\n',
                   'filename', 'zip_file') # _____________  , 'zip_command'
            system(zip_command)
          }
          # copy the zip file to the file argument
          file.copy(zip_file, filename)
          # remove all the files created
          try(file.remove(zip_file))
        }
      })



    ## Download points
    output$ptsDwn <- downloadHandler(
      filename = paste0('points',
                        #"_",
                        #ifelse(!is.null(rv$inPtsSessID), rv$inPtsSessID, rv$sessionID),
                        '.zip'),
      content = function(filename) {
        if(!is.null( rv$pts) ){
          # rv <- list(tempFolder = '/data/temp/O2023090713414105file522721b3f66/', sessionID = 'O2023090713414105file522721b3f66')
          #filename <- paste0('points_', rv$inPointsSessID , '.zip')
          zip_file <- paste0(tempfile(), '_tempPts.zip')
          #'',
          #file.path(rv$tempFolder ,
          #filename)
          #)
          #zip_file <- file.path(tempdir(), filename)
          print(paste0('Making temp zip file: ', zip_file))

          shp_files <- list.files(path = rv$tempFolder,
                                  pattern = "out_simpts_", full.names = TRUE)
          print(paste0('incluiding: ', paste0(zip_file, collapse = ' ')))


          if (os == 'Windows'){
            zip(zipfile = zip_file, files = paste0(shp_files), flags = '-r9X')
          } else {

            # the following zip method works for me in linux but substitute with whatever method working in your OS
            zip_command <<- paste("zip -j",
                                  zip_file,
                                  paste(shp_files, collapse = " "))

            pdebug(devug=devug,sep='\n',pre='\n---- WritePTS\n',
                   'filename', 'zip_file') # _____________  , 'zip_command'

            system(zip_command)
          }
          # copy the zip file to the file argument
          file.copy(zip_file, filename)
          # remove all the files created
          try(file.remove(zip_file))
        }
      })

    ## Download prio
    output$priDwn <- downloadHandler(
      filename = paste0('prioritization',
                        '',
                        #ifelse(!is.null(rv$inPriSessID), rv$inPriSessID, rv$sessionID),
                        '.zip'),
      content = function(filename) {
        if(!is.null( rv$pritif ) & !is.null( rv$prishp ) ){
          # rv <- list(tempFolder = '/data/temp/O2023090713414105file522721b3f66/', sessionID = 'O2023090713414105file522721b3f66')
          #filename <- paste0('points_', rv$inPointsSessID , '.zip')
          # zip_file <- gsub(tempdir(), '',
          #                  #file.path(rv$tempFolder ,
          #                            filename)
          #                  #)
          zip_file <- paste0(tempfile(), '_tempPrio.zip')

          #zip_file <- file.path(tempdir(), filename)


          # rv$pritif <- out_pri$tif
          # rv$prishp <- out_pri$shp
          # rv <- list(tempFolder = '/data/temp/OW2024011618275905filebfa315b6584', inpriSessID = 'IZ2024011618510605filebfa6271938d')

          shp_files <- list.files(path = rv$tempFolder,
                                  pattern = rv$inpriSessID, full.names = TRUE)

          if (os == 'Windows'){
            zip(zipfile = zip_file, shp_files, flags = '-r9X')
          } else {

            # the following zip method works for me in linux but substitute with whatever method working in your OS
            zip_command <<- paste("zip -j",
                                  zip_file,
                                  paste(shp_files, collapse = " "))

            pdebug(devug=devug,sep='\n',pre='\n---- WritePTI\n',
                   'filename', 'zip_file') # _____________  , 'zip_command'

            system(zip_command)
          }
          # copy the zip file to the file argument
          file.copy(zip_file, filename)
          # remove all the files created
          try(file.remove(zip_file))
        }
      })

    ## Download surface
    output$tifDwn <- downloadHandler(
      filename =  paste0('sur_',
                         ifelse(!is.null(rv$inSurSessID), rv$inSurSessID, rv$sessionID),
                         '.tif'),
      content = function(filename) {
        if(!is.null( rv$tif) ){
          terra::writeRaster(rv$tif_sp,
                             filename=filename,
                             #options="INTERLEAVE=BAND",
                             #format="GTiff",
                             overwrite=TRUE)
        }
      })

    ## Download csv
    output$csvDwn <- downloadHandler(
      filename =  paste0('cdmat_', rv$inDistSessID , '.csv'),
      content = function(filename) {
        if(!is.null( rv$cdm) ){
          write.csv(rv$cdm_sp,
                    file = filename, row.names = FALSE, quote = FALSE)
        }
      })

    ## Download lcc
    output$lccDwn <- downloadHandler(
      filename =  paste0('lcc_', rv$inLccSessID , '.tif'),
      content = function(filename) {
        if(!is.null( rv$lcc) ){
          terra::writeRaster(rv$lcc_sp,
                             filename=filename,
                             #options="INTERLEAVE=BAND",
                             #format="GTiff",
                             overwrite=TRUE)
        }
      })

    ## Download crk
    output$crkDwn <- downloadHandler(
      filename =  paste0('crk_', rv$inCrkSessID , '.tif'),
      content = function(filename) {
        if(!is.null( rv$crk) ){
          terra::writeRaster(rv$crk_sp,
                             filename=filename,
                             #options="INTERLEAVE=BAND",
                             #format="GTiff",
                             overwrite=TRUE)
        }
      })

    ## Download edit surface
    output$editifDwn <- downloadHandler(
      filename =  paste0('scenarioSurfRes_', rv$inEdiSessID2 , '.tif'),
      content = function(filename) {
        if(!is.null( rv$tif) ){
          terra::writeRaster(rv$tif_sp,
                             filename=filename,
                             #options="INTERLEAVE=BAND",
                             #format="GTiff",
                             overwrite=TRUE)
        }
      })
  }
}

if (FALSE){
  ## Debug
  if (TRUE){

    #leafletOutput("mymap", width = "100%", height = "100%")
    # proxy <- leafletProxy("mymap")
    # proxy %>% addCircleMarkers(lng = cmlng, lat = cmlat, group = "draw")


    # tempFolder <- '/data/temp/WH2023090803011605filecd6476238cc/'
    # setwd(tempFolder)
    # list.files()
    #
    # rv <<- list(
    #   tif_sp =terra::rast('in_points_MB2023090803012405filecd647c8b04d.tif'), # path
    #   pts_sp = sf::read_sf('out_simpts_VW2023090803012405filecd6cf4cf6d.shp'), # path
    #   lcc_sp =terra::rast('out_lcc_AL2023090803014605filecd6226e357c.tif'), # spatial object
    #   crk_sp =terra::rast('out_crk_.tif') # spatial object
    # )
    # rv$hs_sp <- 100 - rv$tif_sp
    #
    # rv$hs_rng <- range(rv$tif_sp[], na.rm = TRUE)
    # rv$tif_rng <- range(rv$tif_sp[], na.rm = TRUE)
    # rv$lcc_rng <- range(rv$lcc_sp[], na.rm = TRUE)
    # rv$crk_rng <- range(rv$crk_sp[], na.rm = TRUE)
    #
    #   ## Color pal
    # rv$hs_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
    #                           domain = rv$hs_rng, na.color = "transparent")
    # rv$tif_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
    #                            domain = rv$tif_rng, na.color = "transparent")
    # rv$lcc_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
    #                            domain = rv$lcc_rng, na.color = "transparent")
    # rv$crk_pal <-leaflet::colorNumeric(palette = "magma", reverse = TRUE,
    #                             domain = rv$crk_rng, na.color = "transparent")
    # rv$hs <- rv$tif <- rv$lcc <- rv$crk <- 1
    # makeLL()

    # rv$sessionID <- sessionID; rv$tempFolder <- tempFolder
    #
    # llmap <- leaflet::leaflet() %>%leaflet::addTiles() %>%
    #    leaflet::addLayersControl(baseGroups = c("OpenStreetMap", "Esri.WorldImagery"),
    #                    options =  leaflet::layersControlOptions(collapsed = FALSE)) %>%
    #    leaflet::addProviderTiles( "Esri.WorldImagery", group = "Esri.WorldImagery" ) %>%
    #   setView(lng = 25, lat = -21, zoom = 5) %>%
    #   leaflet.extras::addDrawToolbar(targetGroup='draw', polylineOptions = FALSE,
    #                                  rectangleOptions = FALSE, circleOptions = FALSE,
    #                                  markerOptions = FALSE, circleMarkerOptions = FALSE,
    #                                  editOptions = leaflet.extras::editToolbarOptions())
    # rv$llmap0 <- rv$llmap <- llmap
  }
}


{
  #  >> UI ---------------------------------------------------------------------------
  # https://cran.r-project.org/web/packages/dashboardthemes/vignettes/using_dashboardthemes.html

  css <- '.nav-tabs>li>a {
  font-family: "Lucida Sans", sans-serif;
  color: red;
}'

  ui <- shinydashboard::dashboardPage(
    # useShinyjs(),
    ####  title ----
    header = shinydashboard::dashboardHeader(
      title = "ConnectingLandscapes v0"
      #,enable_rightsidebar = TRUE, rightSidebarIcon = "info-circle"
    ),


    ####  sidebar ----
    sidebar =
      shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
          id = "sidebarid",
          shinydashboard::menuItem("Home", tabName = "tab_home", icon = icon("house-user")),
          #HTML(paste("Habitat suitability <>", "resistance surface", sep="<br/>"))

          # https://fontawesome.com/search?q=edit&o=r&m=free

          shinydashboard::menuItem(HTML(paste("Habitat suitability <>", "  resistance surface", sep="<br/>")),
                                   tabName = "tab_surface", icon = icon("right-left")),

          conditionalPanel( 'input.sidebarid == "tab_surface"',
                            div(style = "margin-top: -10px"),
                            textInput('name_tif_sur', label = '', value = "",
                                      width = NULL, placeholder = 'Habitat suitability name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_sur_tif', 'Load Suitability',
                                             buttonLabel = 'Search', placeholder = 'No file',
                                             accept=c('.tif'),
                                             #accept= '.zip',
                                             multiple=FALSE),
                            div(style = "margin-top: -30px")
          ),

          shinydashboard::menuItem(HTML(paste("Customize resistance surface", sep="<br/>")),
                                   tabName = "tab_edit", icon = icon("pencil")),
          conditionalPanel( 'input.sidebarid == "tab_edit"',
                            div(style = "margin-top: -10px"),
                            textInput('name_tif_edi', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_edi_tif', 'Load Resistance',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'),
                                             #accept= '.zip',
                                             multiple=FALSE),
                            div(style = "margin-top: -30px"),
                            shiny::fileInput('in_edi_shp', 'Load polygon', buttonLabel = 'Search',
                                             placeholder = 'INC SHP, DBF, SHX and PRJ ',
                                             accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                             multiple=TRUE),
                            div(style = "margin-top: -30px")
          ),

          shinydashboard::menuItem("Create source points", tabName = "tab_points", icon = icon("map-pin")),
          conditionalPanel( 'input.sidebarid == "tab_points"',
                            div(style = "margin-top: -10px"),
                            textInput('name_hs_pts', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -20px"),
                            shiny::fileInput('in_points_hs', 'Suitability TIF:',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'),
                                             #accept= '.zip',
                                             multiple=FALSE),
                            div(style = "margin-top: -30px"),
                            # ),
                            # conditionalPanel( 'input.sidebarid == "tab_points"',
                            div(style = "margin-top: -10px"),

                            textInput('name_tif_pts', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_points_tif', 'Resistance TIF:',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'),
                                             #accept= '.zip',
                                             multiple=FALSE),
                            div(style = "margin-top: -30px"),
          ),

          shinydashboard::menuItem("Cost distance matrix", tabName = "tab_distance", icon = icon("border-all")),
          conditionalPanel( 'input.sidebarid == "tab_distance"',
                            div(style = "margin-top: -10px"),
                            textInput('name_tif_dst', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_dist_tif', 'Load Resistance',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'), multiple=FALSE),
                            div(style = "margin-top: -30px"),
                            textInput('name_pts_dst', label = '', value = "",
                                      width = NULL, placeholder = 'Points name:'),
                            shiny::fileInput('in_dist_shp', 'Load points files', buttonLabel = 'Search',
                                             placeholder = 'INC SHP, DBF, SHX and PRJ ',
                                             accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                             multiple=TRUE),
                            div(style = "margin-top: -30px"),
                            #actionButton("dist_shp", "Load points!"),

          ),

          shinydashboard::menuItem("CDPOP", tabName = "tab_cdpop", icon = icon("hippo")),

          #shinydashboard::menuItem(HTML(paste("Landscape genetics", "mapping tools", sep="<br/>")),
          #          tabName = "tab_genetics", icon = icon("route")),

          shinydashboard::menuItem(HTML(paste("Connectivity", "dispersal kernels", sep="<br/>")),
                                   tabName = "tab_kernels", icon = icon("bezier-curve")),
          conditionalPanel( 'input.sidebarid == "tab_kernels"',
                            div(style = "margin-top: -10px"),
                            textInput('name_tif_crk', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_crk_tif', 'Load Resistance',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'), multiple=FALSE),
                            div(style = "margin-top: -50px"),
                            textInput('name_pts_crk', label = '', value = "",
                                      width = NULL, placeholder = 'Points name:'),
                            shiny::fileInput('in_crk_shp', 'Load points files (all)', buttonLabel = 'Search',
                                             placeholder = 'INC SHP, DBF, SHX and PRJ',
                                             accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                             multiple=TRUE),
                            div(style = "margin-top: -30px"),
                            #actionButton("dist_shp", "Load points!"),
          ),

          shinydashboard::menuItem("Connectivity - corridors", tabName = "tab_corridors", icon = icon("route")),
          conditionalPanel( 'input.sidebarid == "tab_corridors"',
                            div(style = "margin-top: -10px"),
                            textInput('name_tif_lcc', label = '', value = "",
                                      width = NULL, placeholder = 'Surface resistance name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_lcc_tif', 'Load Resistance',
                                             buttonLabel = 'Search TIF', placeholder = 'No file',
                                             accept=c('.tif'), multiple=FALSE),
                            div(style = "margin-top: -50px"),
                            textInput('name_pts_lcc', label = '', value = "",
                                      width = NULL, placeholder = 'Points name:'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_lcc_shp', 'Load points files (all)', buttonLabel = 'Search',
                                             placeholder = 'INC SHP, DBF, SHX and PRJ ',
                                             accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                             multiple=TRUE),
                            div(style = "margin-top: -30px"),
                            #actionButton("dist_shp", "Load points!"),

          ),


          #shinydashboard::menuItem("Plotting", tabName = "tab_plotting", icon = icon("image")),
          #shinydashboard::menuItem("Mapping", tabName = "tab_maping", icon = icon("map")),
          shinydashboard::menuItem("Connectivity - prioritization",
                                   tabName = "tab_priori", icon = icon("trophy")),

          shinydashboard::menuItem("Compare results", tabName = "tab_compare", icon = icon("clone")),
          conditionalPanel( 'input.sidebarid == "tab_compare"',
                            # div(style = "margin-top: -30px"),
                            div(style = "margin-top: -10px"),
                            textInput('field_shp_com', label = '', value = "",
                                      width = NULL, placeholder = 'Shapefile name'),
                            div(style = "margin-top: -10px"),
                            shiny::fileInput('in_com_shp', 'Load polygon', buttonLabel = 'Search',
                                             placeholder = 'INC SHP, DBF, SHX and PRJ ',
                                             accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                             multiple=TRUE),
                            div(style = "margin-top: -30px")
          ),

          shinydashboard::menuItem("Assign coords", tabName = "tab_coords", icon = icon("globe")),
          shinydashboard::menuItem("PDF", tabName = "tab_pdf", icon = icon("code-fork")),
          shinydashboard::menuItem("Run locally", tabName = "tab_local", icon = icon("code-fork"))
          #,
          #menuItem("Local paths", tabName = "tab_paths", icon = icon("python"))



          #shinydashboard::menuItem("Page 1", tabName = "page1"),
          # conditionalPanel(
          #   'input.sidebarid == "page1"',
          #   sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 30),
          #   selectInput("title", "Select plot title:", choices = c("Hist of x", "Histogram of x"))
          # ),

          # conditionalPanel(
          #   'input.sidebarid <> "pagex"',
          #   verbatimTextOutput("outext")
          # ),
          #shinydashboard::menuItem("Page 2", tabName = "page2")
          # tabhome tabsurface tab_points tab_distance tab_cdpop
          # tab_corridors tab_kernels tab_plotting tab_Mapping tab_priori tab_genetics tablocal
        )
      ),

    ####  body ----
    body =
      shinydashboard::dashboardBody(
        dashboardthemes::shinyDashboardThemes(
          theme = "grey_dark"
        ),
        tags$style(HTML("
    .tabbable > .nav > li > a {background-color: grey;  color:white;}
  ")),
        shinydashboard::tabItems(

          #### UI Tabs  ----

          # tab_home tab_surface tab_points tab_distance tab_cdpop
          # tab_corridors tab_kernels tab_plotting tab_Mapping tab_priori tab_genetics tablocal
          shinydashboard::tabItem(
            'tab_home',
            fluidPage(
              #includeMarkdown("md_intro.md")
              tabsetPanel(
                type = "pills",
                tabPanel(
                  "Home",
                  includeMarkdown(
                    system.file(package = 'cola', 'docs/md_intro.md')
                  )),
                tabPanel(
                  "How it works",
                  includeMarkdown(
                    system.file(package = 'cola', 'docs/md_use.md')
                  )),

                tabPanel(
                  "Performance",

                  tabsetPanel(
                    type = "pills",
                    tabPanel(
                      "Graphs",
                      br(),
                      paste('These results are the comparisson of the developed functions with existing software.',
                            'You can see the results for several scenarios (# of pixels, # of points), in terms ',
                            'of RAM (GB) and time (minutes) spent for both softwares'),
                      # factx logx logy xaxis: npix spix size | 'Total pixels', 'Side-pixels', 'Size order'
                      # " c('Least cost path', 'Kernel density', 'Distance matrix')"
                      fluidRow(
                        column(3,
                               selectInput('xaxis', 'X-axis',
                                           choices = c('Total pixels', 'Side-pixels', 'Size order'),
                                           selected = 'Total pixels', multiple = FALSE,
                                           selectize = TRUE, width = NULL)),
                        column(3,
                               selectInput('soft', 'Software:',
                                           choices =c('Least cost path', 'Kernel density', 'Distance matrix'),
                                           selected = 'Total pixels', multiple = FALSE,
                                           selectize = TRUE, width = NULL)),

                        column(2,
                               checkboxInput('factx', 'Factor X-axis', value = FALSE, width = NULL)),
                        column(2,
                               checkboxInput('logx', 'Log X-axis', value = FALSE, width = NULL)),
                        column(2,
                               checkboxInput('logy', 'Log Y-axis', value = FALSE, width = NULL))
                      ),
                      fluidRow(
                        column(6,
                               highcharter::highchartOutput('hcout1' #, height = "800px"
                               ) %>%shinycssloaders::withSpinner(color="#0dc5c1")),
                        column(6,
                               highcharter::highchartOutput("hcout2"
                                                            #, height = "800px"
                               ) %>%shinycssloaders::withSpinner(color="#0dc5c1"))
                      ),
                      fluidRow(

                      )
                    ),
                    tabPanel(
                      "Scenarios table",
                      h3(' Scenario table'),
                      fluidRow(
                        DT::dataTableOutput(outputId =  "scetable")
                      )
                    ),
                    tabPanel(
                      "Results table",
                      h3(' Full details'),
                      DT::dataTableOutput(outputId =  "perftable")
                    )
                  )
                ),
                # https://stackoverflow.com/questions/61284247/is-there-a-way-to-display-a-gif-file-in-r-shiny
                tabPanel("Showcase",
                         includeMarkdown(
                           system.file(package = 'cola', 'docs/md_showcase.md')),
                         fluidRow(
                           # img(src=file.path(rootPath, 'showcase.gif')
                           #     #, align = "left",height='250px',width='500px'
                           # )
                           fluidPage(
                             leaflet::leafletOutput("ll_map_show", height = "600px") %>%
                               shinycssloaders::withSpinner(color="#0dc5c1")
                           )
                         )
                ),
                tabPanel("ShowcasePriv",
                         fluidRow(
                           fluidRow(
                             column(6,
                                    shiny::fileInput('in_priv_rdata', 'Load your R spatial objects',
                                                     buttonLabel = 'Search RDATA', placeholder = 'No file',
                                                     accept=c('.RData'), multiple=FALSE)),
                             column(6, )),

                           # img(src=file.path(rootPath, 'showcase.gif')
                           #     #, align = "left",height='250px',width='500px'
                           # )
                           fluidPage(
                             leaflet::leafletOutput("ll_map_showPriv", height = "600px") %>%
                               shinycssloaders::withSpinner(color="#0dc5c1")
                           )
                         )
                )
              )
            )),



          #### UI HS 2 SR  ----
          shinydashboard::tabItem(
            tabName = 'tab_surface',
            fluidRow(
              column(3,
                     tags$td(style = "width: 25%", align = "top",
                             h3(' Create surface resistance',
                                style="text-align: center;vertical-align: top")
                     )),
              column(7, verbatimTextOutput("vout_h2r") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_h2r{overflow-y:scroll; max-height: 70px}"))),
              column(1,
                     actionButton("h2rsample", HTML("Load<br/>sample data"), icon = icon("upload")))
            ),

            fluidRow(
              column(2,
                     textInput("in_sur_3", "Min. val:", '0'),
                     textInput("in_sur_4", "Max. val:", '100')
              ),
              column(2,
                     textInput("in_sur_5", "Max-resistance:", '100'),
                     textInput("in_sur_6", "Shape:", '1')
              ),

              column(6,
                     fluidRow(
                       column(3,
                              textInput("in_sur_7", "No Data:", '-9999')
                       ),
                       column(6,
                              selectInput("in_pts_hs", "Source layer:", '50', choices = '')
                       ),
                       column(3,
                              tags$table(
                                style = "width: 100%", align = "left",
                                tags$tr(
                                  tags$td(style = "width: 25%", align = "center",
                                          htmlOutput(outputId = 'out_par_surA',  fill = TRUE))
                                ))
                       ),
                     ),
                     textInput('name_sur', label = 'New layer name:', value = "",
                               width = '100%', placeholder = 'NameOfNewLayertoCreate')
                     #tags$tr(tags$td(style = "width: 20%", align = "center",),
              ),

              column(2,

                     tags$table(
                       style = "width: 100%", align = "left",
                       tags$tr(
                         tags$td(style = "width: 25%", align = "center",
                                 #
                         )),
                       tags$tr(
                         tags$td(style = "width: 25%", align = "center",
                                 br(),
                                 actionButton("h2r", HTML("Get Res\nSurf"), icon = icon("play")),
                         )),
                       tags$tr(
                         tags$td(style = "width: 25%", align = "center",
                                 h6(""),
                                 downloadButton('tifDwn', 'Download')
                         ))
                     ),

              )
            ),

            # tags$tr(tags$td(style = "width: 20%", align = "center", ), ),
            #     tags$td(#style = "width: 60%",
            #       align = "center", ),
            #     tags$td(align = "left", ),
            #     tags$td(align = "left", )
            #     ) ) ),



            # fluidPage(
            #   column(2, textInput("in_sur_3", "Min-grid:", '0')),
            #   column(2, textInput("in_sur_4", "Max-grid:", '100')),
            #   column(2, textInput("in_sur_5", "Max-resistance:", '100')),
            #   column(1, textInput("in_sur_6", "Shape:", '1')),
            #   column(2, textInput("in_sur_7", "No Data:", '-9999')),
            # ),
            fluidPage(
              leaflet::leafletOutput("ll_map_h2r", height = "600px") %>%
                shinycssloaders::withSpinner(color="#0dc5c1")
            ),
            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Habitat 2 resistance parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box
          ),
          #### UI EDIT ----

          shinydashboard::tabItem(
            tabName = 'tab_edit',
            fluidRow(
              column(4, h2(' Customize surface resistance', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_edi") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_edi{overflow-y:scroll; max-height: 70px}"))
              )
            ),

            # fluidPage(
            #   column(1, numericInput("in_edi_val", label = "Value to add/replace:", value = 0),
            #          numericInput("in_edi_wid", label = "Pixel width:", value = 1)),
            #   column(3,
            #          h6(paste(
            #            #"Draw only one geometry type at the time.",
            #            #"Only last type of polygon(s) will be used.",
            #            "Use a positive or negative single value other than 0.",
            #            "Please remove existing polygons before running again. ")),
            #    checkboxInput("in_edi_che", "All pix. touched", FALSE)),
            #   column(2, selectInput("in_edi_rs", "Source layer:", '50', choices = ''),
            #          textInput('name_edi', label = 'New layer name:', value = "",
            #                    width = '100%', placeholder = 'NameOfNewLayertoCreate')),
            #   column(3, actionButton("edi", HTML("Add vals"), icon = icon("play")),
            #          br(),
            #          actionButton("rpl", HTML("Replace vals"), icon = icon("play"))
            #          ),
            #   column(1,
            #          htmlOutput(outputId = 'out_par_ediA',  fill = TRUE),
            #          downloadButton('editifDwn', 'Download')
            #
            #   )
            # ),

            fluidPage(
              column(2,
                     h6(paste("Use a positive or negative single value other than 0.",
                              "Please remove existing polygons before running again. "))),

              column(1, textInput("in_edi_val", label = "Value:", value = 0)), # to add/replace
              column(1, numericInput("in_edi_wid", label = "Pixel width:", value = 1)),
              column(1, checkboxInput("in_edi_che", "All pix. touched", FALSE)),
              #column(2, selectInput("in_edi_rs", "Source layer:", '50', choices = '')),
              # column(2, textInput('name_edi', label = 'New layer name:', value = "",
              #                     width = '100%', placeholder = 'NameOfNewLayertoCreate')),
              column(1, actionButton("edi", HTML("Add vals"), icon = icon("play"))),
              column(1, actionButton("rpl", HTML("Replace vals"), icon = icon("play"))),

              column(2,
                     tags$table(
                       style = "width: 100%", align = "left",
                       tags$tr(
                         tags$td(style = "width: 25%", align = "center",
                                 htmlOutput(outputId = 'out_par_ediA',  fill = TRUE)),
                       ))),
              column(1, downloadButton('editifDwn', 'Download'))

            ),

            fluidPage(
              leaflet::leafletOutput("ll_map_edi", height = "600px") %>%
                shinycssloaders::withSpinner(color="#0dc5c1")),

            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Editing resistance parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box
          ),

          #### UI POINTS ----
          shinydashboard::tabItem(
            tabName = 'tab_points',
            fluidRow(
              column(4, h2(' Create points', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_points") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_points{overflow-y:scroll; max-height: 70px}"))
              )
            ),

            fluidPage(
              column(1, textInput("in_points_3", "Min-grid:", '2')),
              column(1, textInput("in_points_4", "Max-grid:", '95')),
              column(1, textInput("in_points_5", "# of points:", '50')),
              column(3, selectInput("in_points_ly", "Source layer:", '50', choices = '')),
              column(2, textInput('name_pts', label = 'New layer name:', value = "",
                                  width = '100%', placeholder = 'Name new layer')),
              # column(2, textInput('name_edi', label = 'New layer name:', value = "",
              #                     width = '100%', placeholder = 'NameOfNewLayertoCreate')),
              column(1, htmlOutput(outputId = 'out_par_ptsA',  fill = TRUE)),
              column(1, actionButton("points_py", "Create points", icon = icon("play"))),
              column(width = 1, br()),
              column(1, downloadButton('ptsDwn', 'Download'))
            ),
            leaflet::leafletOutput("ll_map_points", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),

            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Create points parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box
          ),

          ##> vout_points; ll_map_points; points_py; in_points_3 -- 5


          #### UI DISTANCE ----
          shinydashboard::tabItem(
            tabName = 'tab_distance',
            fluidRow(
              column(4, h2(' Create Distance Matrix', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_dist") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_dist{overflow-y:scroll; max-height: 70px}"))
              )
            ),

            fluidRow(
              column(12,
                     fluidRow(
                       column(2,
                              tags$table(
                                style = "width: 100%", align = "left",
                                tags$tr(
                                  tags$td(style = "width: 25%", align = "center",
                                          htmlOutput(outputId = 'out_par_distA',  fill = TRUE)
                                  ),
                                  tags$td(style = "width: 25%", align = "center",
                                          htmlOutput(outputId = 'out_par_distB',  fill = TRUE))
                                ))),
                       column(3, textInput("in_dist_3", "Distance threshold (meters):", '250000')),
                       column(2, textInput('name_dst', label = 'New CSV name:', value = "",
                                           width = '100%', placeholder = 'Name new CSV')),
                       column(1, actionButton("dist_py", "Get matrix", icon = icon("play"))),
                       column(2,
                              tags$table( style = "width: 100%", align = "center",
                                          tags$tr(
                                            tags$td(style = "width: 25%", align = "center",
                                                    htmlOutput(outputId = 'out_par_distC',  fill = TRUE))
                                          ))
                       ),
                       column(1, downloadButton('csvDwn', 'Download'))
                     )
              )
              # , column(3,
              #        fluidRow(
              #          #column(9,
              #          #shinydashboard::valueBoxOutput("dist_box1") %>%shinycssloaders::withSpinner(color="#0dc5c1")
              #          #)
              #        ))
            ),
            leaflet::leafletOutput("ll_map_dist", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Distance matrix 2 resistance parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box
          ),
          ##> vout_dist; ll_map_dist; dist_py; in_distance_3, in_distance_shp in_dist_tif

          #### UI CDPOP  ----

          shinydashboard::tabItem(
            'tab_cdpop',
            fluidPage(
              #includeMarkdown("md_intro.md")
              fluidRow(
                column(3,
                       tags$td(style = "width: 25%", align = "top",
                               h3(' Population dynamics',
                                  style="text-align: center;vertical-align: top")
                       ) # , h4(' Run CDPOP module. ')
                ),
                column(9, verbatimTextOutput("vout_cdp") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                       tags$head(tags$style("#vout_cdp{overflow-y:scroll; max-height: 70px}")))
              ),

              fluidRow(
                column(width = 4,
                       tags$table(
                         style = "width: 100%", align = "left",
                         tags$tr(
                           tags$td(style = "width: 25%", align = "center",
                                   htmlOutput(outputId = 'out_par_cdpoA',  fill = TRUE)),
                           tags$td(style = "width: 25%", align = "center",
                                   htmlOutput(outputId = 'out_par_cdpoB',  fill = TRUE))
                         ))),
                column(width = 2,
                       actionButton("run_cdpop", 'Run CDPOP')),
                column(width = 4,
                       selectizeInput(inputId = 'cdpop_ans_yy', 'Year to plot:', choices = c(''), )),
                column(width = 2,
                       actionButton("mapcdpop", "Load year map"),
                       downloadButton('cdpDwn', 'Download'))
              ),

              fluidRow(
                # tableOutput("cdpop_files"),
                # shinydashboard::valueBoxOutput("cdpop_box2"),
                # actionButton("cdpop_check2", "Check files"),
                # uiOutput("out_cdpop_files"),
                DT::dataTableOutput(outputId =  "out_cdpop_filestable"),
                column(width = 6,
                  leaflet::leafletOutput("ll_map_cdp", height = "600px") %>% shinycssloaders::withSpinner(color="#0dc5c1")
                       ),
                column(width = 6,
                       tabsetPanel(
                         type = "pills",
                         tabPanel( "Population", highcharter::highchartOutput("hccdpop1")
                                   %>%shinycssloaders::withSpinner(color="#0dc5c1")),
                         tabPanel( "Alleles", highcharter::highchartOutput("hccdpop2") #, height = "800px")
                                   %>%shinycssloaders::withSpinner(color="#0dc5c1")),
                         tabPanel( "Heterozygosity", highcharter::highchartOutput("hccdpop3") #, height = "800px")
                                   %>%shinycssloaders::withSpinner(color="#0dc5c1"))
                         )
                       )
              ),

              br(),

              fillRow(
                #column(width = 12,
                shinydashboard::box(

                  width = 12, solidHeader = T, collapsible = T,
                  title = "Advanced parameters", status = "primary", collapsed = TRUE
                  ,

                  fluidRow(
                    column(width = 4,
                           #h1('Hola'),
                           fileInput("in_cdpop_xy", "XY CSV File", accept = ".csv")),
                    column(width = 4,
                           #h1('Hallo'),
                           fileInput("in_cdpop_age", "Ages CSV File", accept = ".csv")),
                    column(width = 4,
                           #h1('Hey'),
                           fileInput("in_cdpop_cd", "CDmatrix CSV File", accept = ".csv"))
                  ),

                  fluidRow(
                    column(
                      width = 3,
                      fileInput("in_cdpop_par", "Choose CSV File", accept = ".csv"),
                      checkboxInput("header", "Header", TRUE),
                    ),
                    column(
                      width = 5,
                      textInput("textbox", label="Input the value to replace:"),
                      actionButton("replacevalues", label = 'Replace values'),
                      #uiOutput("selectUI"),
                      actionButton("deleteRows", "Delete Rows"),
                      actionButton("addcolumn", "Add Column"),
                      actionButton("removecolumn", "Remove last column"),
                      actionButton("Splitcolumn", "SplitColumn"),
                      actionButton("Undo", 'Undo')
                    ),
                    column( width = 4,
                            actionButton("cdpop_check1", "Check files"),
                            shinydashboard::valueBoxOutput("cdpop_box1")
                    ),
                  ),
                  fluidRow(
                    column(width = 12,
                           DT::dataTableOutput(outputId =  "table1"),
                           plotOutput("cdpop_params"))
                  )
                )
              )
            )
          ),

          #### UI CRK ----
          shinydashboard::tabItem(
            tabName = 'tab_kernels',
            fluidRow(
              column(4, h2(' Create kernels', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_crk") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_crk{overflow-y:scroll; max-height: 70px}"))
              )
            ),
            fluidPage(
              column(1, htmlOutput(outputId = 'out_par_crkA',  fill = TRUE)),
              column(1, htmlOutput(outputId = 'out_par_crkB',  fill = TRUE)),
              column(2, textInput("in_crk_4", "Max. dispersal distance (meters):", '125000')),
              column(2, selectInput(inputId = "in_crk_5", label = "Kernel shape:",
                                    choices =  c( 'linear', 'gaussian'), # 'RH',
                                    selected = 'linear')),
              column(1, textInput("in_crk_6", "Kernel volume:", '1')),
              column(2, selectInput("in_crk_sr", "Source layer:", '', choices = '')),
              column(2, textInput('name_crk', label = 'New layer name:', value = "",
                                  width = '100%', placeholder = 'Name new layer')),
              column(1, actionButton("crk", "Get kernels", icon = icon("play")),
                     downloadButton('crkDwn', 'Download')),
            ),

            leaflet::leafletOutput("ll_map_crk", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            # ll_map_crk crk vout_crk in_crk_3 4

            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Kernels parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4, h5('Hola')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            )
          ),


          #### UI LCC ----
          shinydashboard::tabItem(
            tabName = 'tab_corridors',
            fluidRow(
              column(4, h2(' Create corridors', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_lcc") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_lcc{overflow-y:scroll; max-height: 70px}"))
              )
            ),
            fluidPage(
              column(1, htmlOutput(outputId = 'out_par_lccA',  fill = TRUE)),
              column(1, htmlOutput(outputId = 'out_par_lccB',  fill = TRUE)),

              column(1, textInput("in_lcc_4", "Max. dispersal distance (meters):", '150000')),
              column(1, textInput("in_lcc_5", "Corridor smoothing factor:", '0')),
              column(2, textInput("in_lcc_6", "Corridor tolerance (meters):", '5')),
              column(2, selectInput("in_lcc_sr", "Source layer:", '50', choices = '')),
              column(2, textInput('name_lcc', label = 'New layer name:', value = "",
                                  width = '100%', placeholder = 'Name new layer')),
              column(1, actionButton("lcc", "Get corridors", icon = icon("play")),
                     actionButton("lcc2", "Get corridors (heavy)", icon = icon("play")),
                     downloadButton('lccDwn', 'Download'))
            ),
            leaflet::leafletOutput("ll_map_lcc", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            # ll_map_corr lcc vout_corr in_lcc_3 4 5
            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Corridors parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4, h5('Hola')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            )
          ),


          #shinydashboard::tabItem('tab_plotting',
          #         h2(' Create plots'),
          #         leaflet::leafletOutput("ll_map_plot", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1")
          # ),
          #
          #shinydashboard::tabItem('tab_mapping',
          #         h2(' Create maps'),
          #         leaflet::leafletOutput("ll_map_map", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1")
          # ),

          #### UI PRIOR ----
          shinydashboard::tabItem(
            tabName = 'tab_priori',
            fluidRow(
              column(4, h2(' Priorization', style="text-align: center;")),
              column(8, verbatimTextOutput("vout_pri") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
                     tags$head(tags$style("#vout_pri{overflow-y:scroll; max-height: 70px}"))
              )
            ),
            h6(paste("Only run this tool if you have more than two isolated pathces on the kernels layers. ")),

            fluidPage(
              column(1, htmlOutput(outputId = 'out_par_prioA',  fill = TRUE)),
              column(1, htmlOutput(outputId = 'out_par_prioB',  fill = TRUE)),
              column(1, htmlOutput(outputId = 'out_par_prioC',  fill = TRUE)),
              column(2, textInput("in_pri_5", "Threshold (quantile: 0-1)", '0.5')),
              column(2, selectInput("in_pri_lcc", "Source layer:", '50', choices = '')),
              column(2, selectInput("in_pri_crk", "Source layer:", '50', choices = '')),
              column(2, textInput('name_pri', label = 'New layer name:', value = "",
                                  width = '100%', placeholder = 'Name new layer')),
              # column(4, textInput("in_pri_6", "Corridor tolerance:", '1000')),
              column(1, actionButton("pri", "Prioritize", icon = icon("play")),
                     downloadButton('priDwn', 'Download')),
            ),

            leaflet::leafletOutput("ll_map_pri", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            br(),
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Priorization parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box
          ),

          #### UI PDF ----
          shinydashboard::tabItem(
            tabName = 'tab_pdf',
            mainPanel(
              tags$div(
                class = "container",
                rowx(
                  colx(3, textInput("pdfurl", "PDF URL"))
                ),
                rowx(
                  #col(6, htmlOutput('pdfviewer')),
                  colx(6, tags$iframe(style="height:600px; width:100%",
                                      #src="http://localhost/ressources/pdf/R-Intro.pdf"
                                      #src="/home/shiny/connecting-landscapes/R/pdf_logoA.pdf"
                                      src="pdf.pdf"
                  )
                  )
                )
              )
            )
          ),

          #### UI GENETICS ----

          shinydashboard::tabItem(
            tabName = 'tab_genetics',
            h2(' Landscape genetics'),
            h6('    Comming soon ... stay tuned')
          ),

          #### UI LOCAL --------
          shinydashboard::tabItem(
            tabName = 'tab_local',
            h2(' Running this locally'),
            #h6('    Comming soon ... stay tuned'),
            includeMarkdown(
              system.file(package = 'cola', 'docs/md_cola_install.md')
            ),
            h2(' '),
            fluidRow(
              column(width = 10,
                     selectizeInput(width = "100%",
                                    "sel_error", "Select log file:",
                                    choices = list.files(path = path_error,
                                                         pattern = 'cola|connec'
                                    )),
              ),
              column(width = 2,
                     br(),
                     actionButton("read_error", HTML("Read file")))
            ),
            verbatimTextOutput("outerror"),

            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "", status = "primary", collapsed = TRUE,
              column(4, textInput("sessInput", "", '')),
              column(4, actionButton("sessPath", "Run")),
              column(4, verbatimTextOutput("sessLog"))
            )
          ),

          #### UI COMPARE ------
          shinydashboard::tabItem(
            tabName = 'tab_compare',
            fluidRow(
              column(3, h2(' Comparing results', style="text-align: center;")
              ),
              column(3,
                     selectInput("in_com_ly", label = "Layers:", selected = '',
                                 choices = c('', #'Surface resistance',
                                             'Dispersal kernels', 'Least cost path corridos'))
              ),
              column(3,
                     verbatimTextOutput("vout_com")
              ),

              column(1,
                     actionButton("com_py", "Compare", icon = icon("play"))
              ),
              column(1,
                     downloadButton('comDwn', 'Download')
              )
              # ,
              # column(8, verbatimTextOutput("vout_com") , # %>%shinycssloaders::withSpinner(color="#0dc5c1")
              #        tags$head(tags$style("#vout_crk{overflow-y:scroll; max-height: 70px}"))
              #)

            ) ,
            # fluidRow(
            #   column(2, br()),
            # ),
            # fluidPage(
            #   column(3, textInput("in_crk_4", "Distance threshold (meters x cost):", '25000')),
            #   column(3, selectInput(inputId = "in_crk_5", label = "Kernel shape:",
            #                         choices =  c( 'linear', 'gaussian'), # 'RH',
            #                         selected = 'linear')),
            #   column(3, textInput("in_crk_6", "Kernel volume (meters x cost):", '1')),
            #   column(3, actionButton("crk", "Get kernels", icon = icon("play")),
            #          downloadButton('crkDwn', 'Download')),
            # ),
            #
            #   fluidRow(
            #     column(6,
            #            tabsetPanel(
            #              type = "pills",
            #              tabPanel(title = 'Map',
            #                       #h2(' Comparing results', style="text-align: center;"),
            #                       leaflet::leafletOutput("ll_map_com", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1")
            #
            #                       #leaflet::leafletOutput("ll_map_cprA", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            #                       #leaflet::leafletOutput("ll_map_cprB", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            #              ),
            #              tabPanel(title = 'Slider',
            #                       h2(' Comparing results', style="text-align: center;")
            #                       #leaflet::leafletOutput("ll_map_cprA", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            #                       #leaflet::leafletOutput("ll_map_cprB", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            #              ),
            #
            #              tabPanel(title = 'Synced',
            #                       fluidPage(
            #                         h2(' Comparing results', style="text-align: center;")
            #                         #column(6, leaflet::leafletOutput("ll_map_cpr1", height = "600px") %>% shinycssloaders::withSpinner(color="#0dc5c1")),
            #                         #column(6, leaflet::leafletOutput("ll_map_cpr2", height = "600px") %>% shinycssloaders::withSpinner(color="#0dc5c1"))
            #                       )
            #              ))
            #     ),
            #     column(6,
            #            tabsetPanel(
            #              type = "pills",
            #              tabPanel(title = 'PNG',
            #                       imageOutput("png1"),
            #                       imageOutput("png2")),
            #              tabPanel(title = 'Chart',
            #                       highcharter::highchartOutput("hccomp1"#, height = "800px"
            #                       ) %>%shinycssloaders::withSpinner(color="#0dc5c1"),
            #                       highcharter::highchartOutput("hccomp2"
            #                                                    #, height = "800px"
            #                       ) %>%shinycssloaders::withSpinner(color="#0dc5c1")
            #
            #                       )
            #            )
            #     )
            #   )
            # ),

            fluidRow(
              # fluidRow(
              #   column(4, selectInput("in_com_rasA", "First scenario:", '', choices = c(''))),
              #   column(4, selectInput("in_com_rasB", "Second scenario:", '', choices = c(''))),
              #   column(4, selectInput("in_com_rasC", "Scenario diff:", '', choices = c('')))
              # ),
              column(6, leaflet::leafletOutput("ll_map_com", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1")),
              column(6,

                     tabsetPanel(
                       type = "pills",
                       tabPanel( "Abosolute", highcharter::highchartOutput("hccomp1")  #, height = "800px")
                                 %>% shinycssloaders::withSpinner(color="#0dc5c1")),
                       tabPanel( "Relative", highcharter::highchartOutput("hccomp2") %>%
                                   shinycssloaders::withSpinner(color="#0dc5c1")),
                       tabPanel( "Abs. PNG", imageOutput("png1")),
                       tabPanel( "Rel. PNG", imageOutput("png2"))
                     )
              )

          ),

          # , br(),
          #   fluidRow(
          #     column(3, imageOutput("png1")),
          #     column(3, imageOutput("png2")),
          #     column(3, highcharter::highchartOutput("hccomp1"#, height = "800px"
          #     ) %>%shinycssloaders::withSpinner(color="#0dc5c1")),
          #     column(3, highcharter::highchartOutput("hccomp2"
          #                                            #, height = "800px"
          #     ) %>%shinycssloaders::withSpinner(color="#0dc5c1"))
          #   )
          # )


          shinydashboard::box(
            width = 12, solidHeader = T, collapsible = T,
            title = "Comparisson info", status = "primary", collapsed = TRUE
            ,

            #fluidRow(
            column(width = 4,
                   h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
            column(width = 4, h5('Hallo')),
            column(width = 4, h5('Hello'))
            #)
          ) # close box

          )
        ,

          #### UI COORDS ----
          shinydashboard::tabItem(
            tabName = 'tab_coords',
            h2(' Assigning projection to your points or raster'),
            fluidRow(
              column(6, shiny::fileInput(
                'in_uncrs_tif', 'Load ASCII or RSG file',
                buttonLabel = 'Search', placeholder = 'No file',
                accept=c('.asc', '.rsg'), multiple=FALSE)
              ),
              column(6, shiny::fileInput(
                'in_uncrs_pts', 'Load CSV or XY file',
                buttonLabel = 'Search', placeholder = 'No file',
                accept=c('.xy', '.csv'), multiple=FALSE)
              )),
            fluidRow(
              column(3,
                     selectizeInput("sel_crs", "Select", choices = NULL), #
              ),
              column(3,
                     actionButton("coo_tif",
                                  HTML("Assign raster projection"), icon = icon("play")),
              ),
              column(3,
                     selectizeInput("sel_crs2", "Select", choices = NULL), #
              ),
              column(3,
                     actionButton("coo_pts",
                                  HTML("Assign raster proyection"), icon = icon("play")),
              ),
            ),
            # https://shiny.posit.co/r/articles/build/selectize/
            leaflet::leafletOutput("ll_coord", height = "600px") %>%shinycssloaders::withSpinner(color="#0dc5c1")
            ,
            shinydashboard::box(
              width = 12, solidHeader = T, collapsible = T,
              title = "Coordinates parameters info", status = "primary", collapsed = TRUE
              ,

              #fluidRow(
              column(width = 4,
                     h6('Min-grid *Change to Min. value
-Minimum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 Max-grid *Change to Max. value
-Maximum suitability value. This is automatically derived from the input file. To change it, type in a new value.

                 ')),
              column(width = 4, h5('Hallo')),
              column(width = 4, h5('Hello'))
              #)
            ) # close box

          )
          #shinydashboard::tabItem('tab_example',
          #         fluidPage( )
          # ),
          # tab_home tab_surface tab_points tab_distance tab_cdpop
          # tab_corridors tab_kernels tab_plotting tab_Mapping tab_priori tab_genetics tablocal
        )
      )
  )
}

## Run the APP
shinyApp(ui, server)
