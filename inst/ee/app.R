### COLA web app.
### Ivan Gonzalez - ig299@nau.edu | gonzalezgarzonivan@gmail.com
### Patrick Jantz - Patrick.Jantz@nau.edu | jantzenator@gmail.com


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
  library(shinyalert)
  library(shinyBS)
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
(cat('\n\n >>>> Welcome to CoLa 2.2.1 \n >>>> getwd(): ', getwd(), '\n'))


## Init A
{
  ## Initials ---
  
  os <- Sys.info()[c("sysname")]
  os <<- os
  
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
  
  devug <<- FALSE
  
  (showcasePath <<- base::paste0(rootPath, '/sampledata')); dir.exists(showcasePath)
  
  uper <- cola::uper ## Unique performance table
  per <- cola::per ## Performance table
  crs_df <- cola::crs_df ## CRS available
  
  
  ## Showcase -----
  sh_object <- base::paste0(rootPath, '/docs/showcase/showcase.RData')
  
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



# >> SERVER ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  shinyalert(html = TRUE, #type = "info",
             # imageHeight = 1000,
             imageUrl = 'https://github.com/connectingLandscapes/cola/blob/main/inst/docs/logoA_bgNA.jpg?raw=true',
             #system.file(package = 'cola', 'docs/logoA_bgNA.jpg'),
             
             title = paste0("Welcome to CoLa<br><br>", sessionID ),
             text = paste0('Please save this sessionID in your records')
  )
  
  shinyjs::disable("ee_fc_upload")
  
  eeparamscov <- read.csv('C:/cola/params_covs.csv', row.names = 1)
  eeparamsccd <- read.csv('C:/cola/params_ccdc.csv', row.names = 1)
  
  output$table2 <- DT::renderDataTable(
    dat <- datatable(eeparamscov, editable = TRUE,
                     options = list(
                       paging = FALSE , pageLength = nrow(eeparamscov)
                     )))
  
  output$ee_extcovs_table <- DT::renderDataTable(
    if (exists('out_eeextcovstable')){
      dat <- datatable(out_eeextcovstable,
                       options = list( paging =TRUE,
                                       pageLength = nrow(out_eeextcovstable)
                       ))
    } else {
      dat <- data.frame()
    }
  )
  
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
  
  isolate(observeEvent(input$run_cdpop, {
    
    # tempFolder = '/tmp/RtmpYiPPnn/colaGBF2024072213145005'; setwd(tempFolder)
    # rv <- list(pts = 'out_simpts_MFU2024072213183205.shp')
    
    # tempFolder = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305'; setwd(tempFolder)
    # rv <- list(pts = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305/out_simpts_ZFE2024080118415505.shp',
    #      cdm = '/mnt/c/tempRLinux/RtmpaEKrMb/colaKZF2024080118414305/out_cdmatrix_XUV2024080118421705.csv')
    # pref <- 'SLW'
    
    ## no params provided
    if (is.null(rv$data)){
      cat('   // No existing CDPOP invars, using sample file\n')
      cdpop_invars <- read.csv(system.file(package = 'cola', 'sampledata/invars.csv'))
    } else {
      cat('   // Existing CDPOP invars\n')
      cdpop_invars <- as.data.frame(t(rv$data))
      colnames(cdpop_invars) <- colnames(read.csv(system.file(package = 'cola', 'sampledata/invars.csv')))
    }
  })
  
  output$table3 <- DT::renderDataTable(
    dat <- datatable(eeparamsccd, editable = TRUE,
                     options = list(
                       paging = FALSE , pageLength = nrow(eeparamsccd)
                     )))
  
  
  ## Upload PTS ---------
  isolate(observeEvent(input$ee_fc_upload, {
    
    (py <- Sys.getenv("COLA_PYTHON_PATH"))
    (ee_scr_path <- system.file(package = 'cola', 'sat_ts_fusion'))
    if(ee_scr_path == ''){
      (ee_scr_path <- ('C:/Users/gonza/AppData/Local/R/win-library/4.5/cola/sat_ts_fusion'))
    }
    file.exists(ee_scr_path)
    # param1 = sys.argv[1] # project name
    # param2 = sys.argv[2] # shapefile path
    # param3 = sys.argv[3] # ee asset
    
    # input <- list(ee_project = 'gonzalezivan', local_file = 'C:/cola/Anoa/Anoa_present_ardianti.shp', ee_pts = 'cola/anoa')
    
    
    cmdee <- paste0(cola::adaptFilePath(py), ' ', ee_scr_path,'/ee_uploadFeature.py ', 
                    input$ee_project, ' ',
                    cola::adaptFilePath(input$local_file), ' ', input$ee_ptspath, ' 2>&1')
    #cmdee <- '/home/shiny/.local/share/r-miniconda/envs/cola/bin/python /srv/shiny-server/cola2/ee_connect.py gonzalezivan colaHRI2025081304123905 2>&1'
    # C:\\Users\\gonza\\AppData\\Local\\r-miniconda\\envs\\cola\\python.exe C:\\cola\\cola2\\ee_connectEE.py C:\\cola\\colaHRI202508130412390.csv
    cat(' Uploading points EE:\n')
    cat(cmdee, '\n')
    
    intCMD <- tryCatch(
      capture.output(
        system( cmdee , ignore.stdout = FALSE, 
                ignore.stderr = FALSE,
                intern = TRUE)),
      error = function(e) e$message)
    
    cat(intCMD)
    
    cond <- !any(grep('ERROR', intCMD))
    
    if(cond){
      taskid <- grep('Upload submitted', intCMD, value = TRUE)
      
      shinyalert(html = TRUE, type = "success",
                 title = paste0("Earth engine connected"),
                 text = paste0(' Task ID submitted. Please check your earth engine console <br>', taskid)
      )
      
    } else {
      shinyalert(html = TRUE, type = "error",
                 title = paste0("Task  not submitted"),
                 text = paste0(intCMD)
      )
    }
    #  ee.Authenticate()
    # ee.Initialize(project='gonzalezivan')
    #
  }
  ))
  
  ## EE Connect ---------
  isolate(observeEvent(input$ee_connect, {
    
    (py <- cola::adaptFilePath(Sys.getenv("COLA_PYTHON_PATH")))
    # (py <- 'C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe')
    # (ee_scr_path <- system.file(package = 'cola', 'sat_ts_fusion'))
    (ee_scr_path <- ('C:/cola/cola2'))
    
    # file.exists(ee_scr_path)
    
    # param1 = sys.argv[1] # project name
    # param2 = sys.argv[2] # shapefile path
    # param3 = sys.argv[3] # ee asset
    # input <- list(ee_project = 'gonzalezivan')
    
    cmdee <- paste0(py, ' ', ee_scr_path,'/cml_connectEE.py ', 
                    input$ee_project, ' ', tempFolder, ' ')
    #cmdee <- '/home/shiny/.local/share/r-miniconda/envs/cola/bin/python /srv/shiny-server/cola2/ee_connect.py gonzalezivan /data/cola/colaXGY2026030209460105 2>&1'
    # cmdee <- 'C:/Users/gonza/AppData/Local/r-miniconda/envs/cola/python.exe C:/cola/cola2/cml_connectEE.py gonzalezivan C:/cola/colaIJD2026030509435605/ 2>&1'
    cat(' Cola2: Connecting to EE\n')
    cat(cmdee, '\n')
    
    intCMD <- tryCatch(
      #capture.output(
      system( cmdee , ignore.stdout = FALSE, 
              ignore.stderr = FALSE,
              intern = TRUE)
      #)
      ,
      error = function(e) e$message)
    
    print(intCMD)
    
    cond <- !any(grep('ERROR', intCMD))
    
    if(cond){
      output$box1 <- shinydashboard::renderValueBox({
        valueBox(width = 12, "EE connected", "Ready", icon = icon("thumbs-up", lib = "glyphicon"),
                 color = "green"
        )
      })
      
      shinyalert(html = TRUE, type = "success",
                 title = paste0("Earth engine connected"),
                 text = paste0(input$ee_project, ' project found')
      )
      
    } else {
      shinyalert(html = TRUE, type = "error",
                 title = paste0("Earth engine not connected"),
                 text = paste0(intCMD)
      )
    }
    #  ee.Authenticate()
    # ee.Initialize(project='gonzalezivan')
    #
  }
  ))
  
  output$box1 <-  shinydashboard::renderValueBox({
    valueBox( "EE not conected", 'Not ready', icon = icon("thumbs-down", lib = "glyphicon"),
              # icon = icon("thumbs-up", lib = "glyphicon")
              color = "red"
    )
  })
}



# >> UI ---------------------------------------------------------------------------
# https://cran.r-project.org/web/packages/dashboardthemes/vignettes/using_dashboardthemes.html

css <- '.nav-tabs>li>a {
 font-family: "Lucida Sans", sans-serif;
 color: red;
}'


ui <- shinydashboard::dashboardPage(
  
  ## Activate enable / disable buttons
  shinyjs::useShinyjs(),
  
  
  
  #### title ----
  header = shinydashboard::dashboardHeader(
    title = "ConnectingLandscapes v.1"
    #,enable_rightsidebar = TRUE, rightSidebarIcon = "info-circle"
  ),
  
  
  #### SIDEBAR ----
  sidebar =
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        id = "sidebarid",
        #HTML(paste("Habitat suitability <>", "resistance surface", sep="<br/>"))
        shinydashboard::menuItem("EE", tabName = "tab_ee", icon = icon("map"))
        
        #shinydashboard::menuItem("Page 2", tabName = "page2")
        # tabhome tabsurface tab_points tab_distance tab_cdpop
        # tab_corridors tab_kernels tab_plotting tab_Mapping tab_priori tab_genetics tablocal
      )
    ),
  
  #### BODY ----
  body =
    shinydashboard::dashboardBody(
      tags$style(
        HTML('
             #buttons {
             background-color:yellow; position:fixed; margin-bottom:50px; opacity:1; height:50px; z-index:5;
             display: flex;
             align-items: center;
             justify-content: center;
             }

             #fluidrow1 {
            height:50px;
             }
             ')
      ),
      
      dashboardthemes::shinyDashboardThemes( theme = "grey_dark"),
      tags$style(HTML(".tabbable > .nav > li > a {background-color: grey; color:white;}")),
      
      shinydashboard::tabItems(
        
        ### TAB EE --------
        shinydashboard::tabItem(
          'tab_ee',
          fluidPage(
            tabsetPanel(
              type = "pills",
              tabPanel(
                "Check + Upload EE",
                
                
                h2(' Check EE '),
                  shinydashboard::box( # open box ABC
                    width = 12, solidHeader = T, collapsible = T,
                    title = "Connect to Earth Engine", status = "info", collapsed = FALSE
                    ,
                    
                    fluidRow(
                      column(width = 6,
                             textInput(width = "100%",
                                       value = 'gonzalezivan',
                                       placeholder = 'EE project',
                                       label =  'EE project',
                                       'ee_project'),
                             
                             actionButton(width = "100%",
                                          label = 'Check EE',
                                          'ee_connect')),
                      column(width = 6,
                             shinydashboard::valueBoxOutput("box1", width = 12))
                    ) # end box ABC
                ), 
                
                
                # verbatimTextOutput("vout_ee"),
                
                

                  
                  h2('  Upload Points '),
                  
                  shinydashboard::box( # open box ABC
                    width = 12, solidHeader = T, collapsible = T,
                    title = "Create feature collection and bounding box", status = "info", collapsed = FALSE
                    ,
                    column(width = 6,
                           textInput(width = "100%",
                                     #value = 'projects/gonzalezivan/assets/cola/name',
                                     placeholder = 'projects/USER/assets/LAYER',
                                     label =  'Input points', inputId = 'ee_ptspath')),
                    
                    column(width = 3, 
                           shiny::fileInput('local_file', 'Load points files', buttonLabel = 'Search',
                                            placeholder = 'INC SHP, DBF, SHX and PRJ ',
                                            accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj", '.zip', '.gpkg', '.SQLite', '.GeoJSON', '.csv', '.xy'),
                                            multiple=TRUE)
                    ),
                    
                    column(width = 3,
                           tags$td(style = "width: 25%", align = "top",
                                   actionButton(width = "100%", class = "btn-primary btn-lg", 
                                                label = 'Upload points', 'ee_fc_upload',
                                                style="text-align: center;vertical-align: center")
                           )
                    )
                ),
                
                
                br(),
                h2('Extract covariates '),
                
                shinydashboard::box( # open box ABC
                  width = 12, solidHeader = T, collapsible = T,
                  title = "Connect to Earth Engine", status = "info", collapsed = FALSE
                  ,
                  
                  column(width = 4,
                         fileInput("in_eeex_par", "Params CSV File", accept = ".csv")),
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Load params',
                                      'ee_load_covs')),
                  
                  column(width = 6,
                         selectizeInput(
                           multiple = TRUE,
                           choices = c(
                             'True desert', 
                             'Semi-arid', 
                             'Dense short vegetation', 
                             'Tree cover', 
                             'Wetland Salt pan', 
                             'Wetland Sparse vegetation', 
                             'Wetland Dense short vegetation', 
                             'Wetland Tree cover', 
                             'Cropland', 
                             'Built-up'	
                           ),
                           selected = c(
                             'True desert', 
                             'Semi-arid', 
                             'Dense short vegetation', 
                             'Tree cover', 
                             'Wetland Salt pan', 
                             'Wetland Sparse vegetation', 
                             'Wetland Dense short vegetation', 
                             'Wetland Tree cover', 
                             'Cropland', 
                             'Built-up'	
                           ),
                           label =  'Land cover types',
                           inputId = 'ee_cov_lc')),
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Load sample table',
                                      'ee_load_covssample')),
                  
                  DT::dataTableOutput(outputId = "ee_extcovs_table")
                  
                ),
                
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Upload model CSV',
                                      'ee_push_params')),
                  column(width = 4,
                           selectizeInput(
                             label = 'Algorithm', 
                             choices = c('CCDC', 'Landtrndr'),
                             selected =  'CCDC',
                             inputId = 'ee_ccdc_type')),
                
                br(),
                h2('Extract covariates '),

                
                fluidRow(

                  
                  
                  column(width = 6,
                         selectizeInput(
                           multiple = TRUE,
                           choices = c(
                             'blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL', 'QA_RADSAT'
                           ),
                           label =  't Mmask Bands',
                           inputId = 'ee_cov_list'))
                ),
                
                shinydashboard::box( # open box ABC
                  width = 12, solidHeader = T, collapsible = T,
                  title = "Extract covariates", status = "info", collapsed = TRUE
                  ,
                  fluidRow(
                    column(
                      width = 4,
                      actionButton("in_eeext_pardef", "Load default parameters "),
                    ),
                    
                    column(
                      width = 4,
                      actionButton("in_eeext_save", "Save edited params file"),
                    )
                  ),
                  fluidRow(
                    column(width = 12,
                           DT::dataTableOutput(outputId = "table2"),
                           plotOutput("eesext_params"))
                  )
                ), # end box ABC
                
                h2('Run model'),
                fluidRow(
                  column(width = 4,
                         selectizeInput(
                           label = 'Algorithm',
                           choices = c(
                             'CCDC', 'Landtrndr'
                           ),
                           selected =  'CCDC',
                           inputId = 'ee_ccdc_type')),
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Get model CSV',
                                      'ee_get_csv')),
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Upload model CSV',
                                      'ee_push_params'))
                ),
                shinydashboard::box( # open box
                  width = 12, solidHeader = T, collapsible = T,
                  title = "Run CCDC ", status = "info", collapsed = TRUE
                  ,
                  fluidRow(
                    column(
                      width = 4,
                      actionButton("in_eeccd_pardef", "Load default parameters "),
                    ),
                    column(
                      width = 4,
                      fileInput("in_eeccd_par", "Params CSV File", accept = ".csv"),
                    ),
                    column(
                      width = 4,
                      actionButton("in_eeccd_save", "Save edited params file"),
                    )
                  ),
                  fluidRow(
                    column(width = 12,
                           DT::dataTableOutput(outputId = "table3"),
                           plotOutput("eeccd_params"))
                  )
                ), # end box
                
                
                br(),
                fluidRow(
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Check covariates',
                                      'ee_check_covar'))
                ),
                
                
              ), # tab panel
              
              tabPanel(
                "Make params",
                h2(' ')
                
              ), # tab panel
              
              tabPanel(
                "Select covariates",
                h2(' ')
                
              ), # tab panel
              
              tabPanel(
                "Model params",
                h2(' ')
                ,
                fluidRow(
                  column(width = 4,
                         selectizeInput(
                           label = 'Algorithm',
                           choices = c(
                             'RF', 'NN'
                           ),
                           selected =  'RF',
                           inputId = 'ee_model_type')),
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Check params',
                                      'ee_check_params'))
                )
              ), # tab panel
              
              tabPanel(
                "Run model",
                h2(' ')
                , fluidRow(
                  column(width = 4,
                         selectizeInput(
                           label = 'Model type',
                           choices = c(
                             'RF', 'NN'
                           ),
                           selected =  'RF',
                           inputId = 'ee_model_type')),
                  
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Check inputs',
                                      'ee_check_model')),
                  column(width = 4,
                         actionButton(width = "100%",
                                      label = 'Run model',
                                      'ee_run_model'))
                )
              ) # tab panel
              
            )
          )
        )
      )
    )
)


## Run the APP
shinyApp(ui, server)
