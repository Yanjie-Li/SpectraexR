list.of.packages <- c("librarian","Biobase","shinydashboard","tictoc",'BiocManager','quickPlot')
 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "https://cloud.r-project.org")

if (!require("EBImage", quietly = TRUE))
 BiocManager::install('EBImage')

library("EBImage")
library("tictoc")
tic()
librarian::shelf("shinythemes",'RCSF','DT',"shinydashboard",'stars', 'sfheaders', 'sf','exactextractr', 'lidR', "shiny",'tidyverse','RStoolbox','viridis',
                 'raster','rdrop2','tools','rasterVis','data.table',quiet = T)
toc()

tic()

# write_rds(data,'G:/data.rds')
# data <- readr::read_rds('e:/My package/shiny/data_for shiny.rds')
data <- readr::read_rds('data.rds')
# saveRDS(data,'data2.rds')
toc()
print('data read time')
# data  <- data %>%
#   drop_na(Chl) 
#   group_by(Fam,month)%>%  
#   dplyr::select(8,9,10,13,1:7,11,12)
raster_read  <- function(url) {
  lapply(url, function(urll){
    imag <- list.files(
      path = urll,
      pattern = '.tif',
      all.files = T,
      full.names = T,
      no.. = T
    )
    imag <- list(imag)
    # imag1 <-imag[-c(2,9)]
    lapply(imag, function(z)
      expr <- tryCatch({
        library(raster)
        # p_dsm <- raster(z[[1]])
        p_blue <- raster(z[[1]])
        p_green <- raster(z[[2]])
        p_red <- raster(z[[4]])
        p_redege <- raster(z[[5]])
        p_nir <- raster(z[[3]])
        # extent(p_dsm) <- extent(p_blue)
        # rc1 <- extend(p_nir,p_dsm )
        # p_dsm <- raster::resample(p_dsm, p_nir,method = 'ngb')
        # tictoc::tic()
        # rrrr <- gdal_resample(p_dsm, p_blue)
        # tictoc::toc()
        
        trainx <- list(p_red,p_blue,p_green,p_redege,p_nir )
        # names(trainx) <- c('red',"blue", "green",'redege','nir','dsm')
        return(trainx)
        
      },
      error = function(e) {
        message('Caught an error!')
        cat("ERROR :", conditionMessage(e), "\n")
        print(e)
      }, 
      print(paste("processing Loop", z, sep = "_"))))})}

###################


multi_rasl <- function(las_list,dsf_list,kwsindice){
  lapply(las_list, function(x){
    expr <- tryCatch({
      library("lidR")
      library("rgdal")
      library(sfheaders)
      library(stars)
      library(raster)
      library(tidyverse)
      library(sf)
      library(data.table)
      message( paste0(names(x@header@VLR$Extra_Bytes$`Extra Bytes Description`)))
      tictoc:: tic("processing las file")
      las = classify_ground(x,csf(cloth_resolution = 0.5, class_threshold = 0.15, rigidness = 1), last_returns = FALSE)
      # pologon1 <- readRDS('e:/OneDrive - Business/宋钊颖/pologon1.rds')
      # proj4string(pologon1) <- CRS("+init=epsg:32650")
      subset = las#lidR::clip_roi(las,pologon1)
      dtm <- grid_terrain(subset, res = 0.5, algorithm = knnidw(k=5,p = 0.5), use_class = c(1L, 2L),keep_lowest = F)
      nlas_dtm  <- subset - dtm
      chm    <- grid_canopy(nlas_dtm , res =0.5, p2r())
      ttops  <- find_trees(nlas_dtm , lmf(ws= kwsindice , hmin=2.6, shape = "circular"))
      algo   <- dalponte2016(chm, ttops )
      lass   <- segment_trees(nlas_dtm, algo, uniqueness ='incremental')
      crown_polo  <- delineate_crowns(lass, func = .stdtreemetrics)
      directions_longlat <- spTransform(crown_polo,sp:: CRS("+proj=longlat +datum=WGS84 +no_defs"))
      sf  <- st_as_sf(directions_longlat)
      sf33 <- as.data.frame(crown_polo)
      sf33 <- sf33%>%dplyr:: select(treeID,ZTOP,convhull_area)%>% as.data.frame()
      dtm2 <- projectRaster(dtm, crs='+proj=longlat +datum=WGS84 +no_defs')
      chm2 <- projectRaster(chm, crs='+proj=longlat +datum=WGS84 +no_defs')
      
      library(data.table) 
      fd1 <- dsf_list[names(dsf_list) ==names(x@header@VLR$Extra_Bytes$`Extra Bytes Description`)]
    
      chm23 <-  terra::resample(chm2,fd1[[1]][[1]], method = 'ngb')
       
      # plot(chm23)
      
      
      library(exactextractr)
      prec_chm <- exactextractr::exact_extract(chm23, sf, include_xy=T) %>%
        setNames(sf$treeID ) %>%
        invoke(rbind,.)%>%
        dplyr:: select(1:3) %>%
        as.data.frame()
      names(prec_chm)[1] <- 'chm'
      prec_chm <-  prec_chm %>%
        mutate(treeID =  sapply(strsplit( rownames(prec_chm),'[.]'), function(x){
          y=x[1]
        }))
      # tictoc::toc()
      # message(paste(names(fd1)))
      # tictoc::tic('spectra extraction')
      tesst <-  lapply(fd1, function(x11){
        x22 <-list(unlist(x11))
        lapply(x22, function(ls11){
          lapply(ls11, function(ls222){
            tryCatch({
              library(tictoc)
              tic("for loop start")
              message( paste0(names(ls222)))
              print(ls222)
              plot(ls222)
              plot(sf,add=T, alpha=0.6,col=rainbow(1))
              library(exactextractr)
              prec_dfs <- exactextractr::exact_extract(ls222, sf, include_xy=T) %>%
                setNames(sf$treeID ) %>%
                invoke(rbind,.)%>%
                dplyr:: select(1) %>%
                as.data.frame()
              names(prec_dfs)  <- names(ls222)
              print("finished")
              toc()
              return(prec_dfs)
            },
            error = function(e) {
              message('Caught an error!')
              cat("ERROR :", conditionMessage(e), "\n")
              print(e)}

            )})
        })
      }) %>% 
        # map_dfr(cbind) %>% 
        unlist(recursive = F)  %>%as.data.frame()
      tictoc::toc()
   names(tesst)
      dat_tes <- cbind(tesst,prec_chm)
      # dat_tes <- as.data.frame(dat_tes)
      # library(tictoc)
      # tic("Extract fam")
      # library(exactextractr)
      # fam_can <- exact_extract(canopy1, sf, include_xy=T) %>%
      #   setNames(sf$treeID ) %>%
      #   invoke(rbind,.)%>% dplyr:: select(1) %>%
      #   as.data.frame()
      # # end timer for second subsection
      # toc()
      # tic("Extract fam2")
      # fam_can$treeID <- sapply(strsplit(rownames(fam_can),'[.]'), function(x){
      #   y=x[1]
      # })
      # toc()
      # tic("Extract fam3")
      # fam_can2 <- exact_extract(canopy2, sf, include_xy=T) %>%
      #   setNames(sf$treeID ) %>%
      #   invoke(rbind,.)%>% dplyr:: select(1) %>%
      #   as.data.frame()
      # fam_can2$treeID <- sapply(strsplit(rownames(fam_can2),'[.]'), function(x){
      #   y=x[1]
      # })
      # toc()
      # tic("Extract fam4")
      # combin1 <- fam_can%>% inner_join(matou_1, by = c ('value'= 'num')) %>% as_tibble()
      # combin1 <- combin1 %>% dplyr::select(1,2,3,4,5,8,9)
      # names(combin1)[3] <- 'Fam'
      # toc()
      # tic("Extract fam5")
      # trid_fam1 <- combin1   %>%
      #   group_by(Fam ) %>%
      #   summarise_at(vars(c(treeID)), funs(min(., na.rm=TRUE)))
      # mon11_fam1 <- dat_tes %>% left_join(trid_fam1, by = c('treeID'))
      # mon11_fam1$treeID <- as.numeric(mon11_fam1$treeID)
      # toc()
      # tic("Extract fam6")
      # combin2 <- fam_can2%>% inner_join(matou_2, by = c ('value'= 'num')) %>% as_tibble()
      # combin2 <- combin2 %>% dplyr::select(1,2,3,4,5,8,9)
      #
      # names(combin2)[3] <- 'Fam'
      # trid_fam <- combin2   %>%
      #   group_by(Fam ) %>%
      #   summarise_at(vars(c(treeID)), funs(min(., na.rm=TRUE)))
      # mon11_fam2 <- dat_tes %>% left_join(trid_fam, by = c('treeID'))
      # mon11_fam2$treeID <- as.numeric(mon11_fam2$treeID)
      # mon11_fam  <- rbind.data.frame(mon11_fam1,mon11_fam2)
      # mon11_fam <- mon11_fam %>% left_join(sf33, by = c('treeID'))
      # mon11_fam <-mon11_fam[!is.na(mon11_fam$Fam),]
      # print("finished")
      # toc()
      # mon11_fam

      return(dat_tes)
      
    },error = function(e) {
      message('Caught an error!')
      cat("ERROR :", conditionMessage(e), "\n")
      print(e)},
    print(paste("processing Loop", names(las_list), sep = "_"))) 
  })
}

 

ras_im_alin <- function(monthi,fami){
  lapply(monthi, function(m){
    m2 <-list(unlist(m))
    lapply(m2, function(m22){
      expr <- tryCatch({ # message( paste0(x,'_',f))
        # lapply(monthi , function(x){
        library(tidyverse)
        library(raster)
        library(EBImage)
        library(tools)
        nir <- filter(data, month == m22[1] )
        lapply(fami, function(f){
          f2 <-list(unlist(f))
          lapply(f2, function(f22){
            nir2 <- filter(nir, Fam == f22[1]  )
            nir3 <- nir2[,-c(1:4)]
            chl <- nir2[,4]
            library(tidyverse)
            matou_vis <- nir3 %>% dplyr:: mutate(
              ndvi=  ((result_NIR - result_Red) / (result_NIR + result_Red)),
              osavi = ((result_NIR-result_Red)*(1+0.16)) / (result_NIR + result_Red + 0.16),
              gndvi = (result_NIR-result_Green)/(result_NIR+result_Green),
              savi = ((result_NIR - result_Red)*(1+0.5))/((result_NIR + result_Red+0.5)),
              msavi = (2*result_NIR+1-sqrt((2*result_NIR+1)^2-8*(result_NIR-result_Red)))/2,
              gci = result_NIR/result_Green-1,
              RECI = result_NIR/result_RedEdge-1,
              LCI = (result_NIR-result_RedEdge)/(result_NIR+result_Red),
              GRVI =(result_Green-result_Red)/(result_Green+result_Red),
              MGRVI =(result_Green^2-result_Red^2)/(result_Green^2+result_Red^2 ),
              RGBVI =(result_Green^2-result_Red*result_Blue)/(result_Green^2+result_Red*result_Blue),
              NDRE= (result_NIR-result_RedEdge)/(result_NIR+result_RedEdge),
              MACI= result_NIR/result_Green,
              ARI= result_Green/result_NIR,
              MARI=(result_Green^(-1)-result_RedEdge^(-1))/result_NIR
              
            )
            # library(raster)
            # # # # set up an 'empty' raster, here via an extent object derived from your data
            # xy <- as.matrix(nir3[,1:2])
            tryCatch({
              library(raster)
              nir2 <- sapply(matou_vis[,-c(1:2,8:9)], function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))
              
              matou_vis2 <- cbind.data.frame(matou_vis[, c(1:2)], nir2)  
              
              # create spatial points data frame
              spg <- matou_vis2
              coordinates(spg) <- ~ x + y
              # coerce to SpatialPixelsDataFrame
              gridded(spg) <- TRUE
              # coerce to raster
              rasterDF <- stack(spg)
              rasterDF2 <-scale(rasterDF)
              # plot(rasterDF)
              library(RStoolbox)
              library(rasterVis)
              library(viridis)
              # myPal <- RColorBrewer::brewer.pal('Blues', n=9)
              # myTheme <- rasterTheme(region = viridis(1500))
              
              # plotRGB(rasterDF, stretch = "lin",axes=T,main= paste0(m22[1],"_", f22[1] ),
              #         r =1, g =5, b =3)
              df <-  ggRGB(rasterDF2, r=1, g=5, b=3, stretch = 'lin')+ggtitle(paste0(m22[1],"_", f22[1] ))
              print(df)
              imgg <- as.array(rasterDF2 )
              ds3 <- EBImage::as.Image(imgg)
              ds3 <- ds3[,,-c(1:2)]
              ds4 <- resize(ds3,30,30)
              f_name  <- list(ds4)
              names(f_name) <- paste0(m22[1],'_',f22[1])
              f_namechl  <- list(chl)
              names(f_namechl) <- paste0(m22[1],'_chl_',f22[1])
              message( paste0(m22[1],'_',f22[1]))
              return(list(f_namechl,f_name))
            },error = function(e) {NULL})
            # e <- extent(xy)
            
            # # # # you need to provide a function 'fun' for when there are multiple points per cell
            
            
          })})},error = function(e) {
            message('Caught an error!')
            cat("ERROR :", conditionMessage(e), "\n")
            print(e)})})})}



month <- (unique(data$month))
fam <- (unique(data$Fam))
library(DT)







ui <- navbarPage( h2("Forestry Phenomics"),

           
                  tabPanel(h4("Segmentation"),
                           fluidPage(
                             
                             theme = shinytheme("yeti"),
                             # Application title
                             titlePanel("Individual tree segmentation"),
                             # # Input: Select the random distribution type ----
                             # radioButtons("dist", "Distribution type:",
                             #              c("Normal" = "norm")),
                             
                           
                             # Sidebar with a slider input for number of bins 
                             sidebarLayout(
                               sidebarPanel(width = 4,
                                            # Input: Select the random distribution type ----
                                            
                                            
                                            selectInput("status", label = h3("fam data"),
                                                        choices = fam,
                                                        selected = fam[1], multiple = F),
                                            selectInput("status2", label =h3("Months data"),
                                                        choices = month,
                                                        selected = month[1], multiple = F),
                                          
                                            
                                            fluidRow(
                                               numericInput("width_png","Width of PNG", value = 1600) ,
                                               numericInput("height_png","Height of PNG", value = 1200 ),
                                               numericInput("resolution_PNG","Resolution of PNG", value = 144 ),
                                               style = "margin-top: 25px;",
                                                     downloadButton('downloadPlotPNG','Download Layer plot PNG'),
                                               downloadButton('downloadPlotPNG2','Download singletree PNG') 
                                            ) 
                                            
                                            
                               ),
                               
                               # Show a plot of the generated distribution
                               mainPanel(
                                 
                                 plotOutput('distPlot',width = "100%", height = "800px"),
                                
                                 # Output: Tabset w/ plot, summary, and table ----
                                 tabsetPanel(type = "tabs",
                                             tabPanel("image info", verbatimTextOutput("summary")),
                                             tabPanel("Data", verbatimTextOutput("summary2")),
                                              tabPanel("Layer plot",plotOutput('predictPlot',width = "100%", height = "1000px"))   
                                             
                                             
                                 )
                               )
                             )
                           )
                  ),
                  
                  tabPanel( h4("VIs graphs generation"),
                           fluidPage(
                             
                             theme = shinytheme("yeti"),
                             # Application title
                             titlePanel("Plot VIs"),
                           
                             sidebarLayout(
                               sidebarPanel(width = 4,
                                          
                                            column(
                                              4,
                                              fileInput(
                                                "red1",
                                                "Please upload an red images ",
                                                multiple = T,
                                                accept = c("tif",
                                                           ".tif")
                                              )
                                            ),
                                            column(
                                              4,
                                              fileInput(
                                                "gree1",
                                                "Please upload a green images ",
                                                multiple = T,
                                                accept = c("tif",
                                                           ".tif")
                                              )
                                            ),
                                            column(
                                              4,
                                              fileInput(
                                                "blue1",
                                                "Please upload a blue images ",
                                                multiple = T,
                                                accept = c("tif",
                                                           ".tif")
                                              )
                                            ),
                                            column(
                                              5,
                                              fileInput(
                                                "redege1",
                                                "Please upload an redege images ",
                                                multiple = T,
                                                accept = c("tif",
                                                           ".tif")
                                              )
                                            ),
                                            column(
                                              5,
                                              fileInput(
                                                "NIR1",
                                                "Please upload a NIR images ",
                                                multiple = T,
                                                accept = c("tif",
                                                           ".tif")
                                              )
                                            ),
                                            
                                            column(12, wellPanel(style = "background: white",
                                              
                                              shinyWidgets::radioGroupButtons(
                                                "testh2o",
                                                h3("Plot Map:"),
                                                direction = "horizontal",
                                                individual = TRUE,
                                                width =
                                                  '100%',
                                                justified =F,
                                                
                                                c(
                                                  "Red" = "Red",
                                                  "Green" = "Green",
                                                  "Blue" = "Blue",
                                                  "Rededage" = "Rededage",
                                                  'NIR' = 'NIR',
                                                  "ndvi" = "ndvi",
                                                  "osavi" = "osavi",
                                                  "gndvi" = "gndvi",
                                                  "savi" = "savi",
                                                  "msavi" = "msavi",
                                                  "gci" = "gci",
                                                  "RECI" = "RECI",
                                                  "LCI" = "LCI",
                                                  "GRVI" = "GRVI",
                                                  "MGRVI" = "MGRVI",
                                                  "NDRE" = "NDRE",
                                                  "MACI" = "MACI",
                                                  "ARI" = "ARI",
                                                  "MARI" = "MARI"
                                                ),
                                                selected ="ndvi",
                                                checkIcon = list(yes = icon("ok",
                                                                            lib = "glyphicon"))
                                              )
                                              
                                            )),
                                            fluidRow(
                                              column(4, numericInput("width_png2","Width of PNG", value = 1600)) ,
                                              column(4, numericInput("height_png2","Height of PNG", value = 1200 )),
                                              column(4, numericInput("resolution_PNG2","Resolution of PNG", value = 144 )),
                                              column(4, numericInput("width_pdf","Width of pdf", value = 16)) ,
                                              column(4, numericInput("height_pdf","Height of pdf", value = 12 )),
                                              
                                              style = "margin-top: 25px;",
                                              downloadButton('downloadPlotPNG11','Download single layer plot PNG'),
                                              downloadButton('downloadPlotPNG22','Download RGB plot PDF')
                                            )
                                            
                                            
                               ),
                               
                               # Show a plot of the generated distribution
                               mainPanel( splitLayout(
                                                        style = "border: 1px solid silver:",
                                                        # cellWidths = c(500, 500),
                                                        plotOutput("plotgraph1" , width = "100%", height = "800px"),
                                                        plotOutput("plotgraph2", width = "100%", height = "800px")
                                                        # plotOutput("plotgraph3 )
                                                      )       
                                  
                               )
                             )
                           )
                  ),
                  
                  tabPanel(h4( "Tree identify and spetral extraction"),
                             fluidPage(
                             
                             theme = shinytheme("yeti"),
                             # Application title
                             titlePanel("Segmentation"),
                             sidebarLayout(
                               sidebarPanel(width = 4,
                               # textInput("caption", "Multi Spec", 'e:/shiny/rasterimage/'),
                               # textInput("caption2", "LAS data ", 'E:/shiny/'),
                                
                               # Input: Select a file ----
                              
                                 column(6,fileInput("file1", "Choose las File",
                                         multiple = FALSE,
                                         accept = c("las",
                                                    ".las"))),

                                 column(6,fileInput("file2", "Choose raster File",
                                         multiple = T,
                                         accept = c("tif",
                                                    ".tif"))) ,
                             
                               
                                 fluidRow( 
                                   column(4,downloadButton("tif_data", label = "Download spectral data")),
                                   column(4, downloadButton("cleaned_data", label = "Download las data"))) ,
                               # Horizontal line ----
                               tags$hr(),
                               
                              wellPanel(fluidRow(
                                 # column(12, sliderInput("select2",
                                 #             "TreeID number:",
                                 #             min = 1,
                                 #             max = 1650,
                                 #             step =1,
                                 #             value = 3)), 
                                 column(12,sliderInput(inputId = 'heightdata',
                                           label = 'height limits:',
                                           value = 0,
                                           step =0.1,
                                           min =0,
                                           max =10))
                                 # column(12,  sliderInput(inputId = 'wscontro',
                                 #                      label = 'find tree ws control:',
                                 #                      value = 6,
                                 #                      step =0.1,
                                 #                      min =0 ,
                                 #                      max =10)) 
                               
                               
                             )),
                             fluidRow(
                               column(4, numericInput("width_png3","Width of PNG", value = 1600)) ,
                               column(4, numericInput("height_png3","Height of PNG", value = 1200 )),
                               column(4, numericInput("resolution_PNG3","Resolution of PNG", value = 144 )),
                               column(4, numericInput("width_pdf3","Width of pdf", value = 16)) ,
                               column(4, numericInput("height_pdf3","Height of pdf", value = 12 )),
                               
                               style = "margin-top: 25px;"
                              
                             )
                             ) ,
                             
                             
                             
                             mainPanel(
                               h1("Introducing Shiny"),
                               p("Shiny is a new package from RStudio that makes it ", 
                                 em("incredibly easy "), 
                                 "to build interactive web applications with R."),
                               fluidRow( column(3,   sliderInput(inputId = 'wscontro',
                                                                    label = 'find tree ws control:',
                                                                    value = 6,
                                                                    step =0.1,
                                                                    min =0 ,
                                                                    max =10)),
                                         column(3, downloadButton('downloadrgball',
                                                                  'Download plot with selected treeID',class = "butt1"),
                                                ) , 
                                         column(3,    sliderInput("select2",
                                                                "TreeID number:",
                                                                min = 1,
                                                                max = 1650,
                                                                step =1,
                                                                value = 3)),
                                         column(3, downloadButton('downloadsfall','Download SF with selected treeID') ) 
                                         # column(4,  downloadButton('downloadsfall','Download SF with selected treeID'))
                                         ) ,
                               fluidRow( column(12, splitLayout(
                                 style = "border: 1px solid silver:",
                                 # cellWidths = c(500, 500),   column(12,  sliderInput(inputId = 'wscontro',
                                 
                                
                                 plotOutput("predictPlo" , width = "100%", height = "800px"),
                               
                                 plotOutput("predictPlot5", width = "100%", height = "800px")
                                 # plotOutput("plotgraph3 )
                               ))),
                               
                               # plotOutput('distPlo',width = "100%", height = "800px"),
                               
                               tableOutput("contents2"),
                               # tabPanel("distPlot2",plotOutput('distPlot2',width = "100%", height = "1000px")),
                               # Output: Tabset w/ plot, summary, and table ----
                               tabsetPanel(type = "tabs",
                                           tabPanel("Point cloud information", verbatimTextOutput("summar")),
                                           # tabPanel("final data",  DT::dataTableOutput("mytable")),
                                           tabPanel("LAS polt information",
                                                    
                                                      sidebarPanel(
                                                    fluidRow( 
                                                      actionButton("dodo", "plot las point cloud"), 
                                                    actionButton("dodo1", "plot single tree point cloud")
                                                    )),
                                                    mainPanel(tableOutput("contents"),
                                                              tableOutput("contents22"))),
                                          
                                           tabPanel("Spectral information", verbatimTextOutput("summar2")),
                                           
                                           # tabPanel("Polygon plot",plotOutput('predictPlot5',width = "100%", height = "1000px")),
                                           tabPanel("Magic",plotOutput('distPlo',width = "100%", height = "1000px")),
                                           tabPanel("Layer plot",plotOutput('predictPlot3',width = "100%", height = "1000px")),
                                           tabPanel("Single tree RGB plot",plotOutput('predictPlot4',width = "100%", height = "1000px")),
                                          
                                           tabPanel("Final data output", 
                                                    
                                                    sidebarPanel(
                                                      fluidRow( 
                                                        column(4, downloadButton("tife_data", label = "Downloadc final data"))
                                                      )),
                                                    mainPanel(verbatimTextOutput("su2")))
                                             
                               )
                             )
                           )
                           )
                           
                           ) 

                  
                  
                  
                  
)




library(tools)
# Define server logic required to draw a histogram

options(shiny.maxRequestSize=1000*1024^2)
# Define server logic required to draw a histogram
server <- function(input, output) {

  trend_data <- reactive({
    
    sele <-input$status
  })
  trend_data2 <- reactive({
    
    # `validate()` is additional; to prepare friendly message for error
    selec <-input$status2
  })
  
  
  #Generate a summary of the dataset ----
  output$summary <- renderPrint({
    
    chl_tree <- ras_im_alin(fami= trend_data(),monthi=trend_data2())
    
    gfdgh1 <- chl_tree %>%unlist(recursive = F)%>%unlist(recursive = F)%>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    fffghj<-gfdgh1%>%unlist(recursive = F)%>%unlist(recursive = F)
    library(data.table)
    
    fffghj2 <- (fffghj[!names(fffghj) %like% 'chl']) 
    fffghj2[sapply(fffghj2, is.null)] <- NULL
    fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
    print(fffghj2[[1]])
  }) 
  
  # dat243 <- reactive({
  #   select <- switch(input$select,
  #                    "Sep image data"= 23,
  #                    "Oct image data"=20,
  #                    "Nov image data"= 1,
  #                    "May image data"= 24,
  #                    "June image data"= 28, 
  #                    "July image data"= 29
  #                    
  #   )
  #   
  # })
  # 
  output$distPlot <- renderPlot({
    
    chl_tree <-sigletree()
    
    chl_tree
    
  })
  
  
  output$summary2  <- renderPrint({
    chl_tree <- ras_im_alin(fami= trend_data(),monthi=trend_data2())
    
    gfdgh1 <- chl_tree %>%unlist(recursive = F)%>%unlist(recursive = F)%>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    fffghj<-gfdgh1%>%unlist(recursive = F)%>%unlist(recursive = F)
    library(data.table)
    stOdds <- (fffghj[names(fffghj)%like% 'chl'])%>% invoke(rbind,.)  
    spl <- strsplit(rownames(stOdds),'_|[.]')
    stOdds$name <- sapply(spl, function(x){
      y=(paste0(x[1],'_',x[3],'_',x[4],'_',x[5]))
    })
    stOdds
    
    
    
  })
  output$predictPlot <- renderPlot({
    chl_tree <- ras_im_alin(fami= trend_data(),monthi=trend_data2())
    
    gfdgh1 <- chl_tree %>%unlist(recursive = F)%>%unlist(recursive = F)%>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    fffghj<-gfdgh1%>%unlist(recursive = F)%>%unlist(recursive = F)
    library(data.table)
    
    fffghj2 <- (fffghj[!names(fffghj) %like% 'chl']) 
    fffghj2[sapply(fffghj2, is.null)] <- NULL
    # EBImage::display(fffghj2[[1]],method = 'raster',all=T )
    # EBImage::display(ds4,method = 'raster',all=T)
    fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
    y <- brick(fffghj2[[1]])
    library(viridis)
   plot(y,col=viridis(20))
    
    
  })
  
  # output$distPlo <- renderPlot({
  #   select <-input$caption
  #   library("lidR")
  #   library("rgdal")
  #   library(raster)
  #   library(tidyverse)
  #   dsf1 <- raster_read(select)
  #   dsf1 <- dsf1 %>% unlist(recursive = F)%>%  unlist(recursive = F)
  #   
  #   dew <- dsf1%>% unlist(recursive = F) %>% stack()
  #   dfe <-  ggRGB(dew, r=1, g=3, b=2, stretch = 'lin')
  #   print(dfe)
  #   
  #   # par(mfrow = c(2, 3))
  #   # library(viridis)
  #   # library("quickPlot")
  #   # message('plot raster images')
  #   # lapply(dsf1, function(x)
  #   #   plot(x, col= viridis(20)))
  # })


  
  
  
  
  
  # #  
  # output$summary <- renderPrint({
  # 
  #      req(input$file1)
  # 
  #   tryCatch(
  #     {
  #       library("lidR")
  #       library("rgdal")
  # las_12 <- readLAS(input$file1$datapath,
  #                   filter = "-change_classification_from_to 1 2",
  #                   select = "xyzirc"
  # 
  # )
  #       print(las_12)
  #     },
  #     error = function(e) {
  #       # return a safeError if a parsing error occurs
  #       stop(safeError(e))
  #     }
  #   )
  # 
  # })
  
  
  df_products_upload <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
        return('please upload las cloud data')
    las_12 <- lapply(inFile$datapath,function(m){
      fdd <- lidR::readLAS(m, filter = "-change_classification_from_to 1 2",
                    select = "xyzirc" )
      las <- lidR::add_lasattribute(fdd, 1,'month', "A new data")
      
    } ) 
    
    las_12
 
  })
 
 
  
  testedd   <- reactive({
    chl_tree <- ras_im_alin(fami= trend_data(),monthi=trend_data2())
    
    gfdgh1 <- chl_tree %>%unlist(recursive = F)%>%unlist(recursive = F)%>%unlist(recursive = F)
    gfdgh1[sapply(gfdgh1, is.null)] <- NULL
    fffghj<-gfdgh1%>%unlist(recursive = F)%>%unlist(recursive = F)
    library(data.table)
    
    fffghj2 <- (fffghj[!names(fffghj) %like% 'chl']) 
    fffghj2[sapply(fffghj2, is.null)] <- NULL
    # EBImage::display(fffghj2[[1]],method = 'raster',all=T )
    # EBImage::display(ds4,method = 'raster',all=T)
    fffghj2[[1]][is.na(fffghj2[[1]])] <- 0
    y <- brick(fffghj2[[1]])
    library(viridis)
    plot(y,col=viridis(20))
    
    
  }) 
  
   
  output$downloadPlotPNG <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ggVolcanoR_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png, height = input$height_png, res = input$resolution_PNG)
      testedd()
      dev.off()
      },
    
    contentType = "application/png" # MIME type of the image
    
  )
  
  sigletree <- reactive({
    
     ras_im_alin(fami= trend_data(),monthi=trend_data2())
    
    
  })
  
  output$downloadPlotPNG2 <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("ggVolcanoR_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      
      png(file, width = input$width_png, height = input$height_png, res = input$resolution_PNG)
      ras_im_alin(fami= trend_data(),monthi=trend_data2())
      dev.off()
    },
    
    contentType = "application/png" # MIME type of the image
    
  )
  
  
   
  # df_raster_upload <- reactive({
  #   inFile1 <- input$file2
  #   if (is.null(inFile1))
  #     return(NULL)
  #  
  #   library("lidR")
  #   library("rgdal")
  #   library(raster)
  #   library(tidyverse)
  #   dsf1 <- raster(inFile1$datapath)
  #    
  #   dsf1
  # }) 
  
  getData <- reactive({
    library(data.table)
    inFile1 <- input$file2
    if (is.null(inFile1))
      return(print('please upload raster images'))
    plr <- lapply(inFile1$datapath, function(m){
         rs <- raster(m )
         rs
        })
    
  for (i in 1:length(plr)) {
    names(plr[[i]]) <- inFile1$name[[i]]
  }
 
   plr
  
  })
  
  
  # output$contents <- renderPrint( {
  #   gg <- getData() 
  #   print(gg) 
  #   
  # }) 
  # 

output$cleaned_data <- downloadHandler(
    filename = "mydata.rds",
    content = function(file) {readr::write_rds(df_products_upload(), file)}
  )
  

  data_las <- reactive({
    select <-input$caption2
    imag <- list.files(
      path = select,
      pattern = '*.las',
      all.files = T,
      full.names = T,
      no.. = T
    )
    library("lidR")
    library("rgdal")
    library(raster)
    library(tidyverse)
    las_12 <-lidR:: readLAS(imag[1],
                            filter = "-change_classification_from_to 1 2",
                            select = "xyzirc" )
  })
  
  output$summar  <- renderPrint({
    las22 <- df_products_upload()
    print(las22)
  })
  


  data_dsf1 <- reactive({
    select1 <-input$caption
    library("lidR")
    library("rgdal")
    library(raster)
    library(tidyverse)
    dsf1 <- raster_read(select1)
    dsf1 <- dsf1 %>% unlist(recursive = F) %>%  unlist(recursive = F)
    
  })
    
  # df_raster_upload <- reactive({
  #   inFile1 <- input$file2
  #   if (is.null(inFile1))
  #     return(NULL)
  #  
  #   library("lidR")
  #   library("rgdal")
  #   library(raster)
  #   library(tidyverse)
  #   dsf1 <- raster(inFile1$datapath)
  #    
  #   dsf1
  # }) 
  
  #  output$tif_data <- downloadHandler(
  #   filename = "raster.rds",
  #   content = function(file) {readr::write_rds(df_raster_upload(), file)}
  # )
  #  

  
  
  
  output$summar2  <- renderPrint({
    dsf1 <-  getData()#df_raster_upload()
    if (is.null(dsf1))
      return(NULL)
    print(dsf1)
  })
  # 
  
  output$distPlo <- renderPlot({
    select23 <-  getData()
    if (is.null(select23))
      return(NULL)
    plot( stack(select23),col= viridis(200))

    # par(mfrow = c(2, 3))
    # library(viridis)
    # library("quickPlot")
    # message('plot raster images')
    # lapply(dsf1, function(x)
    #   plot(x, col= viridis(20)))
  })

  
  data_ext2  <- reactive({
    dsf1 <- getData() #data_dsf1() 
    las_12 <- df_products_upload() #data_las()
    if (is.null(dsf1))
      return('please upload raster images')
    if (is.null(las_12))
      return('please upload las cloud data')
    # las_12 <- add_lasattribute(las_12, 0, "month", "A new data")
    las_list <- las_12
    names(las_list) <- c('month')
    dsf_list <- list(dsf1)
    names(dsf_list) <- c('month')
    par(mfrow = c(2,3))
    data_all <- multi_rasl (las_list,dsf_list,kwsindice = input$wscontro )
   
  })
 

mult  <-  reactive({
  las_12 <- df_products_upload() #data_las()
  las_list <- las_12
  names(las_list) <- c('month')
  
  lasdata <-   lapply(las_list, function(x){
      expr <- tryCatch({
        library("lidR")
        library("rgdal")
        library(sfheaders)
        library(stars)
        library(raster)
        library(tidyverse)
        library(sf)
        library(data.table)
        message( paste0(names(x@header@VLR$Extra_Bytes$`Extra Bytes Description`)))
        las = classify_ground(x,csf(cloth_resolution = 0.5, class_threshold = 0.15, rigidness = 1), last_returns = FALSE)
        dtm <- grid_terrain(las, res = 0.5, algorithm = knnidw(k=5,p = 0.5), use_class = c(1L, 2L),keep_lowest = F)
        nlas_dtm  <- las - dtm
        chm    <- grid_canopy(nlas_dtm , res =0.5, p2r())
        ttops  <- find_trees(nlas_dtm , lmf(ws=input$wscontro, hmin=2.6, shape = "circular"))
        algo   <- dalponte2016(chm, ttops )
        lass   <- segment_trees(nlas_dtm, algo, uniqueness ='incremental')
 
      },error = function(e) {
        message('Caught an error!')
        cat("ERROR :", conditionMessage(e), "\n")
        print(e)},
      print(paste("processing Loop", names(las_list), sep = "_"))) 
    })%>% 
    # map_dfr(cbind) %>% 
    unlist(recursive = F)
  lasdata
  })

crown_pol  <-  reactive({
  crown_polo  <- mult()
  # crown_polo  <- unlist(crown_polo)
  # crown_polo1  <- delineate_crowns(crown_polo, func = .stdtreemetrics)
  crown_polo1  <- lapply(crown_polo, function(x){
    crown_polo  <- delineate_crowns(x, func = .stdtreemetrics)
  })

 
})


sflas  <-  reactive({
  crowtemp  <- mult()
  crown_pow  <- lapply(crowtemp, function(x){
    crowtemp1  <- delineate_crowns(x, func = .stdtreemetrics)
    directions_longlat <- spTransform(crowtemp1,sp:: CRS("+proj=longlat +datum=WGS84 +no_defs"))
    sf  <- st_as_sf(directions_longlat)
  })
  
  crown_pow
 
  
})


sertree  <-  reactive({
  crowte  <- mult()
  crown_se  <- lapply(crowte, function(x){
    tree35_36 <- filter_poi(x, treeID == input$select2)
 
    las.tree35_36 <- filter_poi(tree35_36, Z >input$heightdata)
   
    
  })
  
  crown_se
  
  
})
 


randse <- eventReactive(input$dodo1, {
  sele <- sertree()
  
})  


# output$contents2 <- renderPrint({
#   library("lidR")
#   library("rgdal")
#   library(raster)
#   library(tidyverse)
#   print("plot with RGL device")
#   randese <-  randse()
#   if (is.null(randese))
#     return(NULL)
#   
#   randese1  <- lapply(randese, function(xyr){
#     plot(xyr, size = 8, bg = "white",breaks='quantile' ,  
#          axis = TRUE )
#   })
# })


output$contents  <- renderPrint({
  library("lidR")
  library("rgdal")
  library(raster)
  library(tidyverse)
  print("plot with RGL device")
  randese <-  randse()
  if (is.null(randese))
    return(NULL)
  
  randese1  <- lapply(randese, function(xyr){
    plot(xyr, size = 8, bg = "white",breaks='quantile' ,  
         axis = TRUE )
  })
})






plot2  <- reactive({
  fdff <- crown_pol()
  if (is.null(fdff))
    return(NULL)
  # lapply(fdff,plot) 
  
  dre <- lapply(fdff, function(xx){
    plot(xx)
    text(xx, paste(xx@data$treeID ),
         cex=1,col='blue')
    idnum <- xx[xx@data$treeID == input$select2,]
    plot(idnum,add=T,alpha=0.5)
    text(idnum, paste(idnum@data$treeID ),add=T,
         cex=1,col='red',alpha=0.5)
    # idnum <- xx@data$treeID 
    # lapply(idnum, function(yyy){
    #   
    #   
  })
 
})


plot23  <- reactive({
  fdff <- crown_pol()
  if (is.null(fdff))
    return(NULL)
  # lapply(fdff,plot) 
  
  dre <- lapply(fdff, function(xx){
    plot(xx)
    # text(xx, paste(xx@data$treeID ),
    #      cex=1,col='blue')
    idnum <- xx[xx@data$treeID == input$select2,]
    plot(idnum,add=T,alpha=0.5)
    text(idnum, paste(idnum@data$treeID ),add=T,
         cex=1,col='red',alpha=0.5)
    # idnum <- xx@data$treeID 
    # lapply(idnum, function(yyy){
    #   
    #   
  })
  
})


output$predictPlot5  <- renderPlot({
  fdff2 <-  plot23()
  print(fdff2)
})

output$downloadsfall  <- downloadHandler(
  filename = function() {
    x <- gsub(":", ".", Sys.time())
    paste("spetral_", gsub("/", "-", x), ".png", sep = "")
  },
  content = function(file) {
    
    png(file, width = input$width_png3, height = input$height_png3, res = input$resolution_PNG3)
    print(plot2())
    dev.off()
  },
  
  contentType = "application/png" # MIME type of the image
  
)



rgbplotwithid  <- reactive({
  
  sfff <- sflas()
  select23 <-  getData() 
  if (is.null(sfff))
    return(NULL)
  if (is.null(select23))
    return(NULL)
  se2 <-  stack(select23)
  
  # nlayers(s)
  if(nlayers(se2) < 3){
    plot(se2[[1]] ,col= viridis(200) ) 
    
    lapply(sfff,function(xyxy){
      
      plot(xyxy,add=T,col='white',alpha=0.4)
      
      idnum <- xyxy[xyxy$treeID == input$select2,]
      plot(idnum, add=T,alpha=0.4,col='red')
      text(idnum, paste(idnum$treeID ),
           cex=1.65,col='blue',alpha=0.4) 
      
    }) 
  } else{
    # sfff3 <-  sfff %>%  unlist(recursive = F)%>%  unlist(recursive = F)
    
    # plotRGB(se2, scale=maxValue(se2))
    
    lapply(sfff,function(xyxy){
      library(RStoolbox)
      idnum <- xyxy[xyxy$treeID == input$select2,]
      p <-ggRGB(se2,  stretch = "hist")+ 
        geom_sf(data = xyxy, fill='blue',alpha=0.4   )+
        # geom_sf_label(data = idnum,aes(label = treeID,alpha = 0.7, vjust =1.1)  )
        geom_sf(data = idnum, fill='red'   )+
        ggrepel::geom_label_repel(
          data = idnum,
          aes(label = treeID, geometry = geometry),
          stat = "sf_coordinates",
          min.segment.length = 0,
          colour = "magenta",
          segment.colour = "magenta"
        ) 
        # ggrepel::geom_label_repel(
        #   data = xyxy,
        #   aes(label = treeID, geometry = geometry),
        #   stat = "sf_coordinates",
        #   min.segment.length = 0,
        #   colour = "blue",
        #   segment.colour = "blue"
        # )
      
      print(p)
      
    })
    
  }
  
  
})



output$downloadrgball  <- downloadHandler(
   
  
  filename = function() {
    x <- gsub(":", ".", Sys.time())
    paste("spetral_", gsub("/", "-", x), ".png", sep = "")
  },
  content = function(file) {
    rgbg <- rgbplotwithid()
    png(file, width = input$width_png3, height = input$height_png3, res = input$resolution_PNG3)
    print(rgbg)
    dev.off()
  },
  
  contentType = "application/png" # MIME type of the image
  
)

 

output$predictPlo  <- renderPlot({
  rgbplotwithidw<- rgbplotwithid()
  print(rgbplotwithidw)
 
})





















# output$predictPlo  <- renderPlot({
#  
#   sfff <- sflas()
#   select23 <-  getData() 
#   if (is.null(sfff))
#     return(NULL)
#   if (is.null(select23))
#     return(NULL)
#   se2 <-  stack(select23)
#   
#   # nlayers(s)
#   if(nlayers(se2) < 3){
#     plot(se2[[1]] ,col= viridis(200) ) 
#     
#     lapply(sfff,function(xyxy){
#       
#       # plot(xyxy,add=T)
#       
#       idnum <- xyxy[xyxy$treeID == input$select2,]
#       plot(idnum,col='red',add=T)
#       text(idnum, paste(idnum$treeID ),
#            cex=1.65,col='red') 
#       
#     }) 
#     } else{
#       # sfff3 <-  sfff %>%  unlist(recursive = F)%>%  unlist(recursive = F)
#       
#     # plotRGB(se2, scale=maxValue(se2))
#      
#     lapply(sfff,function(xyxy){
#       library(RStoolbox)
#       idnum <- xyxy[xyxy$treeID == input$select2,]
#      p <-ggRGB(se2,  stretch = "hist")+ geom_sf(data = idnum,fill = "red" )+
#        # geom_sf_label(data = idnum,aes(label = treeID,alpha = 0.7, vjust =1.1)  )
#        
#        ggrepel::geom_label_repel(
#          data = idnum,
#          aes(label = treeID, geometry = geometry),
#          stat = "sf_coordinates",
#          min.segment.length = 0,
#          colour = "magenta",
#          segment.colour = "magenta"
#        )
#  
#   print(p)
#   
#     })
#   
#     }
#   
#  
# })

 

# plot(fdff,main= '70m')
# fdff2 <- fdff %>%  unlist(recursive = F)
# for (i in c(1:length(fdff2@data$treeID))) {
#   crownsPoly1 <- fdff2[fdff2@data$treeID== i,]
#   plot(crownsPoly1, border = "blue", lwd = 0.5,add=T)
#   text(crownsPoly1, paste(i),
#        cex=0.65,col='red') 
# }
# crown_polo1  <- lapply(fdff, function(x){
#   
#   plot(x,main= '70m')
#   for (i in c(1:length(x@data$treeID))) {
#     crownsPoly1 <- x[x@data$treeID== i,]
#     plot(crownsPoly1, border = "blue", lwd = 0.5,add=T)
#     text(crownsPoly1, paste(i),
#          cex=0.65,col='red') 
#   }
#   
#   
# })


 




  
  output$tife_data <- downloadHandler(
    filename = "finaloutput data.rds",
    content = function(file) {readr::write_rds(data_ext2(), file)}
  )
  
 
  # output$mytable  <- DT::renderDataTable({
  #   data_all1 <- data_a()
  #  
  #   head(data_all1)
  #   # data_Jan%>% filter(treeID == '98')%>%
  #   #   ggplot(aes(x,y,fill= result_RedEdge))+geom_tile()
  #   
  # })

  
 finaldata  <- reactive({
    dsf1 <-  data_ext2()
    if (is.null(dsf1))
      return(NULL)
    data_Jan <- dsf1 %>%invoke(cbind,.)
    names(data_Jan) <-  gsub('month|[.]|tif','',names(data_Jan))
    data_Jan <- data_Jan %>% dplyr::select(x,y,treeID,chm, everything())
    data_Jan
  })
  
 
 output$su2  <- renderPrint({
   dsf331 <- finaldata()
   print(dsf331)
 })      
   
              
 
  
 output$predictPlot3  <- renderPlot({
   library(tidyverse)
   library(raster)
   library(EBImage)
   library(tools)
   nir <- filter(finaldata(), treeID == input$select2 )
   nir3 <- nir[,-c(1:4)]
  library(raster)
     nir2 <- sapply(nir3, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))
     
     matou_vis2 <- cbind.data.frame(nir[,c(1:2,4)], nir2)  
     matou_vis2 <- matou_vis2 %>% filter(chm > input$heightdata)
     # create spatial points data frame
     spg <- matou_vis2
     coordinates(spg) <- ~ x + y
     # coerce to SpatialPixelsDataFrame
     gridded(spg) <- TRUE
     # coerce to raster
     rasterDF <- stack(spg)
     # rasterDF2 <-scale(rasterDF)
     # plot(rasterDF)
     library(RStoolbox)
     library(rasterVis)
     library(viridis)
     plot(rasterDF ,col= viridis(200))
     
 })
 
  
 output$predictPlot4  <- renderPlot({
   library(tidyverse)
   library(raster)
   library(EBImage)
   library(tools)
   nir <- filter(finaldata(), treeID == input$select2 )
   nir3 <- nir[,-c(1:4)]
   library(raster)
   nir2 <- sapply(nir3, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm=T)))
   
   matou_vis2 <- cbind.data.frame(nir[,c(1:2,4)], nir2)  
   matou_vis2 <- matou_vis2 %>% filter(chm > input$heightdata)
   # create spatial points data frame
   spg <- matou_vis2
   coordinates(spg) <- ~ x + y
   # coerce to SpatialPixelsDataFrame
   gridded(spg) <- TRUE
   # coerce to raster
   rasterDF <- stack(spg)
   # rasterDF2 <-scale(rasterDF)
   # plot(rasterDF)
   library(RStoolbox)
   library(rasterVis)
   library(viridis)
   
   if (nlayers(rasterDF) < 3) {
     print("at least 3 layers needed for RGB plot")
     
 
   } else{
   
   if (nlayers(rasterDF)>3|nlayers(rasterDF)< 5) {
     df <-  ggRGB(rasterDF,  2, 3, 4, 
                  
                  stretch = 'lin') + ggtitle(paste0('tree ID-', input$select2))
     print(df)
 
   }else{
     
     df <-  ggRGB(rasterDF, 5, 3, 2,
                  
                  stretch = 'lin') + ggtitle(paste0('tree ID-', input$select2))
     print(df)
     
     
   }
   }  


 })


 
 randomVals <- eventReactive(input$dodo, {
   sele <-df_products_upload()
    
 })  
 
 
 output$contents22 <- renderPrint({
   library("lidR")
   library("rgdal")
   library(raster)
   library(tidyverse)
   print("plot with RGL device")
  lapply(randomVals(), plot ) 
 })
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
   
  
  
  
# 
  red2 <-  reactive({
    library(raster)
    sele1 <- input$red1
    dsf1 <- raster::raster(sele1$datapath)
  })
  green2 <-  reactive({
    library(raster)
    sele2 <- input$gree1
    dsf2 <- raster::raster(sele2$datapath)
  })
  blue2 <-  reactive({
    library(raster)
    sele3 <- input$blue1
    dsf3 <- raster::raster(sele3$datapath)
  })
  redege2 <-  reactive({
    library(raster)
    sele4 <- input$redege1
    dsf4 <- raster::raster(sele4$datapath)
  })
  NIR2 <-  reactive({
    library(raster)
    sele5 <- input$NIR1
    dsf5 <- raster::raster(sele5$datapath)
  })

#   # randomVals <- eventReactive(input$go, {
#   #
#   #   st1 <- raster::stack(red2(),green2(),blue2(),redege2(),NIR2())
#   #   st1
#   #
#   # })
# 
  ndviind <- reactive({
    st1 <-  ((NIR2() - red2()) / (NIR2() + red2()))
  })
  osavi <- reactive({
    osavi = ((NIR2() - red2()) * (1 + 0.16)) / (NIR2() + red2() + 0.16)
  })
  gndvi <- reactive({
    gndvi = (NIR2() - green2()) / (NIR2() + green2())
  })
  savi <- reactive({
    savi = ((NIR2() - red2()) * (1 + 0.5)) / ((NIR2() + red2() + 0.5))
  })
  msavi <- reactive({
    msavi = (2 * NIR2() + 1 - sqrt((2 * NIR2() + 1) ^ 2 - 8 * (NIR2() - red2()))) /
      2
  })
  gci <- reactive({
    gci = NIR2() / green2() - 1
  })
  RECI <- reactive({
    RECI = NIR2() / redege2() - 1
  })
  LCI <- reactive({
    LCI = (NIR2() - redege2()) / (NIR2() + red2())
  })
  GRVI <- reactive({
    GRVI = (green2() - red2()) / (green2() + red2())
  })
  MGRVI <- reactive({
    MGRVI = (green2() ^ 2 - red2() ^ 2) / (green2() ^ 2 + red2() ^ 2)
  })
  RGBVI <- reactive({
    RGBVI = (green2() ^ 2 - red2() * blue2()) / (green2() ^ 2 + red2() * blue2())
  })

  NDRE <- reactive({
    NDRE = (NIR2() - redege2()) / (NIR2() + redege2())
  })
  MACI <- reactive({
    MACI = NIR2() / green2()
  })

  ARI <- reactive({
    ARI = green2() / NIR2()
  })

  MARI <- reactive({
    MARI = (green2() ^ (-1) - redege2() ^ (-1)) / NIR2()
  })

  dat243 <- reactive({
    select <- switch(
      input$testh2o,
      "Red" = red2(),
      "Green" = green2(),
      "Blue" = redege2(),
      "Rededage" = blue2(),
      "NIR" = NIR2() ,
      "ndvi" = ndviind(),
      "osavi" = osavi(),
      "gndvi" = gndvi(),
      "savi" = savi(),
      "msavi" = msavi(),
      "gci" = gci(),
      "RECI" = RECI(),
      "LCI" = LCI(),
      "GRVI" = GRVI(),
      "MGRVI" = MGRVI(),
      "NDRE" = NDRE(),
      "MACI" = MACI(),
      "ARI" = ARI(),
      "MARI" = MARI()
    )
  })
# 
  pyt <- reactive({
    st1 <- raster::brick(  blue2(), green2(),red2(), redege2(), NIR2())

  })
  output$plotgraph1 <- renderPlot({

    library(RStoolbox)
    library(raster)
    pyt2 <-  RStoolbox::ggRGB(pyt(),
                              #  r=as.character (input$clusters1),
                              # b= as.character( input$clusters2),
                              # g= as.character (input$clusters3),
                              stretch  = 'hist')
    print(pyt2)

  })
# 
  output$plotgraph2 <- renderPlot({
    library(RStoolbox)
    library(raster)
 
    plot(dat243())
 
  })
# 
  output$downloadPlotPNG11  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {

      png(file, width = input$width_png2, height = input$height_png2, res = input$resolution_PNG2)
      plot(dat243())
      dev.off()
    },

    contentType = "application/png" # MIME type of the image

  )

  output$downloadPlotPNG22  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("spetral_", gsub("/", "-", x), ".pdf", sep = "")
    },

    content = function(file) {

      pdf(file, width  = input$width_pdf, height  = input$height_pdf )
      pf<-  ggRGB(pyt(), stretch  = 'hist')
      print(pf)
      dev.off()
    },

    contentType = "application/pdf" # MIME type of the image

  )

  
}
 


library(shiny)
# Run the application 
app <- shinyApp(ui = ui, server = server)
runApp(app, launch.browser = TRUE)
