library(shiny)
library(plotly)
library(tidyverse)
library(shinythemes)
library(leaflet)
library(leaflet.extras)
library(dplyr)
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(scales)
library(readr)
library(maps)
library(mapproj)

# the ui object has all the information for the user-interface
ui <-(
  navbarPage("Viral SpaceTime", id="main",
            
             tabPanel("Interactive Map", 
                     sidebarLayout(
                       sidebarPanel(
                         #absolutePanel(top = 20, left = 30,checkboxInput("phy", "Derivatives of NextStrain", TRUE)),
                         #absolutePanel(top = 50, left = 30,checkboxInput("heat", "Temp Anomoly", FALSE)),
                         #absolutePanel(top = 80, left = 30,checkboxInput("popdens", "Population Density", FALSE)),
                         #absolutePanel(top = 110, left = 30,checkboxInput("pmlvl", "PM2.5 Levels", FALSE)),
                         selectInput("selectedMapType", 
                                     label = "Choose a dataset",
                                     choices = c("Temperature Anomoly" = "heat",
                                                 "Population Density" = "popdens", 
                                                 "PM2.5 Levels" = "pmlvl"),
                                     selected = "Derivatives of NextStrain"),
                         checkboxInput("phy", "NextStrain", TRUE),
                         width=2
                         
                       ),
                       mainPanel(
                         leafletOutput("mymap", height=700),
                         absolutePanel(top = 20, left = 90,  actionButton("resetMap", "Reset")),
                         width=10
                       )
                     )
             ),
             tabPanel("Phylogenies Data", DT::dataTableOutput("phytable")),
             tabPanel("Test",  sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 # Input: Select a file ----
                 fileInput("file1", "Choose CSV File",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 # Input: Select separator ----
                 radioButtons("sep", "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 
                 # Input: Select quotes ----
                 radioButtons("quote", "Quote",
                              choices = c(None = "",
                                          "Double Quote" = '"',
                                          "Single Quote" = "'"),
                              selected = '"'),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select number of rows to display ----
                 radioButtons("disp", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "head")
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(
                 
                 # Output: Data file ----
                 tableOutput("contents")
                 
               ))),
             
             tabPanel("About",includeMarkdown("app_readme.md"))
  ))


server <- function(input, output, session) {
  #import data
  data <- read.csv("datasets/curated/NOAAGlobalTemp/testdat.csv")
  phydata <- read.csv("outputs/quick_alr/latlon_north_america_phy_divergence_rate.csv")
  
  dsn = "datasets/curated/census-app/"
  us.map <- readOGR(dsn = dsn, layer = "cb_2018_us_county_500k", stringsAsFactors = FALSE)
  ## Remove Alaska(2), Hawaii(15), Puerto Rico (72), Guam (66), Virgin Islands (78), American Samoa (60)
  ##  Mariana Islands (69), Micronesia (64), Marshall Islands (68), Palau (70), Minor Islands (74)
  #us.map <- us.map[!us.map$STATEFP %in% c("02", "15", "72", "66", "78", "60", "69",
  #                                        "64", "68", "70", "74"),]
  
  ave_popdens <- read_csv("datasets/curated/census-app/Average_Household_Size_and_Population_Density_-_County.csv")
  
  # Merge spatial df with downloade ddata.
  
  leafmap <- merge(us.map, ave_popdens, by=c("GEOID")) 
  
  # Format popup data for leaflet map.
  popup_dat <- paste0("<strong>County: </strong>", 
                      leafmap$NAME, 
                      "<br><strong>Value: </strong>", 
                      leafmap$B01001_calc_PopDensity)
  
  
  pal1 <- colorQuantile("YlOrRd", NULL, n = 9)
  
  
  pm <- read_csv("datasets/curated/CDC_PM2.5_Concentrations/ave-county-level-pm2.5-2014.csv")
  #sub[subpm$state Daily_County_Level_PM2_5_Concentrations_2001_2014$statefips/]
  
  #pmlvl <- merge(us.map, pm, by.x = "GEOID", by.y = "FIPs")
  
  #define the legend for temp anomoly
  pal <- colorNumeric(
    palette = c('#0B00FF', '#75FB4C', '#FDFE01', '#EA4426'),
    domain = data$mag)
  
  output$phytable = DT::renderDataTable({
    DT::datatable(phydata, options = list(lengthMenu = c(5, 15, 25, 50, 100), pageLength = 15))
  })
  
  output$mymap <- renderLeaflet({
    leaflet(data) %>% 
      setView(lng = -99, lat = 32, zoom = 4)  %>%
      addProviderTiles("Esri.WorldImagery", group="Satellite Map") %>%
      addProviderTiles("CartoDB.DarkMatter", group="Dark Map") %>%
      addTiles(options = providerTileOptions(noWrap = FALSE), group="Street Map") %>%
      addLayersControl(baseGroups = c("Dark Map","Satellite Map","Street Map"), options = layersControlOptions(collapsed = TRUE))
  })
  
  observe({
    leafmap <- leafletProxy("mymap", data = leafmap)
    heatmap <- leafletProxy("mymap", data = data)
    phymap <- leafletProxy("mymap", data = phydata)
    
    if (input$selectedMapType == "popdens") {
      leafmap %>% clearControls() %>% clearHeatmap() %>%
        addPolygons(fillColor = ~pal1(B01001_calc_PopDensity), 
                    fillOpacity = 0.8, 
                    color = "#BDBDC3", 
                    weight = 1,
                    group = "population",
                    popup = popup_dat) %>%
        leaflet::addLegend("bottomright", pal = pal1, values = ~B01001_calc_PopDensity)
    } else  if (input$selectedMapType == "heat") {
      heatmap %>% clearControls() %>% clearShapes()
      heatmap %>% addHeatmap(lng=~LON-270, lat=~LAT, intensity = ~mag, blur =  65, max = 10, radius = 50) %>%
        leaflet::addLegend("bottomright", pal = pal, values = ~mag)
    }  else {
      phymap %>% clearShapes()
    }
    
    if (input$phy) {
      phymap %>% addCircles(data = phydata, lat = ~lat, lng = ~lon, weight = 3, radius = ~sqrt(rate)*1000, popup = ~as.character(loc), label = ~as.character(paste0("Strain: ", sep = " ", loc)), color = "white", fillOpacity = 0.5)
    }
  })

  
  observeEvent(input$resetMap, {
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% setView(lng = -99, lat = 32, zoom = 4) 
  })
}


shinyApp(ui, server)