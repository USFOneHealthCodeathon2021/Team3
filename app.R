library(shiny)
library(plotly)
library(tidyverse)
library(shinythemes)
library(leaflet)
library(leaflet.extras)
library(dplyr)

# the ui object has all the information for the user-interface
ui <-(
  navbarPage("Viral Space Time", id="main",
             tabPanel("Interactive Map", fluidPage(
               fluidRow(
                 column(2, 
                        absolutePanel(top = 20, left = 30,checkboxInput("phy", "Derivatives of NextStrain", TRUE)),
                        absolutePanel(top = 50, left = 30,checkboxInput("heat", "Temp Anomoly", FALSE))
                 ),
                 column(10,
                        leafletOutput("mymap", height=700),
                        absolutePanel(top = 20, left = 70,  actionButton("resetMap", "Reset")),
                 )
               )
             )),
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
  phydata <- read.csv("outputs/quick_alr/latlon_north_america_phy.csv")
  
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
    proxy <- leafletProxy("mymap", data = phydata)
    proxy %>% clearMarkers()
    if (input$phy) {
      proxy %>% addCircles(data = phydata, lat = ~lat, lng = ~lon, weight = 1, radius = 3, popup = ~as.character(loc), label = ~as.character(paste0("Strain: ", sep = " ", loc)), color = "white", fillOpacity = 0.5)
    }
    else {
      proxy %>% clearShapes()
    }
  })
  
  observe({
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% clearMarkers()  %>% clearControls()
    if (input$heat) {
      proxy %>% addHeatmap(lng=~LON-270, lat=~LAT, intensity = ~mag, blur =  70, max = 10, radius = 40) %>%
        leaflet::addLegend("bottomright", pal = pal, values = ~mag)
    }
    else{
      proxy %>% clearHeatmap() %>% clearControls()
    }
  })
  
  observeEvent(input$resetMap, {
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% setView(lng = -99, lat = 32, zoom = 4) 
  })
}


shinyApp(ui, server)