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
             tabPanel("Map", fluidPage(
               fluidRow(
                 column(2, 
                        absolutePanel(top = 20, left = 30,checkboxInput("phy", "Load Phy", TRUE)),
                        absolutePanel(top = 50, left = 30,checkboxInput("heat", "Load Temp Anomoly", FALSE))
                        ),
                 column(10,
                        leafletOutput("mymap", height=700),
                        absolutePanel(top = 20, left = 70,  actionButton("resetMap", "Reset")),
                        )
               )
             )),
             tabPanel("Overview",includeMarkdown("readme.md")))
)


server <- function(input, output, session) {
  #import data
  data <- read.csv("datasets/curated/NOAAGlobalTemp/testdat.csv")
  phydata <- read.csv("outputs/quick_alr/latlon_north_america_phy.csv")
  #define the color pallate for the magnitidue of the earthquake
  pal <- colorNumeric(
    palette = c('blue', 'deep skyblue', 'cyan', 'orange red', 'red', 'dark red'),
    domain = data$mag)
  
  
  set_phy_color <- function(phyStr) {
    def_color <- "blue"
    if(phyStr == NULL) {
      def_color <- "red"
    }
    return(def_color)
  }
  
  output$mymap <- renderLeaflet({
    leaflet(data) %>% 
      setView(lng = -99, lat = 32, zoom = 4)  %>%
      addProviderTiles("Esri.WorldImagery", group="Satellite Map") %>%
      addProviderTiles("CartoDB.DarkMatter", group="Dark Map") %>%
      addTiles(options = providerTileOptions(noWrap = FALSE), group="Street Map") %>%
      #addCircles(data = data, lat = ~ LAT, lng = ~ LON, weight = 1, radius = 3, popup = ~as.character(mag), label = ~as.character(paste0("Temp Anomoly: ", sep = " ", mag)), color = ~pal(mag), fillOpacity = 0.5) %>%
      #addCircles(data = phydata, lat = ~ lat, lng = ~ lon, weight = 1, radius = 3, popup = "test", label = "Test", color = "blue", fillOpacity = 0.5) %>%
      addLayersControl(baseGroups = c("Dark Map","Satellite Map","Street Map"), options = layersControlOptions(collapsed = TRUE))
      #fitBounds(260, -23, 90, 80) 
  })
  
  observe({
    proxy <- leafletProxy("mymap", data = phydata)
    proxy %>% clearMarkers()
    if (input$phy) {
      proxy %>% addCircles(data = phydata, lat = ~ lat, lng = ~lon, weight = 1, radius = 3, popup = ~as.character(loc), label = ~as.character(paste0("Phy: ", sep = " ", loc)), color = ~set_phy_color(loc), fillOpacity = 0.5)
      }
    else {
      proxy %>% clearShapes()
      }
  })
  
  observe({
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% clearMarkers() 
    if (input$heat) {
      proxy %>%  addHeatmap(lng=~LON-270, lat=~LAT, intensity = ~mag, blur =  70, max = 10, radius = 40)
    }
    else{
      proxy %>% clearHeatmap()
    }
  })
  
  observeEvent(input$resetMap, {
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% setView(lng = -99, lat = 32, zoom = 4)
  })
}


shinyApp(ui, server)
