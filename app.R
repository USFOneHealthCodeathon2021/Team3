library(shiny)
library(plotly)
library(tidyverse)
library(shinythemes)
library(leaflet)
library(leaflet.extras)
library(dplyr)

# increase max R-Shiny user-input file size from 5 to 30 MB
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

# define fonts for plot
f1 <- list(family = "Arial, sans-serif",
           size = 24,
           color = "black")

f2 <- list(family = "Arial, sans-serif",
           size = 20,
           color = "black")

f3 <- list(family = "Arial, sans-serif",
           size = 16,
           color = "black")

# set plot margins
m <- list(
  l = 20,
  r = 20,
  b = 10,
  t = 100,
  pad = 4
)

# the ui object has all the information for the user-interface
ui <-(
  navbarPage("Viral Space Time", id="main",
             tabPanel("Map", leafletOutput("mymap", height=600), 
                      absolutePanel(top = 160, left = 80,checkboxInput("heat", "Heatmap", FALSE)
                      )),
             tabPanel("Overview",includeMarkdown("readme.md")))
)


server <- function(input, output, session) {
  #import data
  data <- read.csv("datasets/curated/NOAAGlobalTemp/testdat.csv")
  
  #define the color pallate for the magnitidue of the earthquake
  pal <- colorNumeric(
    palette = c('blue', 'deep skyblue', 'cyan', 'orange red', 'red', 'dark red'),
    domain = data$mag)
  
  output$mymap <- renderLeaflet({
    leaflet(data) %>% 
      setView(lng = -98.557, lat = 37.26 , zoom = 3) %>%
      # Add two tiles
      addProviderTiles("Esri.WorldImagery", group="Satellite Map") %>%
      addProviderTiles("CartoDB.DarkMatter", group="Dark Map") %>%
      addTiles(options = providerTileOptions(noWrap = TRUE), group="Street Map") %>%
      addCircles(data = data, lat = ~ LAT, lng = ~ LON, weight = 1, radius = 25000, popup = ~as.character(mag), label = ~as.character(paste0("Temp Anomoly: ", sep = " ", mag)), color = ~pal(mag), fillOpacity = 0.5) %>%
      addLayersControl(baseGroups = c("Dark Map","Satellite Map","Street Map"), options = layersControlOptions(collapsed = TRUE))
  })
  
  observe({
    proxy <- leafletProxy("mymap", data = data)
    proxy %>% clearMarkers()
    if (input$heat) {
      proxy %>%  addHeatmap(lng=~LON, lat=~LAT, intensity = ~mag, blur =  3, max = 0.05, radius = 7) 
    }
    else{
      proxy %>% clearHeatmap()
    }
  })
}


shinyApp(ui, server)
