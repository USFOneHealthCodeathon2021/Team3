setwd("C:/Users/myPC/Documents/usf_codeathon")

library(shiny)
library(plotly)
library(tidyverse)
library(shinythemes)
library(leaflet)
library(leaflet.extras)
library(dplyr)

#the ui object has all the information for the user-interface
ui <- ( navbarPage("Viral Space Time", id="main", tabPanel("mymap",
            fluidPage(fluidRow(column(2, absolutePanel(top = 20, left = 30, checkboxInput("phy", "Load Phy", TRUE)),
                        absolutePanel(top = 50, left = 30,  checkboxInput("heat", "Load Temp Anomoly", FALSE))),
                 column(10,leafletOutput("mymap", height=700),absolutePanel(top = 20, left = 70, actionButton("resetMap", "Reset"))))))))

            
# tabPanel("Overview",includeMarkdown("readme.md")))


# r_colors <- rgb(t(col2rgb(colors()) / 255))
# names(r_colors) <- colors()
# 
# ui <- fluidPage(
#   leafletOutput("mymap"),
#   p(),
#   actionButton("recalc", "New points")
# )

server <- function(input, output, session) {
  #import data
  data <- read.csv ("testdata.csv")
  phydata <- read.csv("latlon_north_america_phy.csv")
  
  #define the color pallate for the magnitidue of the earthquake
  pal <- colorNumeric( palette = c('blue', 'deep skyblue', 'cyan', 'orange red', 'red', 'dark red'), domain = data$mag)
  
  
  set_phy_color <- function(phyStr) {
    def_color <- "red"
    if(is.null(phyStr)) {def_color <- "red" } else {'green'}
       return(def_color)}
  
  output$mymap <- renderLeaflet({
    leaflet(data) %>% 
      setView(lng = -99, lat = 32, zoom = 5)  %>%
      addProviderTiles("Esri.WorldImagery", group="Satellite Map") %>%
      addProviderTiles("CartoDB.DarkMatter", group="Dark Map") %>%
      addTiles(options = providerTileOptions(noWrap = TRUE), group="Street Map") %>%
      #addCircles(data = data, lat = data$lat, lng = data$lon, weight = 1, radius = 3, popup = ~as.character(mag), label = ~as.character(paste0("Temp Anomoly: ", sep = " ", mag)), color = ~pal(mag), fillOpacity = 0.5) %>%
      #addCircles(data = phydata, lat = ~ lat, lng = ~ lon, weight = 1, radius = 3, popup = "test", label = "Test", color = "blue", fillOpacity = 0.5) %>%
      addLayersControl(baseGroups = c("Street Map","Dark Map","Satellite Map"), options = layersControlOptions(collapsed = TRUE))
    #fitBounds(260, -23, 90, 80) 
  })
  
  observe({
    proxy <- leafletProxy("mymap", data = phydata) %>%
      clearMarkers()  %>%
  addCircles(data = phydata, lat = phydata$lat, lng = phydata$lon, weight = 9, radius = 30,
                 popup = ~as.character(query), label = ~as.character(paste0("Phy: ", sep = " ", query)),
                 color = ~set_phy_color(query), fillOpacity = 0.5)})
   else {addHeatmap(lng=~lon, lat=~lat, intensity = ~mag, blur =  10, max = 0.05, radius = 15)}

 # observe({
 #  proxy <- leafletProxy("mymap", data = data) %>% proxy %>% clearMarkers()
 #     if (input$heat) { proxy %>%  addHeatmap(lng=~lon, lat=~lat, intensity = ~mag, blur =  10, max = 0.05, radius = 15)} 
 #     else{ proxy %>% clearHeatmap()}
 #  })
 # 
}



 # this is working
#   observe({
#         proxy <- leafletProxy("mymap", data = phydata) %>%
#          clearMarkers() %>% 
#            addCircles(data = phydata, lat = phydata$lat, lng = phydata$lon, weight = 9, radius = 30, 
#                      popup = ~as.character(query), label = ~as.character(paste0("Phy: ", sep = " ", query)), 
#                      color = ~set_phy_color(query), fillOpacity = 0.5)})
#       
#   
# } 
  
   
#   observe({
#     proxy <- leafletProxy("mymap", data = phydata)
#     proxy %>% clearMarkers() %>% 
#     if (input$phy) {
#       proxy %>% addCircles(data = phydata, lat = ~ lat, lng = ~lon, weight = 1, radius = 3, popup = ~as.character(loc), label = ~as.character(paste0("Phy: ", sep = " ", loc)), color = ~set_phy_color(loc), fillOpacity = 0.5)
#     }
#     else {
#       proxy %>% clearShapes()
#     }
#   })
#   
#   observe({
#     proxy <- leafletProxy("mymap", data = data)
#     proxy %>% clearMarkers() 
#     if (input$heat) {
#       proxy %>%  addHeatmap(lng=~LON-270, lat=~LAT, intensity = ~mag, blur =  70, max = 10, radius = 40)
#     }
#     else{
#       proxy %>% clearHeatmap()
#     }
#   })
#   
#   observeEvent(input$resetMap, {
#     proxy <- leafletProxy("mymap", data = data)
#     proxy %>% setView(lng = -99, lat = 32, zoom = 4)
#   })

shinyApp(ui = ui, server = server)

