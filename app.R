library(shiny)
library(leaflet)
library(leaflet.extras)
library(dplyr)
library(DT)

# the ui object has all the information for the user-interface
ui <-(
  navbarPage("Viral SpaceTime", id="main",
             tabPanel("Interactive Map", fluidPage(
               fluidRow(
                 column(2, 
                        absolutePanel(top = 20, left = 30,checkboxInput("phy", "Phylogenies", TRUE)),
                        absolutePanel(top = 50, left = 30,checkboxInput("heat", "Temperature Anomoly", FALSE))
                        ),
                 column(10,
                        leafletOutput("mymap", height=700),
                        absolutePanel(top = 20, left = 70,  actionButton("resetMap", "Reset")),
                        )
               )
             )),
             tabPanel("Phylogenies Data", DT::dataTableOutput("phytable")),
             tabPanel("About",includeMarkdown("app_readme.md"))
             )
)


server <- function(input, output, session) {
  #import data
  data <- read.csv("datasets/curated/NOAAGlobalTemp/testdat.csv")
  phydata <- read.csv("outputs/quick_alr/latlon_north_america_phy.csv", header=TRUE, stringsAsFactors=FALSE)
  str(phydata)
  
  set_phy_color <- function(phyloc) {
    phycolor <- "green"
    # add logic
    return(phycolor)
  }
  
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
      proxy %>% addCircles(data = phydata, lat = ~lat, lng = ~lon, weight = 1, radius = 3, popup = ~as.character(loc), label = ~as.character(paste0("Phy: ", sep = " ", loc)), color = ~set_phy_color(loc), fillOpacity = 0.5)
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
