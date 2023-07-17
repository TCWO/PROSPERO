# Install new packages if missing --------------------------------------------------------

list.of.packages <- c("shiny", 'DT', 'tidyverse')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(tidyverse)
library(shiny)
library(ggplot2)

df_data <- read.csv2(file.path(wd, 'Data', 'df_data.csv'), stringsAsFactors = FALSE)
df_data <- df_data %>% gather(marker,value, BAX:Vimentin)
tmp <- df_data %>% mutate(ID = paste(sample.id, region, treatment, label, sep='_'))
tempdir(check=TRUE)



dirs <- list.dirs(path = "some/path")

library("DT")
##################################################################################
ui <- fluidPage(
  headerPanel('plotting 2D histogram of a selected sample using Sox2 and CD45'),
  sidebarPanel(
    selectInput('sample', 'Select a sample to proceed', unique(tmp$ID),
                selected = unique(tmp$ID)[1]),
    # Button
    downloadButton("downloadData", "Download")
    
  ),
  mainPanel(
    plotOutput("plot1", brush = "plot_brush"),
    verbatimTextOutput("info")
  )
)

server <- function(input, output) {
  selectedData <- reactive({
    tmp %>% filter(ID==input$sample) %>% spread(marker,value) %>% mutate(Sox2h = asinh(Sox2),CD45h = asinh(CD45))
  })
  
  options(width = 100) # Increase text width for printing table
  output$plot1 <- renderPlot({
    ggplot(selectedData(), aes(Sox2h,CD45h), fill = log10(..count..)) + geom_bin2d(bins=c(100,100)) + viridis::scale_fill_viridis() + theme_bw()
  })
  

  result <- reactive({
    brushedPoints(selectedData(), input$plot_brush, allRows = TRUE)
  })
  
  
  output$info <- renderPrint({
    result()
  })
  
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$sample, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(result(), file, row.names = FALSE)
    }
  )
  
}


shinyApp(ui, server)

