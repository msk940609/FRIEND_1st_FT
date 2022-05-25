library(shiny)
library(plotly)
library(dplyr)
library(data.table)
options(shiny.maxRequestSize=300*1024^2)
getOption("digits")
options("digits" = 15)
library(data.table)
library(RColorBrewer)
options(scipen=1000)
library(rsconnect)
library(tidyverse)
library(bit64)
library(shinyFiles)
library(xlsx)
library(writexl)
library(ggplot2)
library(writexl)
library(xlsx)
library(vegan)
options(java.parameters = "- Xmx8192m")

Sys.setlocale("LC_TIME","english")

ui <- fluidPage(
  fluidPage(
    titlePanel("MY CSV FILES MERGER"), ###title of tab panel
    sidebarLayout(
      sidebarPanel(  ###start design side bar panel
        fileInput("file1",
                  "Choose CSV files from directory",
                  multiple = TRUE,
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')), ###select csv file (file format was converted as `xlsx` to `csv`) 
        selectInput(inputId ="k", h4("dimension"),   
                    choices = c(1,2,3,4,5), selected = 2), ####select dimension (now indevelopment)
        selectInput(inputId ="dist", h4("dist method"), 
                    choices = c("euclidean","bray","jaccard"), selected = "bray"), ###select distance method
        downloadButton('downloadData', 'Download') ##download merge file
      ),
      mainPanel(
        plotOutput('plotMDS',width = "600px", height = "400px"), ##display nmds plot
        tableOutput('contents') ####display ordination of nmds plot point
      )
    )
  )
)

server <-  function(input, output) {
  
  getData <- reactive({  ###select file and merge
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }else {
      # browser()
      numfiles = nrow(inFile) 
      dt = list()
      for (i in 1:numfiles)
      {
        
        temp = fread(input$file1[[i, 'datapath']]) ##Read selected file
        
        dt[[i]] = temp
        }
       #browser()
      do.call(rbind, dt) ###merge selected file
      
    }
  })
  
  
  nmds <- reactive({
    
    ft_data=getData() ###load merge data
    
    MDS=melt(ft_data[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
      dcast(Sample~Formula, sum) ###build nmds data

    NMDS=metaMDS(MDS[,-c(1)], k=2, distance = input$dist, trymax = 20) ####Calculate distance and ordination
    
    gnmds=as.data.table(NMDS$points) ###excert nmds point
    gnmds$Sample=MDS$Sample ###add sample information
    
    gnmds=gnmds %>% tidyr::separate("Sample",c("Group","no"),sep="_") ###separate sample to group and no
  })
  
  output$contents <- renderTable( 
    nmds() ###show nmds point table
  )
  
  output$plotMDS <- renderPlot({
    
    p1=ggplot()+
      geom_point(data=nmds(), aes(x=MDS1, y=MDS2, col=Group)) ###display nmds point
    
    p1
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      paste("data-", Sys.Date(), ".csv", sep="") 
    },
    content = function(file) { 
      write.csv(getData(), file, row.names=FALSE)   
    })
}

shinyApp(ui = ui, server = server)
