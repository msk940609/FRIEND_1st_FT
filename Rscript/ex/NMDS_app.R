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
    titlePanel("MY CSV FILES MERGER"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file1",
                  "Choose CSV files from directory",
                  multiple = TRUE,
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        selectInput(inputId ="k", h4("dimension"), 
                    choices = c(1,2,3,4,5), selected = 2),
        selectInput(inputId ="dist", h4("dist method"), 
                    choices = c("euclidean","bray","jaccard"), selected = "bray"),
        downloadButton('downloadData', 'Download')
      ),
      mainPanel(
        plotOutput('plotMDS',width = "600px", height = "400px"),
        tableOutput('contents')
      )
    )
  )
)

server <-  function(input, output) {
  
  getData <- reactive({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }else {
      # browser()
      numfiles = nrow(inFile) 
      dt = list()
      for (i in 1:numfiles)
      {
        
        temp = fread(input$file1[[i, 'datapath']])
        
        dt[[i]] = temp
        }
       #browser()
      do.call(rbind, dt)
      
    }
  })
  
  
  nmds <- reactive({
    
    ft_data=getData()
    
    MDS=melt(ft_data[,c("Sample","Formula","Bromo.Inty")], id.vars = c("Sample","Formula")) %>% 
      dcast(Sample~Formula, sum)

    NMDS=metaMDS(MDS[,-c(1)], k=2, distance = input$dist, trymax = 20)
    
    gnmds=as.data.table(NMDS$points)
    gnmds$Sample=MDS$Sample
    
    gnmds=gnmds %>% tidyr::separate("Sample",c("Group","no"),sep="_")
  })
  
  output$contents <- renderTable( 
    nmds()
  )
  
  output$plotMDS <- renderPlot({
    
    #new_nmds=nmds
    #NMDS=vegan::metaMDS(new_nmds[,-c(1)], k=2, distance = "bray", trymax = 100)
    #gnmds=as.data.table(NMDS$points)
    #gnmds$Sample=new_nmds$Sample
    #gnmds=gnmds %>% tidyr::separate("Sample",c("Group","no"),sep=" ")
    
    p1=ggplot()+
      geom_point(data=nmds(), aes(x=MDS1, y=MDS2, col=Group))
    
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
