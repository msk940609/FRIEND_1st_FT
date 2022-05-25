library(shiny)
library(ggplot2)
library(RCurl)
library(data.table)
library(plotly)


ui <- pageWithSidebar(
  headerPanel('Distribution of assigned mass spectrum'),
  sidebarPanel(
    selectInput('Sample1', 'Select sample', as.character(unique(data$Sample))),
    selectInput('Sample2', 'Select sample', as.character(unique(data$Sample)))
  ),
  mainPanel(
    #Setting the multiplot size as 800x600 for better visibility
    plotlyOutput('plot1',width = "800px", height = "600px"),
    plotlyOutput('plot2',width = "800px", height = "600px")
  )
)


server <- shinyServer(
  function(input, output) {
    
    urlfile<-'https://raw.githubusercontent.com/msk940609/MSneutralloss4W/main/Datafile/FT_merge_4th.csv'
    
    data<-fread(urlfile)
    
    dat1 <- reactive({data[data$Sample == input$Sample1,]})
    dat2 <- reactive({data[data$Sample == input$Sample2,]})
    
   
    output$plot1 <- renderPlotly({
      plotheight <- 600
      fig <- plot_ly(dat1(), color = I("gray20"),
                     width = 780) %>% add_segments(x = ~`Calc m/z`, xend = ~`Calc m/z`, 
                                                   y = 0,yend = ~`Bromo Inty`, showlegend = T,
                                                   hoverinfo = "text", text = ~paste('<br>',  'Formula: ', `Formula`,
                                                                                     '<br>', 'Mass: ', `Calc m/z`,
                                                                                     sep = '')) %>% 
        layout(
          title=input$Sample1,
          margin = list(l = 50, r= 20, b = 70, t = 50),
          xaxis =  list(title = 'mz',
                        showgrid = F,
                        showline = T,
                        zeroline = F,
                        nticks = 5,
                        font = list(size = 8),
                        ticks = "outside",
                        ticklen = 5,
                        tickwidth = 2,
                        tickcolor = toRGB("black")
          )
        )
      
      
      fig <- fig %>% config(
        displaylogo = FALSE
      )
      
      fig
      
      
    })
     
    output$plot2 <- renderPlotly({
      plotheight <- 600
      fig <- plot_ly(dat2(), color = I("gray20"),
                     width = 780) %>% add_segments(x = ~`Calc m/z`, xend = ~`Calc m/z`, 
                                                   y = 0,yend = ~`Bromo Inty`, showlegend = T,
                                                   hoverinfo = "text", text = ~paste('<br>',  'Formula: ', `Formula`,
                                                                                     '<br>', 'Mass: ', `Calc m/z`,
                                                                                     sep = '')) %>% 
        layout(
          title=input$Sample2,
          margin = list(l = 50, r= 20, b = 70, t = 50),
          xaxis =  list(title = 'mz',
                        showgrid = F,
                        showline = T,
                        zeroline = F,
                        nticks = 5,
                        font = list(size = 8),
                        ticks = "outside",
                        ticklen = 5,
                        tickwidth = 2,
                        tickcolor = toRGB("black")
          )
        )
      
      
      fig <- fig %>% config(
        displaylogo = FALSE
      )
      
      fig
      
      
    })
  }
)

shinyApp(ui, server)
