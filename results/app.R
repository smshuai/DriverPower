#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(ggplot2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("View feature selection results"),
  sidebarLayout(
    # Sidebar with a slider input for number of bins
    sidebarPanel(
      selectInput(
        "result",
        "Choose a result:",
        choices = c('Feature selection',
                    'Coding drivers')
      ),
      
      selectInput(
        "tumor",
        "Choose a tumor:",
        choices = c(
          "ColoRect-AdenoCA",
          "Skin-Melanoma",
          "Liver-HCC",
          "Eso-AdenoCa",
          "Lung-SCC",
          "Stomach-AdenoCA",
          "Panc-AdenoCA"
        )
      ),
      
      selectInput('response', "Choose a response:",
                  choices = c("gmean", "count")),
      
      submitButton("Update View")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Table", 
                 h4("Observations"),
                 DT::dataTableOutput('dat')
                 ),
        tabPanel("Plot",
                 selectInput('type', "Choose a type:",
                             choices = c("rndlasso", "lasso", 'rho', 'freg')),
                 submitButton("Update Plot"),
                 plotOutput('barPlot'))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # read meta
  meta = read.table(
    '../feature_meta.tsv',
    header = TRUE,
    row.names = ,
    sep = '\t'
  )
  meta = meta[, c('feature_name' , 'description')]
  # Return the requested dataset
  datasetInput <- reactive({
    if(input$result == 'Feature selection'){
      dat = read.table(paste0("../freeze/", input$tumor, ".fs.", input$response, ".tsv"),
                       header = TRUE)
      dat[, 2:ncol(dat)] = abs(dat[, 2:ncol(dat)])
      dat = merge(dat, meta, by.x = 'fname', by.y = 'feature_name')
      row.names(dat) = dat$fname
      dat = dat[, -1]
      dat 
    } else if (input$result == 'Coding drivers'){
      
    }
  })
  
  # Show the first "n" observations
  output$dat <- DT::renderDataTable({
    datasetInput()
  }, options = list(rownames = TRUE,
                    order = list(4, 'desc')))
  
  # Make a barplot for top features
  output$barPlot <- renderPlot({
    dat = datasetInput()
    N = 10
    dat = head(dat[order(dat[,input$type], decreasing = T), ], N)
    # dat['fname'] = as.character(rownames(dat))
    # ggplot(dat, aes_string(x='fname', y=input$type)) + geom_bar(stat='identity')
    op <- par(mar = c(15,4,4,2) + 0.1)
    barplot(height = dat[,input$type], names.arg = rownames(dat),
            ylab = input$type, las=2, cex.lab=1.5, cex.axis=1.3)
  })
  
  # output$distPlot <- renderPlot({
  #   # generate bins based on input$bins from ui.R
  #   x    <- faithful[, 2]
  #   bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #
  #   # draw the histogram with the specified number of bins
  #   hist(x, breaks = bins, col = 'darkgray', border = 'white')
  # })
}

# Run the application
shinyApp(ui = ui, server = server)
