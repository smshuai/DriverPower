library(shiny)
source('./qqunif.plot.R')
library(ggplot2)
# Define server logic required to draw a histogram
function(input, output, session) {
  output$plot <- renderPlot({
    plot(cars)
  })
  

  # read meta
  meta = read.table(
    './feature_meta.v2.tsv',
    header = TRUE,
    row.names = ,
    sep = '\t'
  )
  meta = meta[, c('feature_name' , 'description')]
  
  # show data table
  ## feature selection
  fsDat <- reactive({
    dat = read.table(
      paste0(
        "./feature_selection/",
        input$tumor,
        ".",
        input$response,
        ".fs.tsv"
      ),
      header = TRUE
    )
    dat[, 2:ncol(dat)] = signif(abs(dat[, 2:ncol(dat)]))
    dat = merge(dat, meta, by.x = 'fname', by.y = 'feature_name')
    row.names(dat) = dat$fname
    dat = dat[,-1]
    dat
  })
  ## results table
  resDat <- reactive({
    dat = read.table(
      paste0(
        './',
        input$type,
        '/driverpower/',
        input$tumor2,
        '.',
        input$func,
        '.',
        input$response2,
        '.res.tsv'
      ),
      sep = '\t',
      header = TRUE
    )
    dat[, 'binID'] = tstrsplit(dat$binID, '::', fixed = T)[[3]]
    dat$MuAdj = signif(dat$MuAdj, 4)
    dat$fscore = signif(dat$fscore, 4)
    dat$Pval = signif(dat$Pval, 4)
    dat$Qval = signif(dat$Qval, 4)
    dat[, c('binID', 'Length', 'nMut', 'nSample', 'MuAdj', 'fscore', 'Pval', 'Qval')]
  })
  # feature selection table
  output$table <- DT::renderDataTable({
    fsDat()
  }, options = list(rownames = TRUE,
                      order = list(4, 'desc')))
  output$res <- DT::renderDataTable({
    resDat()
  }, options = list(rownames = FALSE))

  output$qqplot <- renderPlot({
    dat = resDat()
    qqunif.plot(dat$Pval)
  })
  
  output$info <- renderPrint({
    dat = resDat()
    nearPoints(dat, input$plot_click, threshold = 10, maxpoints = 1,
               addDist = TRUE)
  })
}


