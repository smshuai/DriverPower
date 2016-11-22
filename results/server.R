library(shiny)
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
    sep = '\t'
  )
  meta = meta[, c('feature_name' , 'description')]
  # read project overview
  project = read.table('./project_overview.tsv', header=TRUE, row.names = 'Project_Code')
  
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
    dat = dat[, c('binID', 'Length', 'nMut', 'nSample', 'MuAdj', 'fscore', 'Pval', 'Qval')]
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10( 1:length(dat$o)/length(dat$o))
    dat
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
    genes = dat[dat$Qval<=0.1, ]
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  
  output$summary <- renderText({
    dat = resDat()
    genes = dat[dat$Qval<=0.1, ]
    paste0(
      "Number of samples: ",
      as.character(project[input$tumor2, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$tumor2, 'WG_Mu'])
      )
    
  })
  
  output$pt_dp <- renderPrint({
    dat = resDat()
    
    nearPoints(dat, input$plot_click, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$br_dp <- renderPrint({
    dat = resDat()
    brushedPoints(dat, input$plot_brush)
  })
  
  # oncodriveFML (odf)
  odfDat <- reactive({
    dat = read.table(paste0(input$odfType,'/oncodriveFML/',input$oncodrive,'/',input$odfTumor,'.CDS.oncodriveFML.observed.txt'),
                     header=TRUE, sep='\t')
    dat[, 'element_ID'] = tstrsplit(dat$element_ID, '::', fixed = T)[[3]]
    dat = dat[, c('element_ID', 'p.value', 'q.value')]
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$odfTab <- DT::renderDataTable({
    odfDat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))
  output$odfPlot <- renderPlot({
    dat = odfDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$odfSummary <- renderText({
    dat = odfDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$odfTumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$odfTumor, 'WG_Mu'])
    )
  })
  
  output$odfPt <- renderPrint({
    dat = odfDat()
    
    nearPoints(dat, input$odfClick, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$odfBr <- renderPrint({
    dat = odfDat()
    brushedPoints(dat, input$odfBrush)
  })
  
  # ActiveDriver2 (ad2)
  ad2Dat <- reactive({
    dat = read.table(paste0(input$ad2Type,'/ActiveDriver2/',input$ad2Tumor,'.gc19_pc.cds.ActiveDriver2-burden.observed.txt'),
                     header=TRUE, sep='\t')
    dat[, 'id'] = tstrsplit(dat$id, '::', fixed = T)[[3]]
    dat[,'Qval'] = p.adjust(dat$raw_p, method = 'BH')
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    rownames(dat) = 1:nrow(dat)
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$ad2Tab <- DT::renderDataTable({
    ad2Dat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))  
  output$ad2Plot <- renderPlot({
    dat = ad2Dat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$ad2Summary <- renderText({
    dat = ad2Dat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$ad2Tumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$ad2Tumor, 'WG_Mu'])
    )
  })
  
  output$ad2Pt <- renderPrint({
    dat = ad2Dat()
    
    nearPoints(dat, input$ad2Click, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$ad2Br <- renderPrint({
    dat = ad2Dat()
    brushedPoints(dat, input$ad2Brush)
  })
  
  # compositeDriver (cd)
  cdDat <- reactive({
    dat = read.table(paste0(input$cdType,'/compositeDriver/',input$cdTumor,'.CDS.compositeDriver.observed.txt'),
                     header=TRUE, sep='\t')
    dat[, 'element_ID'] = tstrsplit(dat$element_ID, '::', fixed = T)[[3]]
    dat[,'Qval'] = p.adjust(dat$pvalue, method = 'BH')
    dat = dat[, c('element_ID', 'pvalue', 'Qval')]
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    rownames(dat) = 1:nrow(dat)
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$cdTab <- DT::renderDataTable({
    cdDat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))  
  output$cdPlot <- renderPlot({
    dat = cdDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$cdSummary <- renderText({
    dat = cdDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$cdTumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$cdTumor, 'WG_Mu'])
    )
  })
  
  output$cdPt <- renderPrint({
    dat = cdDat()
    
    nearPoints(dat, input$ad2Click, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$cdBr <- renderPrint({
    dat = cdDat()
    brushedPoints(dat, input$cdBrush)
  })
  
  # ExInAtor (eia)
  eiaDat <- reactive({
    dat = read.table(paste0(input$eiaType,'/ExInAtor/',input$eiaTumor,'.CDS.ExInAtor.observed.txt'),
                     header=TRUE, sep='\t', check.names = FALSE)
    dat[, 'element_ID'] = tstrsplit(dat$element_ID, '::', fixed = T)[[3]]
    dat[,'Qval'] = p.adjust(dat[,'p-value'], method = 'BH')
    dat = dat[, c('element_ID', 'p-value', 'Qval')]
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    rownames(dat) = 1:nrow(dat)
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$eiaTab <- DT::renderDataTable({
    eiaDat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))  
  output$eiaPlot <- renderPlot({
    dat = eiaDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$eiaSummary <- renderText({
    dat = eiaDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$eiaTumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$eiaTumor, 'WG_Mu'])
    )
  })
  
  output$eiaPt <- renderPrint({
    dat = eiaDat()
    
    nearPoints(dat, input$eiaClick, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$eiaBr <- renderPrint({
    dat = eiaDat()
    brushedPoints(dat, input$eiaBrush)
  })
  
  # ncdDetect (ncdd)
  ncddDat <- reactive({
    dat = read.table(paste0(input$ncddType,'/ncdDetect/',input$ncddTumor,'.CDS.ncdDetect.observed.txt'),
                     header=TRUE, sep='\t', check.names = FALSE)
    dat[, 'element_ID'] = tstrsplit(dat$element_ID, '::', fixed = T)[[3]]
    dat[,'Qval'] = p.adjust(dat[,'p_value'], method = 'BH')
    dat = dat[, c('element_ID', 'p_value', 'Qval')]
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    rownames(dat) = 1:nrow(dat)
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$ncddTab <- DT::renderDataTable({
    ncddDat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))  
  output$ncddPlot <- renderPlot({
    dat = ncddDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$ncddSummary <- renderText({
    dat = ncddDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$ncddTumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$ncddTumor, 'WG_Mu'])
    )
  })
  
  output$ncddPt <- renderPrint({
    dat = ncddDat()
    
    nearPoints(dat, input$ncddClick, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$ncddBr <- renderPrint({
    dat = ncddDat()
    brushedPoints(dat, input$ncddBrush)
  })
  
  # ncDriver (ncd)
  ncdDat <- reactive({
    dat = read.table(paste0(input$ncdType,'/ncDriver/',input$ncdTumor,'.gc19_pc.cds.ncDriverConservation.observed.snv.p.comb'),
                     header=TRUE, sep='\t', check.names = FALSE)
    dat[, 'id.element'] = tstrsplit(dat$id.element, '::', fixed = T)[[3]]
    dat[,'Qval'] = p.adjust(dat[,'conservation.p.comb'], method = 'BH')
    dat = dat[, c('id.element', 'conservation.p.comb', 'Qval')]
    colnames(dat) = c('binID', 'Pval', 'Qval')
    # order by pvalue
    dat = dat[order(dat$Pval, decreasing = F),]
    rownames(dat) = 1:nrow(dat)
    dat[, 'o'] = -log10(dat$Pval)
    dat[, 'e'] = -log10(1:length(dat$o)/length(dat$o))
    dat
  })
  output$ncdTab <- DT::renderDataTable({
    ncdDat()
  }, options = list(rownames = TRUE,
                    order = list(3, 'asc')))  
  output$ncdPlot <- renderPlot({
    dat = ncdDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    p = ggplot(data = dat, aes(x=e, y=o)) + geom_point() + geom_abline(intercept=0,slope=1, col="red")
    p = p + xlab(expression(Expected~~-log[10](italic(p)))) +
      ylab(expression(Observed~~-log[10](italic(p)))) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p = p + geom_text(data = genes, aes(x=e, y=o, label = binID), check_overlap = TRUE, hjust = 0, nudge_x = 0.1)
    p
  })
  output$ncdSummary <- renderText({
    dat = ncdDat()
    genes = dat[dat$Qval<=0.1, ]
    genes = na.omit(genes)
    paste0(
      "Number of samples: ",
      as.character(project[input$ncdTumor, 'Num_Donors']),
      "\nNumber of significant (q<0.1) bins: ",
      as.character(nrow(genes)),
      "\nOverall mutation rate (/Mbp): ",
      as.character(project[input$ncdTumor, 'WG_Mu'])
    )
  })
  
  output$ncdPt <- renderPrint({
    dat = ncdDat()
    
    nearPoints(dat, input$ncdClick, threshold = 10, maxpoints = 1,
               addDist = FALSE)
  })
  
  output$ncdBr <- renderPrint({
    dat = ncdDat()
    brushedPoints(dat, input$ncdBrush)
  })
  
}


