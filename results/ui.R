library(shiny)
library(DT)
library(data.table)

navbarPage("Query DriverPower",
           tabPanel("Feature Selection",
                    fluidRow(
                      column(
                        3,
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
                            "Panc-AdenoCA",
                            'Lung-AdenoCA'
                          )
                        )),
                      column(
                        3,
                        radioButtons('response', "Choose a response:",
                                     choices = c("gmean", "count"), selected = 'gmean')
                        )
                    ),
                    DT::dataTableOutput("table")
           ),
           tabPanel("DriverPower",
                    fluidRow(
                      column(3,
                        selectInput(
                          "type",
                          "Choose a category:",
                          choices = c("exon"),
                          selected = 'exon'
                        )),
                      column(3,
                        selectInput(
                          "tumor2",
                          "Choose a tumor:",
                          choices = c(
                            "ColoRect-AdenoCA",
                            "Skin-Melanoma",
                            "Liver-HCC",
                            "Eso-AdenoCa",
                            "Lung-SCC",
                            "Stomach-AdenoCA",
                            "Panc-AdenoCA",
                            'Lung-AdenoCA'
                          ), selected = "ColoRect-AdenoCA"
                        )),
                      column(
                        3,
                        radioButtons('response2', "Choose a response:",
                                     choices = c("gmean", "count"), selected = 'gmean')
                        ),
                      column(
                        3,
                        radioButtons('func', "Choose a functional score:",
                                     choices = c("cadd", "eigen"), selected = 'eigen')
                        )
                    ),
                    tabsetPanel(
                      tabPanel(
                        "Table",
                        DT::dataTableOutput('res')
                      ),
                      tabPanel(
                        "Plot",
                        fluidRow(
                          column(
                            6,
                            plotOutput('qqplot', width = '600px', height = '600px', click = "plot_click", brush = "plot_brush")
                            ),
                          column(
                            6,
                            verbatimTextOutput("summary"),
                            verbatimTextOutput("pt_dp"),
                            verbatimTextOutput("br_dp")
                          )
                        )
                      )
                    )
           ),
           tabPanel(
             'OncodriveFML',
             fluidRow(
               column(3,
                      selectInput(
                        "odfType",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "odfTumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      )),
               column(
                 3,
                 radioButtons('oncodrive', "Choose a method:",
                              choices = c("cadd", "vest3"), selected = 'cadd')
               )
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('odfTab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('odfPlot', width = '600px', height = '600px', click = "odfClick", brush = "odfBrush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("odfSummary"),
                     verbatimTextOutput("odfPt"),
                     verbatimTextOutput("odfBr")
                   )
                 )
               )
             )
           ),
           tabPanel(
             'ActiveDriver2',
             fluidRow(
               column(3,
                      selectInput(
                        "ad2Type",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "ad2Tumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      ))
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('ad2Tab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('ad2Plot', width = '600px', height = '600px', click = "ad2Click", brush = "ad2Brush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("ad2Summary"),
                     verbatimTextOutput("ad2Pt"),
                     verbatimTextOutput("ad2Br")
                   )
                 )
               )
             )
           ),
           tabPanel(
             'compositeDriver',
             fluidRow(
               column(3,
                      selectInput(
                        "cdType",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "cdTumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      ))
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('cdTab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('cdPlot', width = '600px', height = '600px', click = "cdClick", brush = "cdBrush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("cdSummary"),
                     verbatimTextOutput("cdPt"),
                     verbatimTextOutput("cdBr")
                   )
                 )
               )
             )
           ),
           tabPanel(
             'ExInAtor',
             fluidRow(
               column(3,
                      selectInput(
                        "eiaType",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "eiaTumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      ))
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('eiaTab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('eiaPlot', width = '600px', height = '600px', click = "eiaClick", brush = "eiaBrush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("eiaSummary"),
                     verbatimTextOutput("eiaPt"),
                     verbatimTextOutput("eiaBr")
                   )
                 )
               )
             )
           ),
           tabPanel(
             'ncdDetect',
             fluidRow(
               column(3,
                      selectInput(
                        "ncddType",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "ncddTumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      ))
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('ncddTab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('ncddPlot', width = '600px', height = '600px', click = "ncddClick", brush = "ncddBrush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("ncddSummary"),
                     verbatimTextOutput("ncddPt"),
                     verbatimTextOutput("ncddBr")
                   )
                 )
               )
             )
           ),
           tabPanel(
             'ncDriver',
             fluidRow(
               column(3,
                      selectInput(
                        "ncdType",
                        "Choose a category:",
                        choices = c("exon"),
                        selected = 'exon'
                      )),
               column(3,
                      selectInput(
                        "ncdTumor",
                        "Choose a tumor:",
                        choices = c(
                          "ColoRect-AdenoCA",
                          "Skin-Melanoma",
                          "Liver-HCC",
                          "Eso-AdenoCa",
                          "Lung-SCC",
                          "Stomach-AdenoCA",
                          "Panc-AdenoCA",
                          'Lung-AdenoCA'
                        ), selected = "ColoRect-AdenoCA"
                      ))
             ),
             tabsetPanel(
               tabPanel(
                 'Table',
                 DT::dataTableOutput('ncdTab')
               ),
               tabPanel(
                 'Plot',
                 fluidRow(
                   column(
                     6,
                     plotOutput('ncdPlot', width = '600px', height = '600px', click = "ncdClick", brush = "ncdBrush")
                   ),
                   column(
                     6,
                     verbatimTextOutput("ncdSummary"),
                     verbatimTextOutput("ncdPt"),
                     verbatimTextOutput("ncdBr")
                   )
                 )
               )
             )
           )
)
