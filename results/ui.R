library(shiny)
library(DT)
library(data.table)

navbarPage("Query DriverPower",
           tabPanel("Feature Selection",
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
                    ),
                    radioButtons('response', "Choose a response:",
                                choices = c("gmean", "count"), selected = 'gmean'),
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
                        plotOutput('qqplot'),
                        verbatimTextOutput("info")
                      )
                    )
           )
)
