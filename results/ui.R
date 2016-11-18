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
                        "Panc-AdenoCA"
                      )
                    ),
                    radioButtons('response', "Choose a response:",
                                choices = c("gmean", "count"), selected = 'gmean'),
                    DT::dataTableOutput("table")
           ),
           tabPanel("Potential Drivers",
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
                        "Panc-AdenoCA"
                      ), selected = "ColoRect-AdenoCA"
                    ),
                    radioButtons('response2', "Choose a response:",
                                 choices = c("gmean", "count"), selected = 'gmean'),
                    radioButtons('func', "Choose a functional score:",
                                 choices = c("cadd", "eigen"), selected = 'eigen'),
                    tabsetPanel(
                      tabPanel(
                        "Table",
                        DT::dataTableOutput('res')
                      ),
                      tabPanel(
                        "Plot",
                        plotOutput('qqplot')
                      )
                    )
           )
)
