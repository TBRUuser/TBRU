###########################
##### Install Package #####
###########################
if (!require('shiny')) install.packages('shiny')
if (!require('ggvis')) install.packages('ggvis')
if (!require('dplyr')) install.packages('dplyr')
if (!require('ggplot2')) install.packages('ggplot2')
if (!require('cowplot')) install.packages('cowplot')
if (!require('plotly')) install.packages('plotly')
if (!require('DT')) install.packages('DT')
if (!require('reshape2')) install.packages('reshape2')

library(shiny)
library(ggvis)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(DT)
library(reshape2)


######################
##### Data Entry #####
######################
## Dataset2 - Original Function##
path="./data/dataset0.csv"
data <- read.csv(path)
nr <- nrow(data)


####################
##### Shiny UI #####
####################
shinyUI(
  navbarPage("TBRU",
             ## 1.1 --- introduction ## ------------------------------------------------------------------------------------------------
             tabPanel("Introduction",
                      br(),
                      #               HTML("<h1><b><center> T B R U </center></b></h1>"),
                      includeHTML("intro_tx.html")
                      #               img(src='introduction.png')
                      #               imageOutput("Intro1")
             ),
             ## 1.2 --- Data Exploration ## --------------------------------------------------------------------------------------------
             navbarMenu("Data Exploration", 
                        ## 1.2.1 --- Summary ## -------------------------------------------------------------------------
                        tabPanel("Summary",
                                 sidebarPanel(
                                   h2("Research Backgroud"),
                                   br(),
                                   h3("This research ..."),br(),
                                   h3("1. Purpose"),p(),
                                   h3("2. Method"),p(),
                                   h3("3. Result"),p(),
                                   h3("4. Improvement"),br(),
                                   h2("Dataset"),
                                   h4("dataset0 : Test dataset, also use for testing TBRU app"),
                                   h4("dataset2 : Real dataset, the first version"),
                                   h4("dataset3 : Real dataset, the second version (latest)"),
                                   selectInput("dataset","Dataset:",dir("./data")),
                                   downloadButton('download_data','Download')
                                 ),
                                 mainPanel(
                                   tabsetPanel(
                                     ## 1.2.1.1 --- Dataset ## ---------------------------------------------
                                     tabPanel("Dataset", 
                                              br(),       
                                              DT::dataTableOutput("rawdata")
                                     ),
                                     ## 1.2.1.2 --- Level ## -----------------------------------------------
                                     tabPanel("Table",
                                              br(),       
                                              h4("Dilution vs. Visit"),
                                              tableOutput("tb_vis_dilu"),     
                                              br(),
                                              h4("Peptides"),
                                              tableOutput("tb_sample")
                                     ),
                                     ## 1.2.1.3 --- Distribution ## ----------------------------------------
                                     tabPanel("Distribution",
                                              br(),
                                              plotlyOutput("histo2"),
                                              br(),br(),
                                              plotlyOutput("histo1")
                                     ),
                                     ## 1.2.1.4 --- CV Analysis ## -----------------------------------------
                                     tabPanel("Valiability",
                                              h3("Box Plot"),
                                              sidebarPanel(
                                                radioButtons("box","Variability",c("inter-1","inter-2","inter-3"))
                                              ),
                                              mainPanel(
                                                plotlyOutput("box1")
                                              )
                                     )
                                   )
                                 )
                        ),
                        ## 1.2.2 --- Plot ## ----------------------------------------------------------------------------
                        tabPanel("Plot",
                                 tabsetPanel(
                                   ## 1.2.2.1 --- Boxplot ## ---------------------------------------------
                                   tabPanel("Boxplot",
                                            br(),
                                            fluidRow(
                                              column(3,
                                                wellPanel(
                                                  uiOutput("ui_dilution_box2"),
                                                  uiOutput("ui_visit_box2")
#                                                  selectInput("dilution_b2","Dilution value:",levels(factor(data$Dilution))),
#                                                  selectInput("visit_b2","Visit:",levels(factor(data$Visit)))
                                                )
                                              ),
                                              column(9,
                                                plotlyOutput("box2",height = "800px")
                                              )
                                            )
                                   ),
                                   ## 1.2.2.2 --- Variety across Visit ## --------------------------------
                                   tabPanel("Variety across Visit",
                                            br(),
#                                            helpText("This plot can help to see the difference among ID. Some people will 
#                                                     have a huge difference in IFNg when time past, but others may not."),
                                            fluidRow(
                                              column(3,
                                                wellPanel(
                                                  uiOutput("ui_dilution_vav"),
                                                  uiOutput("ui_ID_vav")
#                                                  selectInput("dilution_vav","Dilution value:",levels(factor(data$Dilution))),
#                                                  selectInput("ID_vav","Present ID  (across PEPTIDE):",levels(factor(data$ID)))
                                                  
                                                )
                                              ),
                                              column(9,
                                                plotlyOutput("plot3",height = "800px")
                                              )
                                            )
                                   ),
                                   ## 1.2.2.3 --- Heatmap ## ---------------------------------------------
                                   tabPanel("Heatmap",
                                           br(),
                                           fluidRow(
                                             column(1),
                                             column(3,
                                                    uiOutput("ui_dilution_hmp1")
#                                                    selectInput("dilution_hmp1","Dilution value:",levels(factor(data$Dilution)))
                                             ),
                                             column(3,
                                                    uiOutput("ui_visit_hmp1")
#                                                    selectInput("visit_hmp1","Visit:",levels(factor(data$Visit)))   
                                             )
                                           ),
                                           fluidRow(
                                             plotlyOutput("heatmp1",height = "900px",width = "1800px")
                                           ), 
                                           br(),
                                           fluidRow(
                                             column(1),
                                             column(3,
                                                    uiOutput("ui_dilution_hmp2")
#                                                    selectInput("dilution_hmp2","Dilution value:",levels(factor(data$Dilution)))
                                             ),
                                             column(3,
                                                    uiOutput("ui_visit_hmp2")
#                                                    selectInput("visit_hmp2","Visit:",levels(factor(data$Visit)))   
                                             )
                                           ),
                                           fluidRow(
                                             plotOutput("heatmp2",height = "800px",width = "1500px")
                                           )
                                   ),
                                   ## 1.2.2.4 --- Line Plot ## ---------------------------------------------
                                   tabPanel("Line Plot",
                                            fluidRow(
                                              wellPanel(
                                                radioButtons("across_type","Across:",c("Peptides","Individuals"),inline = T)
                                              )
                                            ),
                                            uiOutput("ui_title_lpall"),
                                            br(),
                                            fluidRow(
                                              column(3,
                                                wellPanel(
                                                  uiOutput("ui_dilution_lpall"),
                                                  uiOutput("ui_others_lpall")
                                                )
                                              ),
                                              column(1),
                                              column(8,
                                                uiOutput("ui_plot_lpall")
                                              )
                                            ),
                                            br(),
                                            uiOutput("ui_title_lpsingle"),
                                            fluidRow(
                                              column(3,
                                                wellPanel(
                                                  uiOutput("ui_dilution_lpsingle"),
                                                  uiOutput("ui_others_lpsingle"),
                                                  br(),
                                                  checkboxInput("select_all","Select All"),
                                                  uiOutput("ui_checkbox_lpsingle")
                                                )
                                              ),
                                              column(1),
                                              column(7,
                                                uiOutput("ui_plot_lpsingle")
                                              )
                                            )
                                            
                                   )  # tabPanel - Line Plot - 1.2.2.4
                                 )  # tabsetPanel 
                        )  # tabPanel - Plot - 1.2.2
             ),  # narbarMenu - DE
             navbarMenu("Data Analysis",
               tabPanel("test",
                 h3("Will update data analysis result here, once the research have progress")
               )
             )
  )  # narbarPanel - TBRU
)  # shinyUI
