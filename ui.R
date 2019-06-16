##---- A shiny webapp to input qPCR raw data and get Transcript Abundance and/or Fold Change

##---- load libraries ----
library(shiny)
library(DT)
library(data.table)
library(tidyverse)
library(magrittr)
library(stringr)
##---- load data ----


##---- UI ----
ui <- navbarPage('qPCR Analysis Tool',
                 ## Input Tab ----                 
                 tabPanel('Input',
                          titlePanel("qPCR Calculator"),
                          sidebarLayout(
                              sidebarPanel(
                                  helpText('Upload a .csv or .txt and select your reference gene to calculate Transcript Abundance and Fold Change',
                                           br(),
                                           'Input file must contain a gene, treatment, replicate and Ct column.'),
                                  fileInput("file",
                                            label = h4("Select a .csv or .txt file with your raw data"),
                                            multiple=FALSE,
                                            accept = c("text/csv",
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  selectInput('genecol',
                                              'Select Gene Column',
                                              choices='',
                                              selected='',
                                              multiple=FALSE),
                                  helpText('Reference gene may not work if names contain spaces'),
                                  selectInput("refgene",
                                              "Select Reference Gene",
                                              choices="",
                                              selected="",
                                              multiple=FALSE),
                                  selectInput('treatcol',
                                              'Select Treatment Column',
                                              choices='',
                                              selected='',
                                              multiple=FALSE),
                                  selectInput("controlsample",
                                              "Select Control Sample",
                                              choices="",
                                              selected="",
                                              multiple=FALSE),
                                  selectInput('repcol',
                                              'Select Replicate Column',
                                              choices='',
                                              selected='',
                                              multiple=FALSE),
                                  selectInput('ctcol',
                                              'Select Ct Column',
                                              choices='',
                                              selected='',
                                              multiple=FALSE),
                                  radioButtons("disp", "Display",
                                               choices = c(Preview = "head",
                                                           All = "all"),
                                               selected = "head"),
                                  br(),
                                  h5('Download the Delta Ct Summary Table'),
                                  downloadButton("dct_table_DL", "Download"),
                                  br(),
                                  h5('Download the Delta-Delta Ct Summary Table'),
                                  downloadButton("ddct_table_DL", "Download"),
                                  br(),
                                  br()
                              ),
                              mainPanel(
                                  tableOutput('tablehead'),
                                  br(),
                                  tableOutput('dct_table'),
                                  br(),
                                  tableOutput('ddct_table'),
                                  br()
                              ),
                          )
                 ),
                 ## Raw Data Tab ----
                 tabPanel('Raw Data Output',
                          titlePanel("Raw Data Table"),
                          sidebarPanel(
                              helpText('Your data will display here once a suitable file has been uploaded.'),
                              radioButtons("dispout", "Display",
                                           choices = c(Preview = "head",
                                                       All = "all"),
                                           selected = "head")),
                          mainPanel(
                              tableOutput('table'))),
                 ## Graph Tab ----
                 tabPanel('Graph Outputs',
                          titlePanel('Delta Ct and Delta-Delta Ct Graphs'),
                          sidebarPanel(
                              helpText('Once you load a data file and select appropriate columns a Transcript Abundance and a Fold Change graph should appear.'),
                              br(),
                              helpText('Download Transcript Abundance Graph'),
                              downloadButton("downloadPlot1", "Download"),
                              helpText('Download Fold Change Graph'),
                              downloadButton("downloadPlot2", "Download")),
                          mainPanel(plotOutput('dctplot'),
                                    plotOutput('ddctplot'),
                                    br(),
                                    br(),
                                    br(),
                                    br(),
                                    br())
                 ),
                 tags$footer("Created by Richard Browne and Charlotte Francois, La Trobe University, 2019. Last updated 09/06/2019.", align = "right", style = "
position:fixed;
min-height: 5vh;
bottom:0;
width:100%;
height:30px;
color: grey;
padding: 10px;
background-color: black;
z-index: 1000;"
                 )
                 
)
