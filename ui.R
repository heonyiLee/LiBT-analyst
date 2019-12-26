library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
#Edit 2019.12.26
# Define UI for app that draws a histogram ----
shinyUI(dashboardPagePlus(
  dashboardHeaderPlus(
    enable_rightsidebar = TRUE,
    rightSidebarIcon = "clock"
  ),
  dashboardSidebar(
    br(),
    gradientBox(
      title = "Input Data Type",
      width = 12,
      #icon = "fa fa-heart",
      gradientColor = "maroon", 
      boxToolSize = "md", 
      closable = F,
      footer= radioButtons("data_type", label="", 
                           choices = list("LFQ" = "LFQ", "iBAQ" = "iBAQ", "TMT" = "TMT"), 
                           selected = "LFQ")
    ),
    gradientBox(
      title = "Input CSV Data",
      width = 12,
      #icon = "fa fa-heart",
      gradientColor = "maroon", 
      boxToolSize = "md", 
      closable = F,
      footer = fluidRow(
        fileInput("input_csv", "Choose CSV File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),
        # Horizontal line ----
        tags$hr(),
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),
        # Input: Select separator ----
        radioButtons("sep", "Separator",
                     choices = c(Comma = ",",
                                 Semicolon = ";",
                                 Tab = "\t"),
                     selected = ","),
        # Input: Select quotes ----
        radioButtons("quote", "Quote",
                     choices = c(None = "",
                                 "Double Quote" = '"',
                                 "Single Quote" = "'"),
                     selected = '"'),
        # Horizontal line ----
        tags$hr(),
        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head"),
        tags$hr(),
        checkboxGroupInput("first_filtering", label="First Filtering (Remove all)", 
                     choices = list("Potential contaminant" = 1, "Reverse" = 2, "Only identified by site" = 3), 
                     selected = T)
      )
    )
  ),
  dashboardBody(
    box(
      width = NULL,
      tableOutput("input_data")
    )
  ),
  rightSidebar(
    background = "dark",
    rightSidebarTabContent(
      id=1,
      active = T,
      icon="",
      uiOutput("timeline")
    )
  ),
  title = "DemoWEB"
  )
)