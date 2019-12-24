library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
# Define UI for app that draws a histogram ----
shinyUI(dashboardPage(
  dashboardHeader(
    
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
      footer= radioButtons("radio", label="", 
                           choices = list("LFQ" = 1, "iBAQ" = 2, "TMT" = 3), 
                           selected = 1)
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
      collapsible = TRUE,
      title = "Reversed Timeline: ideal to include in a box",
      status = "info",
      width = NULL,
      timelineBlock(
        timelineEnd(color = "danger"),
        timelineLabel(2018, color = "teal"),
        timelineItem(
          title = "Item 1",
          icon = "gears",
          color = "olive",
          time = "now",
          footer = "Here is the footer",
          "This is the body"
        ),
        timelineItem(
          title = "Item 2",
          border = FALSE
        ),
        timelineLabel(2015, color = "orange"),
        timelineItem(
          title = "Item 3",
          icon = "paint-brush",
          color = "primary",
          timelineItemMedia(src = "http://placehold.it/150x100"),
          timelineItemMedia(src = "http://placehold.it/150x100")
        ),
        timelineStart(color = "gray")
      )
    ),
    box(
      width = NULL,
      tableOutput("input")
    )
  ), 
  title = "DemoWEB"
  )
)