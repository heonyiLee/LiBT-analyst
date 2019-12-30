library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyalert)

#Edit 2019.12.26
# Define UI for app that draws a histogram ----
ui <- function(request) {shinyUI(
  dashboardPagePlus(
    dashboardHeaderPlus(
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "sliders"
    ),
    sidebar = dashboardSidebar(
      sidebarMenu(
        menuItem("TimeLine", selected = TRUE, icon=icon("clock"),
                 tabName = "timeline", uiOutput("timeline"))
      )
    ), # End of dashboardSidebar
    body = dashboardBody(
      box(
        id = "summary_box",
        solidHeader = T,
        width = 10,
        # tableOutput("uploaded_file_header")
        DT::dataTableOutput("uploaded_file_header")
      ) # End of uploaded file data table
      
    ), # End of dashboardBody
    rightsidebar = rightSidebar(
      useShinyalert(),
      id = "rightsidebar",
      width = 350,
      background = "dark",
      rightSidebarTabContent(
        id=1,
        active = T,
        icon="desktop",
        
        gradientBox(
          title = "Upload Data",
          width = 12,
          gradientColor = "maroon", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            fileInput("fileBrowser", "Choose txt File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            tags$hr(),
            radioButtons("file_type", label="", 
                         choices = list("LFQ" = "LFQ", "iBAQ" = "iBAQ", "TMT" = "TMT"), 
                         selected = "LFQ"),
            
            checkboxGroupInput("first_filtering", label="First Filtering (Remove all)", 
                               choices = list("Potential contaminant" = 1, 
                                              "Reverse" = 2, "Only identified by site" = 3,
                                              "Without Filtering" = 4)) 
            
            
          )
        ), # End of Upload Data box
        actionButton("file_upload_btn", "Upload file")
      ),
      rightSidebarTabContent(
        id=2,
        title = "Input File Format2",
        icon="clock",
        textInput("caption", "Caption", "Data Summary")
      ),
      rightSidebarTabContent(
        id=3,
        title = "Input File Format3",
        icon="gears",
        textInput("caption", "Caption", "Data Summary")
      )
    ), # End of rightSidebar
    title = "DemoWEB"
  )
)}
  