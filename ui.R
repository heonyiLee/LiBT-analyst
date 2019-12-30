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
      tags$head(tags$link(rel="stylesheet",type="text/css",href="body.css")),
      tags$head(tags$link(rel="stylesheet",type="text/css",href="rightsidebar.css")),
      tags$head(tags$link(rel="stylesheet",type="text/css",href="timeline.css")),
      box(
        id = "summary_box",
        solidHeader = T,
        width = 11,
        # tableOutput("uploaded_file_header")
        DT::dataTableOutput("uploaded_file_header")
      ) # End of uploaded file data table
      ,useShinyalert(),
    ), # End of dashboardBody
    rightsidebar = rightSidebar(
      id = "rightsidebar",
      width = 350,
      background = "dark",
      rightSidebarTabContent(
        id=1,
        active = T,
        icon="file-upload",
        
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
        ),
        gradientBox(
            title = "Grouping Data",
            width = 12,
            gradientColor = "maroon", 
            boxToolSize = "md", 
            closable = F,
            footer = ""
        ) # End of Upload Data box
        ,actionButton("file_upload_btn", "Upload file")
      ),
      rightSidebarTabContent(
        id=2,
        icon="gears",
        #active = T,
        gradientBox(
          title = "Transformation",
          width = 12,
          gradientColor = "teal", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            radioButtons("transformation", label="", 
                         choices = list("log2" = "log2", "none" = "none"), 
                         selected = "log2"),
          )
        )
        ,actionButton("preprocess_btn", "Preprocess")
      ),
      rightSidebarTabContent(
        id=3,
        icon="chart-bar",
        textInput("caption", "Caption", "Data Summary")
      )
    ), # End of rightSidebar
    title = "DemoWEB"
  )
)}
  