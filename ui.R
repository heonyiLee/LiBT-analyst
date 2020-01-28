library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyalert)
library(shinycssloaders)

# 2019.12.30
ui <- function(request) {shinyUI(
  dashboardPagePlus(
    dashboardHeaderPlus(
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "gears",
      title = "LiBT-Analyst",
      titleWidth = 300
    ),
    sidebar = dashboardSidebar(
      width=300,
      sidebarMenu(
        menuItem("TimeLine", selected = TRUE, icon=icon("history"),
                 tabName = "timeline", uiOutput("timeline"))
      )
    ), # End of dashboardSidebar
    body = dashboardBody(
      
      tags$head(tags$link(rel="stylesheet",type="text/css",href="body.css"),
                tags$link(rel="stylesheet",type="text/css",href="rightsidebar.css"),
                tags$link(rel="stylesheet",type="text/css",href="timeline.css")),
      tags$head(tags$script(src="body.js")),
      box(
        id = "summary_box",
        solidHeader = T,
        width = 11,
        # tableOutput("uploaded_file_header")
        withSpinner(DT::dataTableOutput("uploaded_file_header"))
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
          useShinyjs(),
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
            
            hidden(radioButtons("TMT_normalize", label="Get Normalized TMT data", 
                         choices = list("YES" = "T", "NO" = "F"), 
                         selected = "T")),
            
            checkboxGroupInput("first_filtering", label="Filtering Options (Remove all)", 
                               choices = list("Potential contaminant" = "potential", 
                                              "Reverse" = "reverse", 
                                              "Only identified by site" = "identified")),
            actionButton("select_all_filtering_btn", "Select All"),
            actionButton("deselect_all_filtering_btn", "Deselect All"),
            tags$hr(),
            actionButton("file_upload_btn", "Upload file")
          )
        ),
        gradientBox(
          title = "Experiment Design",
          width = 12,
          gradientColor = "maroon", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            selectInput("case_group_selection", label="Case samples",
                        choices=c(),selectize=T, multiple=T),
            selectInput("control_group_selection", label="Control samples",
                        choices=c(),selectize=T, multiple=T)
          )
        ) # End of Upload Data box
      ),
      rightSidebarTabContent(
        id=2,
        icon="dna",
        #active = T,
        gradientBox(
          title = "Transformation",
          width = 12,
          gradientColor = "aqua", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            radioButtons("transformation", label="", 
                         choices = list("log2" = "log2", "none" = "none"), 
                         selected = "log2"),
          )
        ),
        gradientBox(
          title = "Filter based on Valid Value",
          width = 12,
          gradientColor = "aqua", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            radioButtons("valid_value", label="Choose % of valid value", 
                         choices = list("30%" = 30, "50%" = 50, "70%"=70, "100%"=100), 
                         selected = 70)
          )
        ),
        gradientBox(
          title = "Deal with Missing Value",
          width = 12,
          gradientColor = "aqua", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            radioButtons("imputation", label="Choose Imputation", 
                         choices = list("Normal distribution" = "nordis", "Constant" = "const", 
                                        "NaN"="nan", "none"="none"), 
                         selected = "nordis")
          )
        ),
        gradientBox(
          title = "Normalization",
          width = 12,
          gradientColor = "aqua", 
          boxToolSize = "md", 
          closable = F,
          footer = fluidRow(
            radioButtons("normalization", label="Choose method of normalization", 
                         choices = list("Width distribution" = "quantile", "none" = "none"), 
                         selected = "quantile")
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
    title = "LiBT-Analyst"
  )
)}