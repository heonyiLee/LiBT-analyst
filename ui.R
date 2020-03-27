library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyalert)
library(shinycssloaders)
library(shinyjqui)
library(sortable)


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
        id="data_table",
        solidHeader = T,
        width = 12,
        withSpinner(DT::dataTableOutput("uploaded_file_header"))
      ), # End of uploaded file data table
      # box(
      #   solidHeader = T,
      #   width = 4,
      #   withSpinner(DT::dataTableOutput("result_table_header"))
      # ),
      box(
        id = "pca_plot_box",
        solidHeader = T,
        width = 6,
        plotOutput("pca_plot"),
        downloadButton("download_pca", "Save_png")
      ),
      box(
        id = "volcano_box",
        solidHeader = T,
        width = 6,
        plotOutput("volcano_plot"),
        downloadButton("download_volcano", "Save_png")
      ),
      box(
        id = "correlation_matrix_box",
        solidHeader = T,
        width = 6,
        plotOutput("correlation_matrix"),
        downloadButton("download_correlation", "Save_png")
      ),
      box(
        id = "heatmap_box",
        solidHeader = T,
        width = 6,
        plotOutput("heatmap"),
        downloadButton("download_heatmap", "Save_png")
      )
      ,useShinyalert(),
    ), # End of dashboardBody
    rightsidebar = rightSidebar(
      id = "entireRightsidebar",
      width = 350,
      background = "dark",
      rightSidebarTabContent(
        uiOutput("sidebar_tab1"),
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
            radioButtons("file_type", label="", 
                         choices = list("LFQ" = "LFQ", "iBAQ" = "iBAQ", "TMT" = "TMT"), 
                         selected = "LFQ"),
            tags$hr(),
            
            fileInput("fileBrowser", "Choose txt File",
                      multiple = FALSE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            hidden(radioButtons("TMT_input_option", label="Get Normalized TMT data", 
                                choices = list("YES" = "T", "NO" = "F"), 
                                selected = "F")),
            
            checkboxGroupInput("nonTMT_input_option", label="Filtering Options (Remove all)", 
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
            pickerInput(inputId = "case_group_selection", label="Case samples", multiple=TRUE,
                        choices=c(), options=list(`actions-box` = TRUE, `selected-text-format` = "count > 1")),
            pickerInput(inputId = "control_group_selection", label="Control samples", multiple=TRUE,
                        choices=c(), options=list(`actions-box` = TRUE, `selected-text-format` = "count > 1")),
            tags$hr(),
            actionButton("exp_design_submit_btn", "Submit")
          )
        ) # End of Upload Data box
      ),
      rightSidebarTabContent(
        uiOutput("sidebar_tab2"),
        id="sidebar_tab2",
        icon="dna",
        active = F,
        gradientBox(
          title = "Details",
          width = 12,
          gradientColor = "aqua",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            hidden(prettyToggle(inputId="exp_design_check", 
                                label_on="submitted", label_off="not_submitted", value=FALSE)),
            bucket_list(
              header = "Drag to select options and set the order",
              add_rank_list(
                text = "Use",
                labels = list(
                  "Use_valid_value" = radioButtons("valid_value", label="Choose % of valid value",
                                                   choices = list("30%" = 0.3, "50%" = 0.5, "70%"= 0.7, "100%"=0)),
                  "Use_imputation" = radioButtons("imputation", label="Choose Imputation",
                                                  choiceNames = list(HTML("QRILC<br/>(quantile regression-based<br/>left-censored function)"),
                                                                     HTML("MinProb<br/>(left-shifted Gaussian distribution)"), 
                                                                     HTML("Man<br/>(manually defined<br/>left-shifted Gaussian distribution)"),
                                                                     HTML("Constant<br/>(replace with zero value)")),
                                                  choiceValues = list("QRILC", "MinProb", "man", "zero")),
                  # choices = list("Normal distribution" = "normal_distribution",
                  #                "Constant" = "constant",
                  #                "NaN"="nan")),
                  "Use_normalization" = radioButtons("normalization", label="Choose method of normalization",
                                                     choices = list("Yes" = "Yes",
                                                                    "No" = "No"))
                  # choices = list("Width distribution" = "quantile",
                  #                "Z-score"="zscore", "none" = "none"))
                ),
                input_id = "use_options"
              ),
              add_rank_list(
                text = "Not to Use",
                labels = NULL,
                input_id = "not_to_use_options"
              )
            ),
            tags$hr(),
            actionButton("preprocess_btn", "Start preprocessing")
          ) # End of Details box's footer
        ) # End of Details box
      ),
      rightSidebarTabContent(
        id=3,
        icon="chart-bar",
        active = F,
        gradientBox(
          title = HTML("Differential <br/>experimental analysis <br/><br/>_Test"),
          width = 12,
          icon = "chart-bar",
          gradientColor = "green",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            uiOutput("dea_case"),
            uiOutput("dea_control"),
            radioButtons("test_method", label="Choose test method", 
                         choices = list("T.Test" = "T.Test", "Ranksum-Wilcoxon" = "Ranksum-Wilcoxon", "edgeR" = "edgeR", "Limma" = "Limma")),
            radioButtons("padj_method", label="Choose p.adj method", 
                         choices = list("FDR" = "fdr", "Bonferroni" = "bonferroni")),
            actionButton("test_btn", "Start DEA")
          )
        ),
        gradientBox(
          title = HTML("Differential <br/>experimental analysis <br/><br/>_Visulization"),
          width = 12,
          icon = "chart-bar",
          gradientColor = "green",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            radioButtons("thres_type", label="Choose threshold type", 
                         choices = list("P.adj" = "P.adj", "P.value" = "P.value", "None" = "none")),
            numericInput("dea_pvalue", label = HTML("Set the threshold <br/>for the adjusted P-value or P-value"),
                         value=0.05),
            numericInput("dea_log2fc", label = HTML("Set the threshold <br/>for the log2 fold change"),
                         value=1.5),
            numericInput("dea_clusterNum", label = HTML("Set a number of cluster <br/>for heatmap"),
                         value=3),
            tags$hr(),
            actionButton("dea_btn", "Start Draw")
          )
        ) # End of DEP box
      )
    ), # End of rightSidebar
    title = "LiBT-Analyst",
    useShinyjs()
  ))
}
