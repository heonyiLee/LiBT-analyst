library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyalert)
library(shinycssloaders)
library(shinyjqui)
library(sortable)
library(V8)
library(httr)

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
                 tabName = "analysis", uiOutput("timeline")),
        menuItem("Analysis", selected = FALSE, icon=icon("chart-pie"),tabName = "analysis"),
        menuItem("User Guide", selected = FALSE, icon=icon("info"),tabName = "user_gide",
                 menuSubItem("Tutorial", tabName = "tutorial"))
      )
    ), # End of dashboardSidebar
    body = dashboardBody(
      tags$head(tags$link(rel="stylesheet",type="text/css",href="body.css"),
                tags$link(rel="stylesheet",type="text/css",href="rightsidebar.css"),
                tags$link(rel="stylesheet",type="text/css",href="timeline.css")),
      tags$head(tags$script(src="body.js")),
      tabItems(
        tabItem(
          tabName = "analysis",
          fluidRow(
            hidden(downloadButton("download_exp_btn","Save Expression Data", "download_exp_btn")),
            hidden(downloadButton("download_dep_info_btn","Save DEP info Data","download_dep_info_btn")),
            box(
              id="data_table",
              solidHeader = T,
              width = 12,
              withSpinner(DT::dataTableOutput("uploaded_file_header"))
            ), # End of uploaded file data table
            tabBox(
              id="plot_tabBox",
              width = 12,
              tabPanel("Results of DEA",
                       fluidRow(
                         # box(
                         #   solidHeader = T,
                         #   width = 4,
                         #   withSpinner(DT::dataTableOutput("result_table_header"))
                         # ),
                         box(
                           id = "pca_plot_box",
                           solidHeader = T,
                           width = 6,
                           prettySwitch(inputId = "show_sampleID", label = "Show SampleID", value = F, width=2, status = "success", inline=T),
                           plotOutput("pca_plot"),
                           downloadButton("download_pca", "Save_png")
                         ),
                         box(
                           id = "correlation_matrix_box",
                           solidHeader = T,
                           width = 6,
                           br(),br(),
                           plotOutput("correlation_matrix"),
                           downloadButton("download_correlation", "Save_png")
                         ),
                         br(),
                         box(
                           id = "volcano_box",
                           solidHeader = T,
                           width = 6,
                           plotOutput("volcano_plot", brush = "volcano_brush"),
                           DT::dataTableOutput("volcano_info"),
                           downloadButton("download_volcano", "Save_png","download_volcano")
                         ),
                         box(
                           id = "heatmap_box",
                           solidHeader = T,
                           width = 6,
                           plotOutput("heatmap"),
                           downloadButton("download_heatmap", "Save_png"),
                           downloadButton("download_gene_cluster","Save_Gene_cluster_info")
                         )
                         ,useShinyalert()
                       )
              ),#tab1 end
              tabPanel("Results of GSA",
                       fluidRow(
                         box(
                           id = "gobp_gsa_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gobp_gsa_plot")),
                           downloadButton("download_gsa_gobp", "Save_GSA_GOBP_all")
                         ),
                         box(
                           id = "gocc_gsa_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gocc_gsa_plot")),
                           downloadButton("download_gsa_gocc", "Save_GSA_GOCC_all")
                         ),
                         box(
                           id = "gomf_gsa_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gomf_gsa_plot")),
                           downloadButton("download_gsa_gomf", "Save_GSA_GOMF_all")
                         ),
                         box(
                           id = "kegg_gsa_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("kegg_gsa_plot")),
                           downloadButton("download_gsa_kegg", "Save_GSA_Kegg_all")
                         )
                         ,useShinyalert()
                       )
              ),#tab2 end
              tabPanel("Results of GSA Pathview",
                       fluidRow(
                         box(
                           title="Top 10 of KEGG pathway",
                           id = "top10_pathview_box",
                           solidHeader = T,
                           width = 12,
                           DT::dataTableOutput("topOfKeggDT")
                         ),
                         box(
                           id = "pathview_select_box",
                           solidHeader = T,
                           width=6,
                           selectInput("pathID_selector", label="Choose a pathway ID",
                                       choices=""),
                           actionButton("render_pathway_btn", "Render pathway"),
                           actionButton("zoom_pathway_btn", "Zoom"),
                           downloadButton("download_pathview", "Save pathway")
                           
                         ),
                         box(
                           id = "pathview_result_box",
                           solidHeader = T,
                           width = 6,
                           plotOutput("pathview_result")
                         )
                         ,useShinyalert()
                       )
              ),#tab3 end
              tabPanel("Results of GSEA",
                       fluidRow(
                         box(
                           id = "gobp_gsea_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gobp_gsea_plot")),
                           downloadButton("download_gsea_gobp", "Save_GSEA_GOBP")
                         ),
                         box(
                           id = "gocc_gsea_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gocc_gsea_plot")),
                           downloadButton("download_gsea_gocc", "Save_GSEA_GOCC")
                         ),
                         box(
                           id = "gomf_gsea_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("gomf_gsea_plot")),
                           downloadButton("download_gsea_gomf", "Save_GSEA_GOMF")
                         ),
                         box(
                           id = "kegg_gsea_box",
                           solidHeader = T,
                           width = 12,
                           withSpinner(plotOutput("kegg_gsea_plot")),
                           downloadButton("download_gsea_kegg", "Save_GSEA_Kegg")
                         )
                         ,useShinyalert()
                       )
              ),#End of GSEA
              tabPanel("Results of GSEA Pathview",
                       fluidRow(
                         box(
                           tags$head(
                             tags$script("function gsea_pathview(d){
                                Shiny.onInputChange('js.gsea.pathview',d.getAttribute('val'));}"
                             )
                           ),
                           id = "gsea_pathview_box",
                           solidHeader = T,
                           width = 12,
                           hidden(imageOutput("gsea_pathview_image",height="1000px")),
                           h3(id="gsea_up_title","Up-regulated:"),
                           DT::dataTableOutput("result_of_gsea_up_regulated"),
                           h3(id="gsea_down_title","Down-regulated:"),
                           DT::dataTableOutput("result_of_gsea_down_regulated"),
                           downloadButton("download_gsea_pathview","Save_GSEA_Pathview Zip", "download_gsea_pathview"),
                         ),
                         useShinyalert()
                       )
              ),
              tabPanel("Result of PPI Network Analysis",
                       fluidRow(
                         box(
                           id = "ppi_box",
                           solidHeader = T,
                           width = 12,
                           downloadButton("ppi_image_download_btn", "Download PPI Network PNG"),
                           downloadButton("ppi_tsv_download_btn", "Download PPI Network tsv"),
                           withSpinner(uiOutput("ppi_image"))
                         )
                         ,useShinyalert()
                       )
              )
            )#End of tabBox
          )#End of fluidRow
        ),#End of tabItem_analysis
        tabItem(
          tabName = "tutorial",
          h2("Tutorial")
        )
      )#End fo tabItems_tutorial
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
                         choices = list("T.Test" = "T.Test", "Wilcoxon-Ranksum" = "Wilcoxon-Ranksum", "edgeR" = "edgeR")),
            radioButtons("padj_method", label="Choose p.adj method", 
                         choices = list("Benjamini-Hocherg" = "BH", "Bonferroni" = "bonferroni")),
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
            numericInput("dea_log2fc", label = HTML("Set the threshold <br/>for the log<sub>2</sub> Fold Change"),
                         value=1.5),
            numericInput("dea_clusterNum", label = HTML("Set a number of cluster <br/>for heatmap"), value=""),
            tags$hr(),
            actionButton("dea_btn", "Start Draw")
          )
        ) # End of DEP box
      ),
      rightSidebarTabContent(
        id=4,
        icon="chart-pie",
        active = F,
        gradientBox(
          title = HTML("Gene Set Analysis <br/><br/>"),
          width = 12,
          icon = "chart-pie",
          gradientColor = "yellow",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            uiOutput("dep_up"),
            uiOutput("dep_down"),
            radioButtons("gsa_set", label="Choose Data set",
                         choices = list("Case-UP" = "caseup", "Ctrl-UP" = "casedown")),
            radioButtons("gsa_input_set", label="Choose input Data",
                         choices = list("Total" = "total", "DEP" = "dep")),
            # radioButtons("gsa_tool", label="Choose GSA Tool", 
            #              choices = list("enrichR" = "enrichR", "DAVID" = "DAVID")),
            numericInput("set_nterm", label = HTML("Set a number of showed Term <br/>for barplot"),
                         value=10),
            actionButton("gsa_btn", "Start GSA")
          )
        ),
        gradientBox(
          title = HTML("Gene Set <br/>Enrichment Analysis <br/><br/>"),
          width = 12,
          icon = "chart-pie",
          gradientColor = "yellow",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            tags$div(
              id="select_genelevel_stats", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="select_genelevel_stats", "Choose stats of gene level"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="log2fcxmlog10padj", checked="checked",
                                             tags$span(HTML("log<sub>2</sub>FC*-log<sub>10</sub>(P.adj)")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="P.adj",
                                             tags$span(HTML("P.adj")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="P.value",
                                             tags$span(HTML("P.value")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="log2fc",
                                             tags$span(HTML("log<sub>2</sub>FoldChange")))
                                )
                       )
              )
            ),
            actionButton("gsea_btn", "Start GSEA")
          )
        ), # End of GSEA box
        gradientBox(
          title = HTML("Protein-Protein Interaction <br/>Network Analysis <br/><br/>_StringDB"),
          width = 12,
          icon = "chart-pie",
          gradientColor = "yellow",
          boxToolSize = "md",
          closable = F,
          footer = fluidRow(
            pickerInput("input_organism_toppi", label=HTML("Choose organism"), choices = c("Homo Sapiens", "Mus musculus")),
            numericInput("input_num_toppi", label = HTML("Choose number of input gene <br/> (Defualt : # of DEP)"),value=100),
            tags$div(
              id="select_ppi_condition", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="select_ppi_condition", "Choose condition of gene select"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="P.adj", checked="checked",
                                             tags$span(HTML("P.adj")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="P.value",
                                             tags$span(HTML("P.value")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="log2fc",
                                             tags$span(HTML("log<sub>2</sub>FoldChange")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="Gene_Name",
                                             tags$span(HTML("Gene Name")))
                                )
                       )
              )
            ),
            hidden(textInput("input_gene_toppi", "Input gene name", placeholder = "Ex)TP53 or TP53, PLK1, RHOB")),
            actionButton("ppi_btn", "Start PPIA")
          )
        ) # End of STRING box
      )
    ), # End of rightSidebar
    title = "LiBT-Analyst",
    useShinyjs()
  ))
}


