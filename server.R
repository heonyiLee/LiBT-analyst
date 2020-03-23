library(dplyr)
library(readr)
library(stringr)
library(DEP)
library(sortable)


source("function.R")
# 2019.12.30

shinyServer(function(input,output, session){
  
  options(shiny.maxRequestSize=100*1024^2)
  
  step <- c()
  info <- c()
  sample_num <- c()
  time <- c()
  
  ############################
  ##      RENDERING         ##
  ############################
  
  observeEvent(input$file_type, {
    if(input$file_type=="TMT"){
      shinyjs::show("TMT_input_option")
      shinyjs::hide("nonTMT_input_option")
      shinyjs::hide("select_all_filtering_btn")
      shinyjs::hide("deselect_all_filtering_btn")
    }else{
      shinyjs::hide("TMT_input_option")
      shinyjs::show("nonTMT_input_option")
      shinyjs::show("select_all_filtering_btn")
      shinyjs::show("deselect_all_filtering_btn")
    }
  })
  
  
  observeEvent(input$fileBrowser, {
    temp <- file_input()
    if(is.null(temp)){
      shinyalert("Check your file type!", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      reset("fileBrowser")
    } else {
      shinyalert("Data will be converted to log2 scale in the following order!", type = "info")
      if(input$file_type=="TMT"){
        state <- colnames(temp)[grep("normalized", colnames(temp), ignore.case = T)]
        updateRadioButtons(session, "TMT_input_option", label="Get Normalized TMT data",
                           choices = list("YES" = "T", "NO" = "F"), selected="F")
        if(length(state)==0){
          shinyjs::disable("TMT_input_option")
        }
      }
    }
  })
  
  observeEvent(input$select_all_filtering_btn, {
    if(input$select_all_filtering_btn == 0) {
      return(NULL)
    } else if(input$select_all_filtering_btn >0){
      updateCheckboxGroupInput(session, "nonTMT_input_option",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c("potential","reverse","identified"))
    }
  })
  
  observeEvent(input$deselect_all_filtering_btn, {
    if(input$deselect_all_filtering_btn == 0) {
      return(NULL)
    } else if(input$deselect_all_filtering_btn >0){
      updateCheckboxGroupInput(session, "nonTMT_input_option",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c())
    }
  })
  
  observeEvent(input$file_upload_btn, {
    if((is.null(input$fileBrowser) || 
        is.null(input$nonTMT_input_option)) && length(input$TMT_input_option)<0){
      shinyalert("Choose option!", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    } else {
      render_df <- main_data()
      if(length(render_df) != 0){
        if(input$file_type != "TMT"){
          info <- paste0("File Type : ", input$file_type,"\n",
                         "Option : ", input$nonTMT_input_option)
          timeLine <<- data.frame(step="Data Input",info=info,
                                  sample_num=as.numeric(nrow(render_df)),
                                  time=as.character(Sys.time()),color="maroon")
          addTimeLine(timeLine)
        } else{
          info <- paste0("File Type : ", input$file_type,"\n",
                         "Is normalized ? : ", input$TMT_input_option)
          timeLine <<- data.frame(step="Data Input",info=info,
                                  sample_num=as.numeric(nrow(render_df)),
                                  time=as.character(Sys.time()),color="maroon")
          addTimeLine(timeLine)
        }
        output$uploaded_file_header <- DT::renderDataTable({
          render_df}, 
          options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
        
        samples <- total_samples()
        print("!!!!!!!!")
        print(samples)
        updatePickerInput(session, "case_group_selection", choices = samples)
      }
    }
  })
  
  
  observeEvent(input$case_group_selection, {
    choices <- control_samples()
    updatePickerInput(session, "control_group_selection",
                      choices = choices, selected = choices)
  })
  
  observeEvent(input$exp_design_submit_btn, {
    updatePickerInput(session, "control_group_selection",
                      label=paste0("Control samples (n=", length(control_samples()), ")"),
                      choices = control_samples(), selected=control_samples())
    updatePickerInput(session, "case_group_selection",
                      label=paste0("Case samples (n=", length(case_samples()), ")"),
                      choices = case_samples(), selected=case_samples())
    
    updatePrettyToggle(session, "exp_design_check", label=NULL, value=TRUE)
    
    
    tmp <- paste0("# of Case : ", length(case_samples()),"\n",
                  "# of Control : ", length(control_samples()))
    info <- paste0(info,tmp,"\n")
    newTL <- data.frame(step="Exp Design Submitted",
                        info=info,
                        sample_num=as.numeric(nrow(main_data())),
                        time=as.character(Sys.time()),color="maroon")
    timeLine <- rbind(timeLine,newTL)
    addTimeLine(timeLine)
  })
  
  observeEvent(input$preprocess_btn, {
    if(input$exp_design_check==TRUE) {
      data_se <- ready_for_dea()
      output$uploaded_file_header <- DT::renderDataTable({
        assay(data_se)}, 
        options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
      
      tmp <- paste0("Valid value : ", paste0(input$valid_value,"%"), "\n",
                    "Imputation : " , str_to_title(input$imputation), "\n",
                    "Normalization : ", str_to_title(input$normalization))
      info <- paste0(info,tmp,"\n")
      newTL <- data.frame(step="Preprocessing",
                          info=info,
                          sample_num=as.numeric(nrow(data_se)),
                          time=as.character(Sys.time()),color="maroon")
      timeLine <- rbind(timeLine, newTL)
      addTimeLine(timeLine)
      
      
      output$dea_case <- renderUI({
        boxPlus(
          title = "Case samples",
          closable = F,
          collapsed = T,
          enable_label = T,
          label_text = length(case_samples()),
          label_status = "danger",
          width = 12,
          solidHeader = F,
          collapsible = T,
          footer = HTML(paste(case_samples(), collapse=",<br/>"))
        )
      })
      
      output$dea_control <- renderUI({
        boxPlus(
          title = "Control samples",
          closable = F,
          collapsed = T,
          enable_label = T,
          label_text = length(control_samples()),
          label_status = "danger",
          width = 12,
          solidHeader = F,
          collapsible = T,
          footer = HTML(paste(control_samples(), collapse=",<br/>"))
        )
      })
      
    } else {
      shinyalert("Please submit experiment design", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
    
  })
  
  
  observeEvent(input$dea_btn, {
    if(!is.null(ready_for_dea())){
      output$pca_plot <- renderPlot({
        pca_input()
      })
      
      output$correlation_matrix <- renderPlot({
        correlation_input()
      })
      
      output$heatmap <- renderPlot({
        heatmap_input()
      })
      
      output$volcano_plot <- renderPlot({
        volcano_input()
      })
      
    }
  })
  
  
  
  
  
  
  
  # output$download_frequency_svg <- downloadHandler(
  #   filename = function() {"frequency_plot.png"},
  #   content = function(file) {
  #     png(file)
  #     print(plot_imputation(data_norm(), data_imp()))
  #     dev.off()
  #   }
  # )
  
  
  
  ##--------------------------------------------------- reactive/EventReactive Section
  file_input <- reactive({NULL})
  file_input <- eventReactive(input$fileBrowser, {
    req(input$fileBrowser)
    
    if(is.null(input$fileBrowser)){
      return(NULL)
    }
    temp_df <- readLines(input$fileBrowser$datapath, n=1)
    
    if(grepl("\t", temp_df)){
      sep <- c("\t")
    } else if(grepl(";", temp_df)) {
      sep <- c(";")
    } else if(grepl(",", temp_df)) {
      sep <- c(",")
    } else {
      sep <- c(" ")
    }
    
    temp_df <- read.table(input$fileBrowser$datapath,
                          header = T, fill = T,
                          sep = sep)
    state <- file_input_test(temp_df, input$file_type)
    if(length(state)==0){
      return(NULL)
    } else {
      return(temp_df)
    }
  })
  
  
  
  main_data <- reactive({NULL})
  main_data <- eventReactive(input$file_upload_btn, {
    temp_df <- file_input()
    
    file_type <- input$file_type
    if(file_type=="TMT"){
      temp_df <- get_main_data_T(temp_df,input$TMT_input_option)
  
    } else{
      checked_option <- input$nonTMT_input_option
      temp_df <- filter_with_option(checked_option, temp_df)
      temp_df <- get_main_data_LiB(temp_df, file_type)
    }
  })
  
  
  total_samples <- reactive({
    df <- main_data()
    print("###########")
    print(colnames(df))  
    file_type <- input$file_type
    if(file_type=="TMT"){
      if(input$TMT_input_option=="T"){
        samples <- make_case_samples_diffT(df, input$TMT_input_option)
      } else {
        samples <- make_case_samples_T(df)
      }
    } else{
      samples <- make_case_samples_LiB(df, file_type)
    }
    return(samples)
  })
  
  case_samples <- reactive({
    case_samples <- input$case_group_selection
    
    return(case_samples)
  })
  
  control_samples <- eventReactive(input$case_group_selection, {
    control_samples <- setdiff(total_samples(), case_samples())
    
    return(control_samples)
  })
  
  # exp_DesignData <- reactive({
  #   case <- case_samples()
  #   ctrl <- control_samples()
  #   exp_design <- make_expDesignData(case,ctrl)
  #   return(exp_design)
  # })
  
  summarized_Data <- reactive({
    case <- case_samples()
    ctrl <- control_samples()
    design <- make_expDesignData(case, ctrl)
    main <- main_data()
    file_type <- input$file_type
    summary <- make_summarizedData(main,file_type,design)
    return(summary)
  })
  
  
  preprocessed_data <- reactive({
    file_type <- input$file_type
    transformation_item <- input$transformation
    filter_item <- input$valid_value
    imputation_item <- input$imputation
    normalization_item <- input$normalization
    
    case <- case_samples()
    control <- control_samples()
    
    data_se <- summarized_Data()
    data_filt <- NULL
    data_norm <- NULL
    data_imp <- NULL
    
    preprocessing_options <- input$use_options
    for(i in 1:length(preprocessing_options)) {
      switch(preprocessing_options[i],
             "Use_valid_value" = {
               data_se <- use_valid_option(data_se, case, control, input$valid_value)
               data_filt <- data_se},
             "Use_imputation" = {
               data_se <- use_imputation_option(data_se, case, control, input$imputation)
               data_imp <- data_se},
             "Use_normalization" = {
               data_se <- use_normalization_option(data_se,input$normalization)
               data_norm <- data_se},
             NULL = {
               data_se <- data_se
             }
      )
    }
    preprocessed_data <- list(data_se, data_filt, data_norm, data_imp)
    return (preprocessed_data) # format : SummarizedExperiment
  })
  
  ready_for_dea <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[1]]
  })
  
  data_filt <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[2]]
  })
  
  data_norm <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[3]]
  })
  
  data_imp <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[4]]
  })
  
  
  dep <- eventReactive(input$dea_btn,{
    req(preprocessed_data())
    
    pvalue <- input$dea_pvalue
    log2fc <- input$dea_log2fc
    
    data_diff <- test_diff(ready_for_dea(), type="control", control="control")
    dep <- add_rejections(data_diff, alpha=pvalue, lfc=log2fc)
  })
  
  pca_input <- reactive({
    plot_pca(dep(), x=1, y=2, n=500, point_size=3, indicate = c("condition"))
  })
  
  correlation_input <- reactive({
    plot_cor(dep(), significant = T, pal="Spectral", lower=-1, upper=1, indicate=c("condition"))
  })
  
  heatmap_input <- reactive({
    col_limit <- length(total_samples())
    plot_heatmap(dep(), type="centered", kmeans=T, k=6,
                 col_limit=4, show_row_names=F, indicate=c("condition"))
  })
  
  volcano_input <- reactive({
    plot_volcano(dep(), contrast="case_vs_control", label_size=2, add_names=T)
  })
  
  data_results <- reactive({
    data_results <- get_results(dep())
  })
  
  
  addTimeLine <-  function(timeLine){
    output$timeline <- renderUI({
      timelineBlock(
        reversed = F,
        timelineEnd(color = "gray"),
        lapply(1:nrow(timeLine), FUN = function(i){
          tagList(
            timelineItem(
              icon = "file-upload",
              color = timeLine$color[i],
              time = timeLine$time[i],
              tags$h4(timeLine$step[i]),
              footer= timeLine$info[i]
              
            )
            ,timelineLabel(paste0("# of proteins : ",timeLine$sample_num[i]), color = "olive")
          )
        }),
        timelineStart(color = "gray")
      )
    })
  } # End of addTimeLine
})