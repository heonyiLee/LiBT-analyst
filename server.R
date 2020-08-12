library(dplyr)
library(readr)
library(stringr)
library(DEP)
library(sortable)
library(SummarizedExperiment)
library(ggplot2)
library(edgeR)
library(enrichR)
library(tibble)
library(phenoTest)
library(pathview)

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
      
      if(input$file_type=="TMT"){
        state <- colnames(temp)[grep("normalized", colnames(temp), ignore.case = T)]
        updateRadioButtons(session, "TMT_input_option", label="Get Normalized TMT data",
                           choices = list("YES" = "T", "NO" = "F"), selected="F")
        if(length(state)==0){
          shinyjs::disable("TMT_input_option")
        }
      }
      
      info = paste0("* File Type : ", input$file_type,"\n")
      timeLine <<- data.frame(step="Start!:)",info=info,
                              sample_num=as.numeric(nrow(temp)),
                              time=as.character(Sys.time()),color="maroon",icon="file-upload")
      addTimeLine(timeLine)
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
      shinyalert("Data will be converted to log2 scale in the following order!", type = "info")
      render_df <- main_data()
      if(length(render_df) != 0){
        if(input$file_type != "TMT"){
          option <- c()
          for(i in 1:length(input$nonTMT_input_option)){
            if(i!=length(input$nonTMT_input_option)){
              option <- paste0(option,input$nonTMT_input_option[i]," / ")  
            } else{
              option <- paste0(option,input$nonTMT_input_option[i])  
            }
          }
          info <- paste0("* Numerical Filter\n  : Peptides = 0\n    Intensity = 0\n",
                         "* Categorical Filter\n", "  : ", option)
          newTL <- data.frame(step="Data Input",info=info,
                              sample_num=as.numeric(nrow(render_df)),
                              time=as.character(Sys.time()),color="maroon",icon="file-upload")
          timeLine <<- rbind(timeLine,newTL)
          addTimeLine(timeLine)
        } else{
          info <- paste0("* Is normalized ? : ", input$TMT_input_option)
          newTL <- data.frame(step="Data Input",info=info,
                              sample_num=as.numeric(nrow(render_df)),
                              time=as.character(Sys.time()),color="maroon",icon="file-upload")
          timeLine <<- rbind(timeLine,newTL)
          addTimeLine(timeLine)
        }
        output$uploaded_file_header <- DT::renderDataTable({
          render_df}, 
          options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
        
        samples <- total_samples()
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
    
    
    condition <- make_condition(case_samples(),control_samples())
    group_name <- unique(condition)
    tmp <- paste0("* Case Group : ", group_name[1], "\n",
                  "* Control Group : ", group_name[2], "\n",
                  "* # of Case : ", length(case_samples()),"\n",
                  "* # of Control : ", length(control_samples()))
    info <- paste0(info,tmp,"\n")
    newTL <- data.frame(step="Exp Design Submitted",
                        info=info,
                        sample_num=as.numeric(nrow(main_data())),
                        time=as.character(Sys.time()),color="maroon",icon="file-upload")
    timeLine <<- rbind(timeLine,newTL)
    addTimeLine(timeLine)
  })
  
  observeEvent(input$preprocess_btn, {
    if(input$exp_design_check==TRUE) {
      data_se <- ready_for_dea()
      # res <- round(assay(data_se),digits=2)
      res <- assay(data_se)
      output$uploaded_file_header <- DT::renderDataTable({
        DT::datatable(res,options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), selection="none") %>% 
          DT::formatRound(colnames(res), digits=2)
      }, server=F) 
      
      tmp <- c()
      preprocessing_options <- input$use_options
      for(i in 1:length(preprocessing_options)) {
        if(i != 1){
          tmp <- paste0(tmp,"\n")
        }
        switch(preprocessing_options[i],
               "Use_valid_value" = {
                 vv <- as.character(as.numeric(input$valid_value)*100)
                 tmp <- paste0(tmp,i,". Valid value : ", paste0(vv,"%"))},
               "Use_imputation" = {
                 tmp <- paste0(tmp,i,". Imputation : " , str_to_title(input$imputation))},
               "Use_normalization" = {
                 tmp <- paste0(tmp,i,". Normalization : ", str_to_title(input$normalization))},
               NULL = {
                 tmp <- tmp
               }
        )
      }
      
      info <- paste0(info,tmp,"\n")
      newTL <- data.frame(step="Preprocessing",
                          info=info,
                          sample_num=as.numeric(nrow(data_se)),
                          time=as.character(Sys.time()),color="aqua",icon="dna")
      timeLine <<- rbind(timeLine, newTL)
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
  
  observeEvent(input$test_btn, {
    if(!is.null(res_test())){
      shinyalert("Complete Tests!", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      info <- paste0("* Test Method : ", input$test_method,"\n",
                     "* P.adj Method : ", input$padj_method)
      newTL <- data.frame(step="DEA_Test",info=info,
                          sample_num=as.numeric(nrow(assay(res_test()))),
                          time=as.character(Sys.time()),color="green",icon="chart-bar")
      timeLine <<- rbind(timeLine,newTL)
      addTimeLine(timeLine)
    }
  })
  
  observeEvent(input$dea_btn, {
    print(dep())
    withProgress(message = 'Plots calculations are in progress',
                 detail = 'Please wait for a while', value = 0, {
                   for (i in 1:10) {
                     incProgress(1/10)
                     Sys.sleep(0.25)
                   }
                 })
    
    sig <- which(rowData(dep())$significant==T)
    if(input$thres_type == "none"){
      info <- paste0("* Threshold type : ", input$thres_type,"\n",
                     "* The number of cluster : ", input$dea_clusterNum)
    } else{
      info <- paste0("* Threshold type : ", input$thres_type,"\n",
                     "* Threshold value for\n  ", 
                     input$thres_type, " : ", input$dea_pvalue,"\n",
                     "  log2FC : ", input$dea_log2fc,"\n",
                     "* The number of cluster : ", input$dea_clusterNum)
    }
    newTL <- data.frame(step="DEA_Visualization",info=info,
                        sample_num=length(sig),
                        time=as.character(Sys.time()),color="green",icon="chart-bar")
    timeLine <<- rbind(timeLine,newTL)
    addTimeLine(timeLine)
    
    if(!is.null(ready_for_dea())){
      output$volcano_plot <- renderPlot({
        volcano_input()
      })
      output$pca_plot <- renderPlot({
        pca_input_noSample()
      })
    }
    
    if(length(sig)==0){
      shinyalert("Ther is no DEP!", "Change threshold value or threshold type", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      shinyjs::hide("correlation_matrix")
      shinyjs::hide("heatmap")
    }else{
      output$correlation_matrix <- renderPlot({
        correlation_input()
      })
      output$heatmap <- renderPlot({
        heatmap_input()
      })
      
      # if(!is.null(heatmap_input())){
      shinyalert("Complete Visualization!", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      shinyjs::show("correlation_matrix")
      shinyjs::show("heatmap")  
      # }
    }
    
    data_results()
  })
  
  observeEvent(input$show_sampleID ,{
    if(input$show_sampleID){
      output$pca_plot <- renderPlot({
        pca_input_Sample()
      })
    } else {
      output$pca_plot <- renderPlot({
        pca_input_noSample()
      })
    }
  })
  
  output$download_pca <- downloadHandler(
    filename = function() {paste0("PCA_plot_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(pca_input())){
        png(file)
        print(pca_input())
        dev.off()
      }
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() {paste0("Heatmap_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(heatmap_input())){
        png(file)
        dev.off()
      }
    }
  )
  
  output$download_gene_cluster <- downloadHandler(
    filename = function() {paste0("Heatmap_gene_cluster_info_",Sys.Date(),".csv")},
    content = function(file){
      filtered = subset(rowData(dep()), cut = significant == T)
      if(!is.null(filtered)){
        res <- get_gene_cluster(assay(dep()), rowData(dep()), input$dea_clusterNum)
        write.csv(res,file,row.names = F, quote = F)
      } 
    }
  )
  
  output$download_correlation <- downloadHandler(
    filename = function() {paste0("Correlation_plot_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(correlation_input())){
        png(file)
        print(correlation_input())
        dev.off()
      }
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {paste0("Volcano_plot_", Sys.Date(), ".png")},
    content = function(file) {s
      if(!is.null(volcano_input())){
        png(file)
        print(volcano_input())
        dev.off()
      }
    }
  )
  
  observeEvent(input$gsa_btn,{
    if(!is.null(dep())){
      shinyalert("Start GSA!","Please wait for a while", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      result_gsa()
      result_gsa_GOBP()
      gobp_plot()
      result_gsa_GOCC()
      gocc_plot()
      result_gsa_GOMF()
      gomf_plot()
      result_gsa_Kegg()
      kegg_plot()
      reverted_kegg()
    }else{
      shinyalert("Ther is no DEP!", "Change threshold value or threshold type", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
  })
  
  output$download_gobp <- downloadHandler(
    filename = function() {paste0("GO_BP_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_GOBP())){
        write.csv(result_gsa_GOBP(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gocc <- downloadHandler(
    filename = function() {paste0("GO_CC_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_GOCC())){
        write.csv(result_gsa_GOCC(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gomf <- downloadHandler(
    filename = function() {paste0("GO_MF_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_GOMF())){
        write.csv(result_gsa_GOMF(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_kegg <- downloadHandler(
    filename = function() {paste0("Kegg_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_Kegg())){
        write.csv(result_gsa_Kegg(),file,row.names=F,quote=F)
      }
    }
  )
  
  observeEvent(input$gsea_btn, {
    if(!is.null(dep())){
      condition <- make_condition(case_samples(),control_samples())
      group_name <- unique(condition)
      
    } else{
      shinyalert("Ther is no Data!", "Chack your data again", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
  })
  
  observeEvent(input$gsa_btn, {
    kegg_info <- reverted_kegg()
    pathway_choices <- kegg_info$Term
    updateSelectInput(session, "pathID_selector",
                      choices = pathway_choices, selected = "")
    
    output$topOfKeggDT <- DT::renderDataTable({
      DT::datatable(top_of_kegg(), options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)),
                    selection ="single") 
    }, server=T)
  })
  
  
  
  
  
  # observeEvent(input$pathID_selector,{
  #   if(!is.null(input$topOfKeggDT_rows_selected)){
  #     dtProxy = DT::dataTableProxy("topOfKeggDT", session=session)
  #     DT::reloadData(dtProxy, clearSelection=c("all"))
  #     
  #     if(input$pathID_selector != ""){
  #       output$pathview_result <- renderImage({
  #         outfile <- pathway_graph()
  #         
  #         list(src=outfile, contentType="image/png",
  #              width="100%", height="100%",
  #              alt="Pathview_graph")}, deleteFile=F)
  #       
  #     }
  #   }
  # })
  # 
  # 
  # observeEvent(input$topOfKeggDT_rows_selected, {
  #   if(!is.null(input$topOfKeggDT_rows_selected)){
  #     output$pathview_result <- renderImage({
  #       outfile <- pathway_graph()
  #       list(src=outfile, contentType="image/png",
  #            width="100%", height="100%",
  #            alt="Pathview_graph")}, deleteFile=F)
  #   }
  # })
  
  observeEvent(input$pathID_selector, {
    if(!is.null(input$topOfKeggDT_rows_selected)){
      dtProxy = DT::dataTableProxy("topOfKeggDT", session=session)
      DT::reloadData(dtProxy, clearSelection=c("all"))
    }
  })
  
  observeEvent(input$render_pathway_btn,{
    if(input$pathID_selector != "") {
      outfile <- pathway_graph()
      
      output$pathview_result <- renderImage({
        list(src=outfile, contentType="image/png",
             width="100%", height="100%",
             alt="Pathview_graph")}, deleteFile=F)
    } 
  })

  observeEvent(input$topOfKeggDT_rows_selected, {
    if(input$pathID_selector!=""){
      req(input$pathID_selector)
      updateSelectInput(session, "pathID_selector",
                        selected = "")
      dtProxy = DT::dataTableProxy("topOfKeggDT", session=session)
      DT::selectRows(dtProxy, selected=input$topOfKeggDT_rows_selected)
    }
    
    outfile <- pathway_graph()
          
    output$pathview_result <- renderImage({
      list(src=outfile, contentType="image/png",
           width="100%", height="100%",
           alt="Pathview_graph")}, deleteFile=F)
  })
  
  output$download_pathview <- downloadHandler(
    filename = function() {
      paste0("Pathview_", pathway_graph())},
    # paste0("Pathview_", selected_pathway(), ".png")},
    content = function(file) {
      file.copy(pathway_graph(), file)
    },
    contentType = "image/png"
  )
  
  
  observeEvent(input$dimension,{
    req(input$zoom_pathway_btn)
    width <- (input$dimension[1])*0.7
    height <- (input$dimension[2])*0.7
    
    showModal(modalDialog(
      renderImage({
        outfile <- pathway_graph()
        # width <- (session$clientData$output_download_pathview_width)*0.7
        # height <- (session$clientData$output_download_pathview_height)*0.7
        
        list(src=outfile, contentType="image/png",
             width=width, height=height,
             alt="Pathview_graph")}, deleteFile=F),
      easyClose = TRUE,
      footer=NULL
    ))
  })
  
  
  
  
  
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
  
  res_test <- eventReactive(input$test_btn,{
    req(preprocessed_data())
    data_diff <- test_diff(ready_for_dea(), type="all")
    if(input$test_method != "Limma"){
      diff_rowData <- rowData(data_diff)
      diff_colData <- colData(data_diff)
      
      data <- assay(data_diff)
      df <- test(data,diff_colData,input$test_method,input$padj_method)
      pval_pos <- grep("p.val",colnames(diff_rowData),fixed = T)
      padj_pos <- grep("_p.adj",colnames(diff_rowData),fixed = T)
      
      diff_rowData[,pval_pos] <- df$pval
      diff_rowData[,padj_pos] <- df$padj
      
      data_diff_edit <- SummarizedExperiment(assays = list(assay(data_diff)), rowData = diff_rowData, colData = diff_colData)
    } else{
      data_diff <- data_diff
    }
  })
  
  ### if we need z-score use these code!!
  # data_add_rejections <- reactive({
  #   type <- input$thres_type
  #   pvalue <- input$dea_pvalue
  #   log2fc <- input$dea_log2fc
  #   data_rejection <- add_rejections(res_test(), alpha=pvalue, lfc=log2fc)
  #   if(type == "pvalue"){
  #     dep_rowData <- rowData(data_rejection)
  #     pval_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
  #     log2fc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
  #     sig_pos <- grep("significant",colnames(dep_rowData),fixed = T)
  #     
  #     pval_df <- as.numeric(dep_rowData[,pval_pos])
  #     log2fc_df <- abs(as.numeric(dep_rowData[,log2fc_pos]))
  #     res_pval <- as.logical(pval_df < pvalue)
  #     res_lfc <-  as.logical(log2fc_df > log2fc)
  #     TT_pos <- which(res_pval=="TRUE" & res_lfc=="TRUE")
  #     
  #     res <- c()
  #     res[1:nrow(dep_rowData)] <- "FALSE"
  #     res[TT_pos] <- "TRUE"
  #     res <- as.logical(res)
  #     
  #     dep_rowData[,sig_pos] <- res
  #     write.csv(dep_rowData,"D:/dep_rowData.csv",quote = F)
  #     data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
  #   }
  #   data_rejection
  # })
  
  dep <- eventReactive(input$dea_btn,{
    req(preprocessed_data())
    
    ### if we need z-score use these code!!
    # first_dep <- data_add_rejections()
    
    # if(input$normalization == "No"){
    #   vsn_data <- rowData(first_dep)
    #   vsn_data <- vsn_data[,c(1,5:9)]
    #   znrom_res_test <- test_diff(ready_for_dea(), type="all")
    #   znorm_data <- rowData(znrom_res_test)
    #   znorm_data <- znorm_data[,c(1:4)]
    #   norms_rowData <- merge(znorm_data,vsn_data,by="name")
    #   first_dep <- SummarizedExperiment(assays = list(assay(znrom_res_test)), rowData = norms_rowData, colData = colData(znrom_res_test))
    # }
    # return(first_dep)
    
    type <- input$thres_type
    pvalue <- input$dea_pvalue
    log2fc <- input$dea_log2fc
    data_rejection <- c()
    
    if(type == "P.adj"){
      data_rejection <- add_rejections(res_test(), alpha=pvalue, lfc=log2fc)
      dep_rowData <- rowData(data_rejection)
      dep_rowData$name <- as.character(dep_rowData$name)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    }
    else if(type == "P.value"){
      data_rejection <- add_rejections(res_test(), alpha=pvalue, lfc=log2fc)
      dep_rowData <- change_Sig(data_rejection,pvalue,log2fc)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    }
    else {
      data_rejection <- add_rejections(res_test(), alpha=1, lfc=0)
      dep_rowData <- rowData(data_rejection)
      dep_rowData$name <- as.factor(dep_rowData$name)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    } 
    return(data_rejection)
  })

  
  volcano_input <- reactive({
    condition <- dep()$condition
    case_name <- condition[1]
    ctrl_name <- condition[length(condition)]
    contrast <- paste0(case_name,"_vs_", ctrl_name)

    dep_rowData <- rowData(dep())
    pv_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
    lfc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
    
    input_vc <- data_frame(name=dep_rowData$name, lfc=dep_rowData[,lfc_pos], 
                           p=-log10(as.numeric(dep_rowData[,pv_pos])), sig=dep_rowData$significant)
    plot_volcano(dep(), contrast=contrast, label_size=2, add_names=F) + 
      geom_point(data = filter(input_vc,sig), aes(lfc, p), color = "red", size= 2)+
      labs(x=expression(paste(log[2],FoldChange)),y=expression(paste(-log[10],P.value)))
  })
  
  pca_input_noSample <- reactive({
    n = nrow(assay(dep()))
    plot_pca(dep(), x=1, y=2, n=n, point_size=3, indicate = c("condition"))+
      ggtitle(paste("PCA plot -", n, "variable proteins", sep=" "))
  })
  
  pca_input_Sample <- reactive({
    n = nrow(assay(dep()))
    pca_df <- get_pca_df(assay(dep()), colData(dep()), n)
    plot_pca(dep(), x=1, y=2, n=n, point_size=3, indicate = c("condition"))+
      ggtitle(paste("PCA plot -", n, "variable proteins", sep=" "))+
      geom_text(data=pca_df,aes(label=rowname),nudge_y = 1.5)
  })
  
  correlation_input <- reactive({
    plot_cor(dep(), significant = T, pal="Blues", lower=-1, upper=1, indicate=c("condition"))
  })
  
  heatmap_input <- reactive({
    plot_heatmap(dep(), type="centered", kmeans=T, k=input$dea_clusterNum,
                 col_limit=2, show_row_names=F, indicate=c("condition"))
    
  })
  
  data_results <- reactive({
    data_results <- get_results(dep())
  })
  
  
  
  selected_pathway <- reactive({
    if(!is.null(input$topOfKeggDT_rows_selected)){
      kegg_info <- top_of_kegg()
      selected_pathway <- kegg_info[input$topOfKeggDT_rows_selected, 1]
    } else if(input$pathID_selector !="") {
      selected_pathway <- input$pathID_selector  
    } 
    return(selected_pathway)
  })
  
  top_of_kegg <- reactive({
    kegg_info <- reverted_kegg()
    top_kegg <- kegg_info[order(kegg_info$Adjusted.P.value, decreasing=F),]
    top_kegg <- head(top_kegg, 10)
    top_kegg <- dplyr::select(top_kegg, -P.value, -matches("Old|Odds|Score"))
    return(top_kegg)
  })
  
  
  pathway_graph <- eventReactive(list(input$render_pathway_btn, input$topOfKeggDT_rows_selected),{
    pathway_name <- selected_pathway()
    kegg_info <- reverted_kegg()
    rowdt <- rowData(dep())

    fc_cols <- colnames(rowdt)[grepl("diff",colnames(rowdt))]
    fc <- rowdt[,fc_cols]
    names(fc) <- rowdt$name

    pathid <- as.character(kegg_info[kegg_info$Term==pathway_name, "kegg_id"])

    outfile <- paste0("hsa", pathid,".pathview.png")

    pathview(fc, pathway.id=pathid, gene.idtype="SYMBOL", species = "hsa",
             kegg.dir="./PATHVIEW/")

    return(outfile)
  })
  
  
  output$volcano_info <- DT::renderDataTable(DT::datatable({
    dep_rowData <- rowData(dep())
    pv_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
    lfc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
    
    input_vc <- data.frame(name=dep_rowData$name, log2FC=dep_rowData[,lfc_pos], 
                           `-log10P.val`=-log10(as.numeric(dep_rowData[,pv_pos])), sig=dep_rowData$significant)
    
    input_vc <- data.frame(input_vc)
    colnames(input_vc) <- c("name", "log<sub>2</sub>FC", "-log<sub>10</sub>P.value", "sig")
    # input_vc <- input_vc %>% filter(input_vc, sig)
    df <- brushedPoints(input_vc, input$volcano_brush, xvar="log<sub>2</sub>FC", yvar = "-log<sub>10</sub>P.value")
    # return(df)
  }, options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), escape = F, selection="none") %>%  DT::formatRound(c(2:3),digits=2))
  
  
  
  ########################### GSA ######################################
  result_gsa <- reactive({
    data <- rowData(dep())
    res_gsa <- gsa(data,input$gsa_set,input$gsa_tool)
    return(res_gsa)
  })
  
  result_gsa_GOBP <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      if(input$gsa_tool == "enrichR"){
        gobp<-res_gsa[["GO_Biological_Process_2018"]]
      }
    }
    return(gobp)
  })
  
  gobp_plot <- reactive({
    gobp <- result_gsa_GOBP()
    gobp$P.value2 <- -log10(gobp$P.value)
    gobp <- gobp[order(-gobp$P.value2),]
    gobp <- gobp[c(1:input$set_nterm),]
    output$gobp_plot <- renderPlot({
      ggplot(data=gobp, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
        geom_bar(stat = "identity",fill="#3c8dbc")+
        labs(title="GO_BP",x=expression(paste(-log[10],P.value)),y="")+
        theme_bw()+
        theme(axis.text=element_text(size=12),title = element_text(size=15,face="bold"))
    })
  })
  
  result_gsa_GOCC <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      if(input$gsa_tool == "enrichR"){
        gocc<-res_gsa[["GO_Cellular_Component_2018"]]
      }
    }
    return(gocc)
  })
  
  gocc_plot <- reactive({
    gocc <- result_gsa_GOCC()
    gocc$P.value2 <- -log10(gocc$P.value)
    gocc <- gocc[order(-gocc$P.value2),]
    gocc <- gocc[c(1:input$set_nterm),]
    output$gocc_plot <- renderPlot({
      ggplot(data=gocc, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
        geom_bar(stat = "identity",fill="#3c8dbc")+
        labs(title="GO_CC",x=expression(paste(-log[10],P.value)),y="")+
        theme_bw()+
        theme(axis.text=element_text(size=12),title = element_text(size=15,face="bold"))
    })
  })
  
  result_gsa_GOMF <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      if(input$gsa_tool == "enrichR"){
        gomf<-res_gsa[["GO_Molecular_Function_2018"]]
      }
    }
    return(gomf)
  })
  
  gomf_plot <- reactive({
    gomf <- result_gsa_GOMF()
    gomf$P.value2 <- -log10(gomf$P.value)
    gomf <- gomf[order(-gomf$P.value2),]
    gomf <- gomf[c(1:input$set_nterm),]
    output$gomf_plot <- renderPlot({
      ggplot(data=gomf, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
        geom_bar(stat = "identity",fill="#3c8dbc")+
        labs(title="GO_MF",x=expression(paste(-log[10],P.value)),y="")+
        theme_bw()+
        theme(axis.text=element_text(size=12),title = element_text(size=15,face="bold"))
    })
  })
  
  result_gsa_Kegg <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      if(input$gsa_tool == "enrichR"){
        kegg<-res_gsa[["KEGG_2019_Human"]]
      }
    }
    return(kegg)
  })
  
  reverted_kegg <- reactive({
    kegg <- result_gsa_Kegg()
    kegg <- changePathwayID(kegg)
    return(kegg)
  })
  
  
  kegg_plot <- reactive({
    kegg <- result_gsa_Kegg()
    kegg$P.value2 <- -log10(kegg$P.value)
    kegg <- kegg[order(-kegg$P.value2),]
    kegg <- kegg[c(1:input$set_nterm),]
    output$kegg_plot <- renderPlot({
      ggplot(data=kegg, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
        geom_bar(stat = "identity",fill="#3c8dbc")+
        labs(title="KEGG",x=expression(paste(-log[10],P.value)),y="")+
        theme_bw()+
        theme(axis.text=element_text(size=12),title = element_text(size=15,face="bold"))
    })
  })
  
  addTimeLine <-  function(timeLine){
    output$timeline <- renderUI({
      timelineBlock(
        reversed = F,
        timelineEnd(icon="hourglass-start", color = "gray"),
        lapply(1:nrow(timeLine), FUN = function(i){
          tagList(
            timelineItem(
              icon = timeLine$icon[i],
              color = timeLine$color[i],
              time = timeLine$time[i],
              tags$h4(timeLine$step[i]),
              footer= timeLine$info[i]
              
            )
            ,timelineLabel(paste0("# of proteins : ",timeLine$sample_num[i]), color = "purple")
          )
        }),
        timelineStart(icon="hourglass-end", color = "gray")
      )
    })
  } # End of addTimeLine
})