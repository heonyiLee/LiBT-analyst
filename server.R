library(dplyr)
library(readr)
library(stringr)
library(DEP)

source("function.R")

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
          print(input$TMT_input_option)
        }
      }
      
    }
    
  })
  
  
  observeEvent(input$file_upload_btn, {
    if(input$file_upload_btn==0) {
      shinyjs::hide("summary_box")
    }
    shinyjs::show("summary_box")
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
        print(samples)
        updateSelectInput(session, "case_group_selection", choices = samples)
      }
      # } else {
      #   shinyalert("Choose proper file type!", type="error", timer = 10000,
      #              closeOnClickOutside = T, closeOnEsc = T)
      # }
      # updateSelectInput(session, "control_group_selection", choices= samples,
      #                   label="Control samples")
    }
  })
  
  
  observeEvent(input$case_group_selection, {
    choices <- control_samples()
    updateSelectInput(session, "control_group_selection",
                      # label=paste0("Control samples (n=", length(choices), ")"),
                      choices = choices, selected = choices)
  })
  
  observeEvent(input$exp_design_submit_btn, {
    print(case_samples())
    print("---------------")
    print(control_samples())
    
    
    updateSelectInput(session, "control_group_selection",
                      label=paste0("Control samples (n=", length(control_samples()), ")"),
                      choices = control_samples(), selected = control_samples())
    
    updateSelectInput(session, "case_group_selection",
                      label=paste0("Case samples (n=", length(case_samples()), ")"),
                      choices = case_samples(), selected = case_samples())
    
    
    tmp <- paste0("# of Case : ", length(case_samples()),"\n",
                  "# of Control : ", length(control_samples()))
    info <- paste0(info,tmp,"\n")
    newTL <- data.frame(step="Exp Design Submit",
                        info=info,
                        sample_num=as.numeric(nrow(main_data())),
                        time=as.character(Sys.time()),color="maroon")
    timeLine <- rbind(timeLine,newTL)
    addTimeLine(timeLine)
    
    # output$sidebar_tab2 <- renderUI({
    #   
    # })
    # shinyjs::show("rightsidebar_2")
    
  })
  
  observeEvent(input$preprocess_btn, {
    temp_df <- preprocessed_data()
    output$preprocessed_data_header <- DT::renderDataTable({
      temp_df}, 
      options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
    
    tmp <- paste0("Transformation : ", str_to_title(input$transformation),"\n",
                  "Valid value : ", paste0(input$valid_value,"%"), "\n",
                  "Imputation : " , str_to_title(input$imputation), "\n",
                  "Normalization : ", str_to_title(input$normalization))
    info <- paste0(info,tmp,"\n")
    newTL <- data.frame(step="Preprocessing",
                        info=info,
                        sample_num=as.numeric(nrow(temp_df)),
                        time=as.character(Sys.time()),color="maroon")
    timeLine <- rbind(timeLine,newTL)
    addTimeLine(timeLine)
  })
  

 
### LOAD DATA
  file_input <- reactive({NULL})
  file_input <- eventReactive(input$fileBrowser, {
    req(input$fileBrowser)
    
    if(is.null(input$fileBrowser)){
      return(NULL)
    }
    temp_df <- readLines(input$fileBrowser$datapath, n=1)
    if(grepl(" ", temp_df)){
      sep <- c(" ")
    } else if(grepl("\t", temp_df)) {
      sep <- c("\t")
    } else if(grepl(";", temp_df)) {
      sep <- c(";")
    } else if(grepl(",", temp_df)) {
      sep <- c(",")
    } 
    
    temp_df <- read.table(input$fileBrowser$datapath,
                          header = T, fill = T,
                          sep = sep)
    state <- file_input_test(temp_df, input$file_type)
    print(ncol(temp_df))
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
      print(input$TMT_input_option)
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
      if(input$TMT_input_option=="F" || is.null(input$TMT_input_option)){
        samples <- make_case_samples_diffT(df, input$TMT_input_option)
      } else {
        samples <- make_case_samples_T(df)
      }
    } else{
      samples <- make_case_samples_LiB(df)
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
  

  preprocessed_data <- eventReactive(input$preprocess_btn, {
    transformation_item <- input$transformation
    filter_item <- input$valid_value
    imputation_item <- input$imputation
    normalization_item <- input$normalization
    
    case_samples <- case_samples()
    control_samples <- control_samples()
    temp_df <- main_data()
    
    preprocessing_option <- list(c(transformation_item, filter_item,
                                   imputation_item, normalization_item))
    
    temp_df <- get_preprocessed_data(temp_df, preprocessing_option, case_samples, control_samples)
    return(temp_df)
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
              footer= timeLine$info[i],
              
            )
            ,timelineLabel(paste0("# of samples : ",timeLine$sample_num[i]), color = "olive")
          )
        }),
        timelineStart(color = "gray")
      )
    })
  } # End of addTimeLine
})