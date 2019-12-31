library(shiny)
library(shinyjs)
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
      updateCheckboxGroupInput(session, "first_filtering",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c("potential","reverse","identified"))
    }
  })
  
  observeEvent(input$deselect_all_filtering_btn, {
    print("deselect")
    if(input$deselect_all_filtering_btn == 0) {
      return(NULL)
    } else if(input$deselect_all_filtering_btn >0){
      updateCheckboxGroupInput(session, "first_filtering",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c())
    }
  })
 
  
  observeEvent(input$file_upload_btn, {
    if(is.null(input$first_filtering) | is.null(input$fileBrowser)){
      shinyalert("Choose option!", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    } else {
      render_df <- filtered_data()
      output$uploaded_file_header <- DT::renderDataTable({
        render_df
      }, options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))

      print(nrow(render_df))
      filter <- paste0(as.character(unlist(input$first_filtering)))
      for(i in 1:length(filter)){
        if(filter[i]=="potential"){
          tmp <- "'Potential contaminant'"
          info <- paste0(info,tmp," remove\n")
        }
        else if(filter[i]=="reverse"){
          tmp <- "'Reverse'"
          info <- paste0(info,tmp," remove\n")
        }else{
          tmp <- "'Only identified by site'"
          info <- paste0(info,tmp," remove\n")
        }
        
      }
      newTL <- data.frame(step="First Filter",
                          info=info,
                          sample_num=as.numeric(nrow(render_df)),
                          time=as.character(Sys.time()),color="maroon")
      timeLine <- rbind(timeLine,newTL)
      addTimeLine(timeLine)
    }
  })
  

  uploaded_data <- reactive({NULL})
  uploaded_data <- eventReactive(input$fileBrowser, {
    req(input$fileBrowser)

    if(is.null(input$fileBrowser)){
      return(NULL)
    }
    tryCatch({
      uploaded_data <- read.table(input$fileBrowser$datapath,
                                  header = T, fill = T,
                                  sep = "\t") 
      
      uploaded_data <- dplyr::filter(uploaded_data, Unique.peptides != 0)#4903
      uploaded_data <- dplyr::filter(uploaded_data, Intensity != 0)#4898

      info <- paste0("File Type : ", input$file_type,"\n","'Unique peptied' == 0 remove\n'Intensity' == 0 remove")
      timeLine <<- data.frame(step="Data Input",info=info,
                             sample_num=as.numeric(nrow(uploaded_data)),
                             time=as.character(Sys.time()),color="maroon")
      addTimeLine(timeLine)

    }, error=function(e){
      stop(safeError(e))
    })
    return(uploaded_data)
  })
  
  
  filtered_data <- reactive({NULL})
  filtered_data <- eventReactive(input$first_filtering, {
    
    temp_data <- uploaded_data()
    checked <- input$first_filtering
    
    if(length(checked)==1 & checked[1] == "potential") { #4868
      temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+")
    } else if(length(checked)==1 & checked[1] == "reverse") { #4868
      temp_data <- dplyr::filter(temp_data, Reverse != "+")
    } else if(length(checked)==1 & checked[1] == "identified") { #4901
      temp_data <- dplyr::filter(temp_data, Only.identified.by.site != "+")
    } else if(length(checked)==2 & 
              (checked[1] == "potential" & checked[2] == "reverse")) { #4817
      temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                   Reverse != "+")
    } else if(length(checked)==2 & 
              (checked[1] == "potential" & checked[2] == "identified")) { # 4850
      temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                   Only.identified.by.site != "+")
    } else if(length(checked)==2 & 
              (checked[1] == "reverse" & checked[2] == "identified")) { #4853
      temp_data <- dplyr::filter(temp_data, Reverse != "+" &
                                   Only.identified.by.site != "+")
    } else if(length(checked)==3 & 
              (checked[1] == "potential" &
               checked[2] == "reverse" & checked[3] == "identified")) { #4805
      temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                   Only.identified.by.site != "+" & Reverse != "+")
    } else if(length(checked)==0) {
      temp_data <- temp_data
    }
    
    return(temp_data)
  })
  
  
  observeEvent(input$preprocess_btn, {
    
    
  })
  

  addTimeLine = function(timeLine){
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