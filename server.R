library(shiny)
library(shinyjs)
# Output: Histogram ----

shinyServer(function(input,output, session){
  
  options(shiny.maxRequestSize=100*1024^2)

  step <- c()
  fType <- c()
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
  

  uploaded_data <- reactive({NULL})
  
  uploaded_data <- eventReactive(input$file_upload_btn, {
    req(input$fileBrowser)
    
    if(is.null(input$fileBrowser)){
      return(NULL)
    }
    tryCatch({
      temp_data <- read.table(input$fileBrowser$datapath,
                                  header = T,
                                  fill = T,
                                  sep = "\t") #4901

      temp_data <- dplyr::filter(temp_data, Unique.peptides != 0)
      temp_data <- dplyr::filter(temp_data, Intensity != 0)
      
      if(input$first_filtering == "potential") { #4868
        temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+")
      } else if(input$first_filtering == "reverse") { #4868
        temp_data <- dplyr::filter(temp_data, Reverse != "+")
      } else if(input$first_filtering == "identified") { #4901
        temp_data <- dplyr::filter(temp_data, Only.identified.by.site != "+")
      } else if(input$first_filtering == "select_all") {
        temp_data <- temp_data
      } 
      
      timeLine <- data.frame(step="Data Input",fType=input$file_type,
                             sample_num=as.numeric(nrow(temp_data)),
                             time=as.character(Sys.time()))
      addTimeLine(timeLine)
      
    }, error=function(e){
      stop(safeError(e))
    })
    return(temp_data)
  })
  
  
  
  observeEvent(input$file_upload_btn, {
    if(is.null(input$first_filtering)){
      shinyalert("Choose filetering option!", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    } else {
      output$uploaded_file_header <- DT::renderDataTable({
        uploaded_data()
      }, options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
      
      print(nrow(uploaded_data()))
    }
    
  })
  
  
  observeEvent(input$preprocess_btn, {
    
    
  })
  
  # observeEvent(input$first_filtering, {
  #   if(input$first_filtering=="select_all"){
  #     print("ALL")
  #     updateCheckboxGroupInput(session, "first_filtering",
  #                              selected = c("potential","reverse","identified","select_all"))
  #   } else if(input$first_filtering=="remove_all"){
  #     updateCheckboxGroupInput(session, "first_filtering",
  #                              selected = c())
  #   } else {
  #     NULL
  #   }
  # 
  # })

  
  
  addTimeLine = function(timeLine){
    output$timeline <- renderUI({
      timelineBlock(
        reversed = F,
        timelineEnd(color = "danger"),
        apply(timeLine, 1, FUN = function(i){
          tagList(
            timelineItem(
              icon = "file-upload",
              color = "olive",
              time = timeLine$time,
              tags$h4(timeLine$step),
              footer= tagList(
                tags$p(paste0("File Type : ", timeLine$fType)),
              )
            )
            ,timelineLabel(timeLine$sample_num, color = "teal")
          )
        }),
        timelineStart(color = "gray")
      )
    })
  } # End of addTimeLine
  
  
  
  
})