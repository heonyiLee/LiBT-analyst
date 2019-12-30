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
                                  sep = "\t")
      
      if(input$first_filtering == 1) {
        temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+")
      } else if(input$first_filtering == 2) {
        temp_data <- dplyr::filter(temp_data, Reverse != "+")
      } else if(input$first_filtering == 3) {
        temp_data <- dplyr::filter(temp_data, Only.identified.by.site != "+")
      } else if(input$first_filtering == 4) {
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
  
  observeEvent(input$first_filtering==4, {
    print("ALL")
    
  })
  
  

  


  
  
  
  
  addTimeLine = function(timeLine){
    output$timeline <- renderUI({
      timelineBlock(
        reversed = F,
        timelineEnd(color = "danger"),
        apply(timeLine, 1, FUN = function(i){
          tagAppendAttributes(
            #timelineLabel(timeLine$sample_num, color = "teal"),
            timelineItem(
              style="width:400px;",
              title = paste0(timeLine$step),
              icon = "gears",
              color = "olive",
              # time = timeLine$time,
              footer = paste0("File Type : \n ", timeLine$fType)
            )
          )
        }),
        timelineStart(color = "gray")
      )
    })
  } # End of addTimeLine
  
  
  
  
})