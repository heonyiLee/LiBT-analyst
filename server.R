library(shiny)
# Output: Histogram ----

shinyServer(function(input,output, session){

  step <- c()
  dType <- c()
  sample_num <- c()
  time <- c()
  
  
  # output$data_type <- renderText({
  #     input$data_type
  # })
  
  output$input_data <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$input_csv)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        df <- read.csv(input$input_csv$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
        
        timeLine <- data.frame(step="Data Input",dType=input$data_type,sample_num=as.numeric(nrow(df)),time=as.character(Sys.time()))
        print(timeLine)
        #step <- c(step,"Data Input")
        addTimeLine(timeLine)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
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
              title = timeLine$step,
              icon = "gears",
              color = "olive",
              time = timeLine$time,
              footer = paste0("Data Type : \n ", timeLine$dType)
            )
          )
        }),
        timelineStart(color = "gray")
      )
    })
  }
})