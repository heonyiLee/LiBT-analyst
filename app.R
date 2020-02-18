ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE,logical.return=TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkg <- c("shiny", "shinyjs", "shinydashboard", "shinydashboardPlus", "shinyalert", 
         "shinycssloaders", "stringr",
         "readr","dplyr", "RColorBrewer", "d3heatmap", "viridis")
ipak(pkg)


shinyApp(ui=ui, server=server)
