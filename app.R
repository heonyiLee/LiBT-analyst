ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE,logical.return=TRUE)
  require_pkg <- sapply(pkg, require, character.only = TRUE)
  isFalse <- as.numeric(which(require_pkg == F))
  if(length(isFalse) != 0){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg[isFalse])
  }
}

#DEP, edgeR, gage, fgsea, pathview install by BiocManager
pkg <- c("shiny","shinyjs","shinydashboard","shinydashboardPlus",
         "shinyWidgets","shinyalert","shinycssloaders","shinyjqui",
         "sortable","dplyr","readr","stringr",
         "DEP","sortable","SummarizedExperiment",
         "ggplot2","edgeR","factoextra","enrichR","tibble",
         "gage","fgsea","pathview","zip")
ipak(pkg)

shinyApp(ui=ui, server=server)
