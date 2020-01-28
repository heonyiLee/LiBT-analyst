filter_with_condition <- function(selected, temp_data) {
  
  if(length(selected)==1 & selected[1] == "potential") { #4844
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+")
    
  } else if(length(selected)==1 & selected[1] == "reverse") { #4844
    
    temp_data <- dplyr::filter(temp_data, Reverse != "+")
    
  } else if(length(selected)==1 & selected[1] == "identified") { #4877
    
    temp_data <- dplyr::filter(temp_data, Only.identified.by.site != "+")
    
  } else if(length(selected)==2 & 
            (selected[1] == "potential" & selected[2] == "reverse")) { #4793
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                 Reverse != "+")
  } else if(length(selected)==2 & 
            (selected[1] == "potential" & selected[2] == "identified")) { # 4826
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                 Only.identified.by.site != "+")
  } else if(length(selected)==2 & 
            (selected[1] == "reverse" & selected[2] == "identified")) { #4829
    
    temp_data <- dplyr::filter(temp_data, Reverse != "+" &
                                 Only.identified.by.site != "+")
  } else if(length(selected)==3 & 
            (selected[1] == "potential" &
             selected[2] == "reverse" & selected[3] == "identified")) { #4805
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                 Only.identified.by.site != "+" & Reverse != "+")
    
  } else if(length(selected)==0) {
    temp_data <- temp_data
  }
  return(temp_data)
}




make_case_samples <- function(data) {
  col <- data.frame(id=c())
  
  temp_col <- colnames(data)
  temp_col <- tolower(temp_col)
  temp_col <- temp_col[grep("total_", temp_col)]
  temp <- strsplit(temp_col, split="total")
  temp <- unlist(temp)
  
  for(i in 1:length(temp)) {
    if(i %% 2 == 0) {
      t <- unlist(strsplit(temp[i], split="\\."))
      t <- data.frame(id=t[1])
      col <- rbind(col, t)
    }
  }
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  col <- paste0("Total", col)
  return(col)

}


