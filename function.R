filter_with_condition <- function(checked, temp_data) {
  
  if(length(checked)==1 & checked[1] == "potential") { #4844
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+")
    
  } else if(length(checked)==1 & checked[1] == "reverse") { #4844
    
    temp_data <- dplyr::filter(temp_data, Reverse != "+")
    
  } else if(length(checked)==1 & checked[1] == "identified") { #4877
    
    temp_data <- dplyr::filter(temp_data, Only.identified.by.site != "+")
    
  } else if(length(checked)==2 & 
            (checked[1] == "potential" & checked[2] == "reverse")) { #4793
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                 Reverse != "+")
  } else if(length(checked)==2 & 
            (checked[1] == "potential" & checked[2] == "identified")) { # 4826
    
    temp_data <- dplyr::filter(temp_data, Potential.contaminant != "+" &
                                 Only.identified.by.site != "+")
  } else if(length(checked)==2 & 
            (checked[1] == "reverse" & checked[2] == "identified")) { #4829
    
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
  
  # col <- foreach(i=1:length(temp), .combine=rbind) %do%{
  #   if(i %% 2 == 0) {
  #     t <- temp [i]
  #   }
  #   return(col)
  #   
  # }
}