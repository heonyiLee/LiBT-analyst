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

get_main_data <- function(data){
  db <- read_csv("./base/uniprotKB_idmapping_human9606_201912ver.csv")
  
  ProteinID <- as.character(data$Majority.protein.IDs)
  Protein_count <- as.character(data$Peptide.counts..all.)
  GeneName <- as.character(data$Gene.names)
  
  multi_pos <- grep(";",GeneName,fixed = T)
  ProteinID_multi <- ProteinID[multi_pos] 
  Protein_count_multi <- Protein_count[multi_pos]
  GeneName_multi <- GeneName[multi_pos]
  
  pt_multi <- c()
  for(i in 1:length(ProteinID_multi)){
    pt_id <- ProteinID_multi[i]
    pt_id <- data.frame(do.call('rbind', strsplit(as.character(pt_id), split = ';', fixed = TRUE)))
    pt_id <- as.character(t(pt_id))
    
    pt_count <- Protein_count_multi[i]
    pt_count <- data.frame(do.call('rbind', strsplit(as.character(pt_count), split = ';', fixed = TRUE)))
    pt_count <- as.numeric(t(pt_count))
    pt_count <- pt_count[1:length(pt_id)]
    
    pt <- data.frame(accid=pt_id,count=pt_count)
    id_pos <- as.numeric(grep("-",pt$accid,fixed = T))
    if(length(id_pos)!=0){
      edit_id <- as.character(pt$accid[id_pos])
      edit_id <- data.frame(do.call('rbind', strsplit(edit_id, split = '-', fixed = TRUE)))
      pt_id[id_pos] <- as.character(edit_id[,1])
      pt$accid <- pt_id
    }
    
    pt_db <- merge(pt,db,by="accid")
    pt_db <- pt_db[pt_db$reviewed=="Y",]
    pt_db <- pt_db[order(pt_db$count,decreasing = T),]
    pt_db <- pt_db[pt_db$count==pt_db$count[1],]
    pt_db <- pt_db[!duplicated(pt_db$gene),]
    
    id <- c()
    gn <- c()
    if(nrow(pt_db)>1){
      for(j in 1:nrow(pt_db)){
        id <- paste0(id,";",pt_db$accid[j])
        gn <- paste0(gn,";",pt_db$gene[j])
      }
      pt_id <- gsub("^;","",id)
      pt_gn <- gsub("^;","",gn)
      gn_multi <- data.frame(ProteinID=pt_id,GeneName=pt_gn)
    }  else if(nrow(pt_db) == 0){
      pt_gn <- NA
      gn_multi <- data.frame(ProteinID=ProteinID_multi[i],GeneName=GeneName_multi[i])
    }  else{
      pt_gn <- pt_db$gene
      gn_multi <- data.frame(ProteinID=pt_db$accid,GeneName=pt_gn)
    }
    
    
    pt_multi <- rbind(pt_multi,gn_multi)
  }
  
  ProteinID[multi_pos] <- as.character(pt_multi$ProteinID)
  GeneName[multi_pos] <- as.character(pt_multi$GeneName)
  
  cn_data <- colnames(data)
  main_pos <- grep("LFQ",cn_data)
  main <- data[,main_pos]
  
  unique_pos <- grep("Unique.peptides",cn_data)
  Unique_peptides <- data[,unique_pos]
  Sequence_coverage <- data$Sequence.coverage....
  
  main_df <- cbind(ProteinID,GeneName,main,Unique_peptides,Sequence_coverage)
  return(main_df)
}


make_case_samples <- function(data,file_type) {
  # col <- data.frame(id=c())
  
  temp_col <- colnames(data)
  #temp_col <- tolower(temp_col)
  temp_col <- temp_col[grep(file_type, temp_col)]
  file_type <- as.character(paste0(file_type,".intensity"))
  #temp <- strsplit(temp_col, split=file_type)
  temp <- data.frame(do.call('rbind', strsplit(as.character(temp_col), split = file_type, fixed = TRUE)))
  
  # temp <- unlist(temp)
  # 
  # for(i in 1:length(temp)) {
  #   if(i %% 2 == 0) {
  #     t <- unlist(strsplit(temp[i], split="\\."))
  #     t <- data.frame(id=t[1])
  #     col <- rbind(col, t)
  #   }
  # }
  col <- data.frame(id=temp[,2])
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  col <- paste0(file_type, col)
  return(col)
  
  # col <- foreach(i=1:length(temp), .combine=rbind) %do%{
  #   if(i %% 2 == 0) {
  #     t <- temp [i]
  #   }
  #   return(col)
  #   
  # }
}