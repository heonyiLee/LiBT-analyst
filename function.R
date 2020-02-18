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
  } else if(checked[1]=="TMT"){
    temp_data <- temp_data
  }
  return(temp_data)
}

get_main_data_LiB <- function(data) {
  rv_db <- read.delim("D:/proteomics/web/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("D:/proteomics/web/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")

  ProteinID <- as.character(data$Majority.protein.IDs)
  Protein_count <- as.character(data$Peptide.counts..all.)
  GeneName <- as.character(data$Gene.names)
  
  multi_pos <- grep(";",Protein_count,fixed = T)
  ProteinID_multi <- ProteinID[multi_pos] 
  Protein_count_multi <- Protein_count[multi_pos]
  
  pt_multi <- c()
  for(i in 1:length(ProteinID_multi)){
    pt_id <- ProteinID_multi[i]
    pt_id <- data.frame(do.call('rbind', strsplit(as.character(pt_id), split = ';', fixed = TRUE)))
    pt_id <- as.character(t(pt_id))
    
    pt_count <- Protein_count_multi[i]
    pt_count <- data.frame(do.call('rbind', strsplit(as.character(pt_count), split = ';', fixed = TRUE)))
    pt_count <- as.numeric(t(pt_count))
    pt_count <- pt_count[1:length(pt_id)]
    
    pt <- data.frame(Entry=pt_id,count=pt_count)
    id_pos <- as.numeric(grep("-",pt$Entry,fixed = T))
    if(length(id_pos)!=0){
      edit_id <- as.character(pt$Entry[id_pos])
      edit_id <- data.frame(do.call('rbind', strsplit(edit_id, split = '-', fixed = TRUE)))
      pt_id[id_pos] <- as.character(edit_id[,1])
      pt$Entry <- pt_id
    }
    
    pt_db <- merge(pt,rv_db,by="Entry")
    if(nrow(pt_db)!=0){
      pt_db <- pt_db[order(pt_db$count,decreasing = T),]
      pt_db <- pt_db[pt_db$count==pt_db$count[1],]
      pt_db <- pt_db[!duplicated(pt_db$Entry),]
    }else{
      pt_nrv <- merge(pt,nrv_db,by="Entry")
      if(nrow(pt_nrv) != 0){
        pt_nrv <- pt_nrv[pt_nrv$count==pt_nrv$count[1],]
        pt_db <- pt_nrv
        pt_db <- pt_db[!duplicated(pt_db$Entry),]
      }else{
        data_pos <- grep(ProteinID_multi[i],data$Majority.protein.IDs)
        pt_db <- data.frame(Entry=data$Majority.protein.IDs[data_pos],Gene.symbol=data$Gene.names[data_pos])
      }
    }
    
    id <- c()
    gn <- c()
    if(nrow(pt_db)>1){
      for(j in 1:nrow(pt_db)){
        id <- paste0(id,";",pt_db$Entry[j])
        gn <- paste0(gn,";",pt_db$Gene.symbol[j])
      }
      pt_id <- gsub("^;","",id)
      pt_gn <- gsub("^;","",gn)
      gn_multi <- data.frame(ProteinID=pt_id,GeneName=pt_gn)
    } else{
      pt_gn <- pt_db$Gene.symbol
      gn_multi <- data.frame(ProteinID=pt_db$Entry,GeneName=pt_gn)
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

get_main_data_T <- function(data,normalization) {
  rv_db <- read.delim("./base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("./base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")

  cn <- colnames(data)
  data$Accession <- as.character(data$Accession)
  data_class <- sapply(data,class)
  rm_pos <- grep("factor",data_class) 
  rm_pos <- c(rm_pos,grep("Ratio",cn))
  
  data <- data[,-rm_pos]
  cn <- cn[-rm_pos]
  main_pos <- grep("Accession",cn)

  TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  for(i in 1:length(TMT_marker)){
    main_pos <- c(main_pos,grep(TMT_marker[i],cn))  
  }
  main_pos <- unique(main_pos)
  
  main_data <- data[,main_pos]
  cn <- colnames(main_data)

  if(normalization == "T"){
    nor_pos <- grep("Normalized",cn)
    main_pos <- c(1,nor_pos)
    main_data <- main_data[,main_pos]
  } else{
    nor_pos <- grep("Normalized",cn)
    main_data <- main_data[,-nor_pos]
  }
  
  main_data <- na.omit(main_data)

  acc <- main_data$Accession
  id_pos <- grep("-",acc)
  acc_id <- acc[id_pos]
  
  acc_new <- c()
  for(j in 1:length(acc_id)){
    tmp <- acc_id[j]
    tmp <- unlist(strsplit(as.character(tmp), split = '-', fixed = TRUE))
    acc_new <- c(acc_new,tmp[1])
  }
  acc[id_pos] <- acc_new
  main_data$Entry <- acc

  pt_data <- merge(main_data,rv_db,by="Entry",all.x = T)
  pt_rv <- pt_data[!is.na(pt_data$Entry.name),]
  pt_rv <- pt_rv[,c(2,18,3:12)]
  
  pt_nrv <- pt_data[is.na(pt_data$Entry.name),]
  pt_nrv <- pt_nrv[,c(1:12)]
  pt_nrv <- merge(pt_nrv,nrv_db,by="Entry",all.x = T)
  
  pt_none <- pt_nrv[is.na(pt_nrv$Entry.name),]
  pt_none <- pt_none[,c(2,18,3:12)]
  
  pt_nrv <- pt_nrv[!is.na(pt_nrv$Entry.name),]
  pt_nrv <- pt_nrv[,c(2,18,3:12)]

  main_df <- rbind(pt_rv,pt_nrv,pt_none)
  return(main_df)
}


make_case_samples_LiB <- function(data,file_type) {
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

make_case_samples_T <- function(data) {
  data <- data[,c(3:ncol(data))]
  temp_col <- colnames(data)
  
  temp <- c()
  TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  for(i in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[i],temp_col)],TMT_marker[i]) 
    temp <- c(temp,substr(temp_col[i],col_pos[1,1],nchar(temp_col)))
  }
  
  col <- data.frame(id=temp)
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  
  return(col)
}

make_case_samples_diffT <- function(data,normalization) {
  cn <- colnames(data)
  data$Accession <- as.character(data$Accession)
  data_class <- sapply(data,class)
  rm_pos <- grep("factor",data_class) 
  rm_pos <- c(rm_pos,grep("Ratio",cn))
  
  data <- data[,-rm_pos]
  cn <- cn[-rm_pos]
  
  main_pos <- grep("Accession",cn)
  
  TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  for(i in 1:length(TMT_marker)){
    main_pos <- c(main_pos,grep(TMT_marker[i],cn))  
  }
  main_pos <- unique(main_pos)
  
  main_data <- data[,main_pos]
  cn <- colnames(main_data)
  if(normalization == "T"){
    nor_pos <- grep("Normalized",cn)
    main_pos <- c(1,nor_pos)
    main_data <- main_data[,main_pos]
  } else{
    nor_pos <- grep("Normalized",cn)
    main_data <- main_data[,-nor_pos]
  }
  
  main_data <- na.omit(main_data)
  main_data <- main_data[,c(2:ncol(main_data))]
  temp_col <- colnames(main_data)

  temp <- c()
  for(j in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[j],temp_col)],TMT_marker[j]) 
    temp <- c(temp,substr(temp_col[j],col_pos[1,1],nchar(temp_col)))
  }
  
  col <- data.frame(id=temp)

  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  print(col)
  return(col)
}