file_input_test <- function(file, file_type) {
  if(file_type=="TMT"){
    cols <- colnames(file)[grep("Abundance", colnames(file))]
  } else {
    cols <- colnames(file)[grep(file_type, colnames(file))]
  }
  
  return(cols)
}



filter_with_option <- function(option, data) {
  data <- dplyr::filter(data, Peptides != 0)#4903
  data <- dplyr::filter(data, Intensity != 0)#4898
  
  if(length(option)==1 & option[1] == "potential") { #4844
    
    data <- dplyr::filter(data, Potential.contaminant != "+")
    
  } else if(length(option)==1 & option[1] == "reverse") { #4844
    
    data <- dplyr::filter(data, Reverse != "+")
    
  } else if(length(option)==1 & option[1] == "identified") { #4877
    
    data <- dplyr::filter(data, Only.identified.by.site != "+")
    
  } else if(length(option)==2 & 
            (option[1] == "potential" & option[2] == "reverse")) { #4793
    
    data <- dplyr::filter(data, Potential.contaminant != "+" &
                            Reverse != "+")
  } else if(length(option)==2 & 
            (option[1] == "potential" & option[2] == "identified")) { # 4826
    
    data <- dplyr::filter(data, Potential.contaminant != "+" &
                            Only.identified.by.site != "+")
  } else if(length(option)==2 & 
            (option[1] == "reverse" & option[2] == "identified")) { #4829
    
    data <- dplyr::filter(data, Reverse != "+" &
                            Only.identified.by.site != "+")
  } else if(length(option)==3 & 
            (option[1] == "potential" &
             option[2] == "reverse" & option[3] == "identified")) { #4805
    
    data <- dplyr::filter(data, Potential.contaminant != "+" &
                            Only.identified.by.site != "+" & Reverse != "+")
    
  } else if(length(option)==0) {
    data <- data
  } 
  return(data)
}



get_main_data_LiB <- function(data, file_type) {
  rv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")
  
  ID <- as.character(data$Majority.protein.IDs)
  Protein_count <- as.character(data$Peptide.counts..all.)
  name <- as.character(data$Gene.names)
  
  multi_pos <- grep(";",Protein_count,fixed = T)
  ID_multi <- ID[multi_pos] 
  Protein_count_multi <- Protein_count[multi_pos]
  
  pt_multi <- c()
  for(i in 1:length(ID_multi)){
    pt_id <- ID_multi[i]
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
        data_pos <- grep(ID_multi[i],data$Majority.protein.IDs)
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
      gn_multi <- data.frame(ID=pt_id,name=pt_gn)
    } else{
      pt_gn <- pt_db$Gene.symbol
      gn_multi <- data.frame(ID=pt_db$Entry,name=pt_gn)
    }
    
    pt_multi <- rbind(pt_multi,gn_multi)
  }
  
  ID[multi_pos] <- as.character(pt_multi$ID)
  name[multi_pos] <- as.character(pt_multi$name)
  
  cn_data <- colnames(data)
  
  if(file_type=="LFQ"){
    main_pos <- grep("LFQ",cn_data)
  } else {
    main_pos <- grep("iBAQ.",cn_data)
  }
  
  main <- data[,main_pos]
  
  unique_pos <- grep("Unique.peptides",cn_data)
  Unique_peptides <- data[,unique_pos]
  Sequence_coverage <- data$Sequence.coverage....
  
  main_df <- cbind(ID,name,main,Unique_peptides,Sequence_coverage)
  main_df <- make_unique(main_df, "name", "ID", delim=";")
  return(main_df)
}



get_main_data_T <- function(data,normalization) {
  rv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")
  
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
    main_data <- main_data[,-nor_pos]
  } else{
    origin_pos <- grep("Abundance", cn)
    nor_pos <- grep("Normalized",cn)
    nor_pos <- setdiff(origin_pos, nor_pos)
    main_pos <- c(1,nor_pos)
    main_data <- main_data[,main_pos]
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


make_case_samples_LiB <- function(data) {
  
  # file_type <- paste(file_type, ".")
  temp_col <- colnames(data)
  temp_col <- temp_col[grep(file_type, temp_col)]
  # file_type <- as.character(paste0(file_type,".intensity."))
  # temp <- data.frame(do.call('rbind', strsplit(as.character(temp_col),
  #                                              split = file_type, fixed = TRUE)))
  # col <- data.frame(id=temp[,2])
  # 
  # col1 <- col[1:(nrow(col)/2),]
  # col2 <- col[((nrow(col)/2)+1):nrow(col),]
  # 
  # col <- union(col1, col2)
  # col <- paste0(file_type, col)
  
  return(temp_col)
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
  
  nor_pos <- grep("Normalized",cn)
  main_pos <- c(1,nor_pos)
  main_data <- main_data[,main_pos]
  # if(normalization == "T"){
  #   nor_pos <- grep("Normalized",cn)
  #   main_pos <- c(1,nor_pos)
  #   main_data <- main_data[,main_pos]
  # } else{
  #   nor_pos <- grep("Normalized",cn)
  #   main_data <- main_data[,-nor_pos]
  # }
  
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
  return(col)
}


use_transformation_option <- function(data, samples) {
  for(i in 1:nrow(data)) {
    data[i,c(samples)] <- apply(data[i,c(samples)],1,function(x){2^(x)})
  }
  return(data)
}

use_valid_option <- function(data, case, control, valid_num) {
  case_min <- length(case) * valid_num
  control_min <- length(control) * valid_num
  cols <- colnames(data)
  
  for(i in 1:nrow(data)) {
    case_sums <- rowSums(is.finite(as.matrix(data[i,c(case)])))
    control_sums <- rowSums(is.finite(as.matrix(data[i,c(control)])))
    
    if((case_sums < case_min) || (control_sums < control_min)){
      data[i, "Keep_or_not"] <- c("Not")
    } else if((case_sums >= case_min) || (control_sums >= control_min)){
      data[i, "Keep_or_not"] <- c("Keep")
    }
  }
  
  data <- subset(data, Keep_or_not=="Keep")
  data <- select(data, cols)
  
  return(data)
}

fit_data_into_DEP_format <- function(file_type, data){
  if(file_type=="TMT"){
    data$Protein.IDs <- as.character(data$Accession)
    data$Description <- gsub("Gene_Symbol=", "GN=", data$Description)
    temp <- data.frame(do.call("rbind",
                               strsplit(as.character(data$Description), split="GN=", fixed=T)))
    temp <- select(temp, X2)
    temp <- data.frame(do.call("rbind",
                               strsplit(as.character(temp$X2), split="PE=", fixed=T)))
    temp <- select(temp, X1)
    temp$X1 <- gsub(" [[:print:]]*", "", temp$X1)
    colnames(temp) <- c("Gene.names")
    data <- cbind(temp, data)
  } 
  
}



get_preprocessed_data <- function(file_type, data, option, case, control) {
  
  log_option <- option[[1]][1]
  valid_option <- option[[1]][2]
  imputation_option <- option[[1]][3]
  normalization_option <- option[[1]][4]
  samples <- c(case, control)
  
  if(log_option=="none") {
    data <- use_transformation_option(data, samples)
  } else {
    data <- data
  }
  
  cols <- paste0(samples, collapse="|")
  cols <- grep(cols, colnames(data))
  # data_se <- make_se_parse(data, cols)
  case_design <- data.frame(label=c(case), condition=rep("Case", length(case)),
                            replicate=c(1:length(case)))
  control_design <- data.frame(label=c(control), condition=rep("Control", length(control)),
                               replicate=c(1:length(control)))
  
  exp_design <- rbind(case_design, control_design)
  exp_design$label <- as.character(exp_design$label)
  exp_design$condition <- as.character(exp_design$condition)
  exp_design$replicate <- as.character(exp_design$replicate)
  
  data_se <- make_se(data, cols, exp_design)
  
  return(data_se)
  # plot_frequency(data_se)
  
  
  
  # data <- switch(valid_option,
  #                   "30" = use_valid_option(data, case, control, 0.3),
  #                   "50" = use_valid_option(data, case, control, 0.5),
  #                   "70" = use_valid_option(data, case, control, 0.7),
  #                   "100" = use_valid_option(data, case, control, 1))
  # 
  # data <- switch(imputation_option,
  #                   "normal_distribution" = use_imputation_option(data, "normal_distribution"),
  #                   "constant" = use_imputation_option(data, "constant"),
  #                   "nan" = use_imputation_option(data, "nan"),
  #                   "none" = use_imputation_option(data, "none"))
  
  # list("Normal distribution" = "normal_distribution", "Constant" = "constant", 
  # "NaN"="nan", "None"="none")
  
}
# case <- c("LFQ.intensity.Total_309B", "LFQ.intensity.Total_445B", "LFQ.intensity.Total_555B",
# "LFQ.intensity.Total_588B", "LFQ.intensity.Total_636B", "LFQ.intensity.Total_667B",
# "LFQ.intensity.Total_764B", "LFQ.intensity.Total_741B", "LFQ.intensity.Total_876B",
# "LFQ.intensity.Total_883B")
# control <- c("LFQ.intensity.Total_309M", "LFQ.intensity.Total_445M", "LFQ.intensity.Total_555M",
# "LFQ.intensity.Total_588M", "LFQ.intensity.Total_636M", "LFQ.intensity.Total_667M",
# "LFQ.intensity.Total_741M", "LFQ.intensity.Total_764M", "LFQ.intensity.Total_876M",
# "LFQ.intensity.Total_883M")
# 
# 
# 
# # data <- tmt
# data$Protein.IDs <- as.character(data$Accession)
# data$Description <- gsub("Gene_Symbol=", "GN=", data$Description)
# temp <- data.frame(do.call("rbind",
#                            strsplit(as.character(data$Description), split="GN=", fixed=T)))
# temp <- select(temp, X2)
# temp <- data.frame(do.call("rbind",
#                            strsplit(as.character(temp$X2), split="PE=", fixed=T)))
# temp <- select(temp, X1)
# temp$X1 <- gsub(" [[:print:]]*", "", temp$X1)
# colnames(temp) <- c("Gene.names")
# data <- cbind(temp, data)
# 
# 
# data <- main_df
# data$name %>% duplicated() %>% any()
# data %>% group_by(name) %>% summarise(frequency = n()) %>%
#   arrange(desc(frequency)) %>% filter(frequency > 1)
# data_unique <- make_unique(data, "name", "ID", delim=";")
# data_unique$name %>% duplicated() %>% any()
# 
# 
# 
# 
# 
# 
# 
# data_filt <- filter_missval(data_se, thr=0)
# plot_numbers(data_filt)
# 
# plot_coverage(data_filt)
# 
# data_norm <- normalize_vsn(data_filt)
# plot_normalization(data_filt, data_norm)
# plot_missval(data_filt)
# plot_detect(data_filt)
# 
# data_imp<-impute(data_norm, fun="MinProb", q=0.01)
# data_imp_man <- impute(data_norm, fun="man", shift=1.8, scale=0.3)
# data_imp_knn <- impute(data_norm, fun="knn", rowmax=0.9)
# plot_imputation(data_norm, data_imp)
# 
# 
# data_diff <- test_diff(data_imp, type="control", control="C")
