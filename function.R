file_input_test <- function(file, file_type) {
  switch(file_type,
         "TMT" = 
           {cols <- colnames(file)[grep("Abundance", colnames(file))]},
         "iBAQ" = 
           {cols <- colnames(file)[grep("iBAQ", colnames(file))]},
         "LFQ" = 
           {cols <- colnames(file)[grep("LFQ", colnames(file))]})
  
  # if(file_type=="TMT"){
  #   cols <- colnames(file)[grep("Abundance", colnames(file))]
  # } else {
  #   cols <- colnames(file)[grep(file_type, colnames(file))]
  # }
  
  return(cols)
}



filter_with_option <- function(option, data) {
  data <- dplyr::filter(data, Peptides != 0)#4903
  data <- dplyr::filter(data, Intensity != 0)#4898
  
  option <- paste0(option, collapse=",")
  switch(option,
         "potential" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+")},
         "reverse" = 
           {data <- dplyr::filter(data, Reverse != "+")},
         "identified" = 
           {data <- dplyr::filter(data, Only.identified.by.site != "+")},
         "potential,reverse" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Reverse != "+")},
         "potential,identified" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Only.identified.by.site != "+")},
         "reverse,identified" = 
           {data <- dplyr::filter(data, Reverse != "+" &
                                    Only.identified.by.site != "+")},
         "potential,reverse,identified" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Only.identified.by.site != "+" & Reverse != "+")},
         {data <- data}
  )
  return(data)
}



get_main_data_LiB <- function(data, file_type) {
  rv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")
  
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
  
  if(file_type=="LFQ"){
    main_pos <- grep("LFQ",cn_data)
  } else {
    main_pos <- grep("iBAQ.",cn_data)
  }
  
  main_data <- data[,main_pos]

  unique_pos <- grep("Unique.peptides",cn_data)
  Unique_peptides <- data[,unique_pos]
  Sequence_coverage <- data$Sequence.coverage....
  
  main_df <- cbind(ProteinID,GeneName,main_data,Unique_peptides,Sequence_coverage)
  
  isDup_Gene <- main_df$GeneName %>% duplicated() %>% any()
  if(isDup_Gene==T){
    data_unique <- make_unique(main_df, "GeneName", "ProteinID")
  } else{
    data_col <- colnames(data)
    data_pos <- which("ProteinID"==data_col | "GeneName"==data_col)
    colnames(data)[data_pos] <- c("ID","name")
    data_unique <- data
  }
  
  isDup_name <- data$name %>% duplicated() %>% any()
  if(isDup_name == F){
    data_unique <- data.frame(ID=data_unique$ID,name=data_unique$name,data_unique[,-c(1:2,(ncol(data_unique)-1):ncol(data_unique))])
    return(data_unique)
  }
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
    main_pos <- c(1,nor_pos)
    main_data <- main_data[,main_pos]
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
  
  isDup_Gene <- main_df$Gene.symbol %>% duplicated() %>% any()
  if(isDup_Gene==T){
    print("true")
    data_unique <- make_unique(main_df, "Gene.symbol", "Accession")
  } else{
    print("false")
    data_col <- colnames(data)
    data_pos <- which("Accession"==data_col | "Gene.symbol"==data_col)
    colnames(data)[data_pos] <- c("ID","name")
    data_unique <- data
  }
  
  isDup_name <- data$name %>% duplicated() %>% any()
  if(isDup_name==F){
    data_unique <- data.frame(ID=data_unique$ID,name=data_unique$name,data_unique[,-c(1:2,(ncol(data_unique)-1):ncol(data_unique))])
    return(data_unique)
  } 
}

make_case_samples_LiB <- function(data, file_type) {
  temp_col <- colnames(data)
  temp_col <- temp_col[grep(file_type, temp_col)]
  
  return(temp_col)
}

make_case_samples_T <- function(data) {
  data <- data[,c(3:ncol(data))]
  temp_col <- colnames(data)
  
  temp <- c()
  TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  for(i in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[i],temp_col)],TMT_marker[i]) 
    temp <- c(temp,substr(temp_col[i],col_pos[1,1],(nchar(temp_col)+1)))
  }
  
  col <- data.frame(id=temp)
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  
  return(col)
}

make_case_samples_diffT <- function(data,normalization) {
  cn <- colnames(data)
  data$ID <- as.character(data$ID)
  data_class <- sapply(data,class)
  rm_pos <- grep("factor",data_class) 
  rm_pos <- c(rm_pos,grep("Ratio",cn))
  
  data <- data[,-rm_pos]
  cn <- cn[-rm_pos]
  
  main_pos <- grep("ID",cn)
  
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
  
  main_data <- na.omit(main_data)
  main_data <- main_data[,c(2:ncol(main_data))]
  temp_col <- colnames(main_data)
  
  temp <- c()
  for(j in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[j],temp_col)],TMT_marker[j]) 
    temp <- c(temp,substr(temp_col[j],col_pos[1,1],(nchar(temp_col)+1)))
  }
  
  col <- data.frame(ID=temp)
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  
  return(col)
}

make_expDesignData <- function(case,ctrl) {
  label <- c(case,ctrl)
  
  condition <- make_condition(case,ctrl)
  
  # condition <- c()
  # condition[1:length(case)] <- group_name[1]
  # condition[(length(case)+1):(length(case)+length(ctrl))] <- group_name[2]
  
  replicate_case <- c()
  for(i in 1:length(case)){
    replicate_case[i] <- i
  }
  
  replicate_ctrl <- c()
  for(j in 1:length(ctrl)){
    replicate_ctrl[j] <- j
  }
  
  replicate <- c(replicate_case,replicate_ctrl)
  
  exp_design <- data.frame(label=label, condition=condition, replicate=replicate)
  
  return(exp_design)
}


make_condition <- function(case,ctrl){
  label_case <- strsplit(case,".",fixed = T)
  label_ctrl <- strsplit(ctrl,".",fixed = T)
  
  case_name <- c()
  ctrl_name <- c()
  for(i in 1:length(label_case)){
    list_case <- label_case[[i]]
    list_ctrl <- label_ctrl[[i]]
    case_name <- c(case_name,list_case[length(list_case)])
    ctrl_name <- c(ctrl_name,list_ctrl[length(list_ctrl)])
  }
  
  case_name <- gsub("\\d+","",case_name)
  ctrl_name <- gsub("\\d+","",ctrl_name)
  
  condition <- c(case_name,ctrl_name)
  
  return(condition)
}


make_summarizedData <- function(main,file_type,design){
  if(file_type == "TMT"){
    data_col <- grep("Sample",colnames(main))
  } else if(file_type == "LFQ"){
    data_col <- grep("LFQ",colnames(main))
    main <- main[,c(1:(2+length(design$label)))]
  } else{
    data_col <- grep("iBAQ",colnames(main))
    main <- main[,c(1:(2+length(design$label)))]
  }
  design$label <- as.character(design$label)
  design$condition <- as.character(design$condition)
  design$replicate <- as.numeric(design$replicate)
  summary_data <- make_se(main,data_col,design)
  return(summary_data)
}

# use_transformation_option <- function(data, case, control, options) {
#   samples <- c(case, control)
# 
#   if(options=="none") {
#     for(i in 1:nrow(data)) {
#       data[i,c(samples)] <- apply(data[i,c(samples)],1,function(x){2^(x)})
#     }
#   } else {
#     data <- data
#   }
# 
#   return(data)
# }

use_valid_option <- function(data, case, control, option) {
  option <- as.numeric(as.character(option))
  valid_num <- length(case) * option
  
  data_filt <- filter_missval(data, thr=valid_num)
  # res <- assay(data_filt)
  # res <- round(res,digits = 2)
  # data_filt <- SummarizedExperiment(assays = list(as.data.frame(res)), rowData = rowData(data_filt), colData = colData(data_filt))
  
  return(data_filt)
}

use_normalization_option <- function(data, option) {
  if(option=="Yes") {
    data_norm <- normalize_vsn(data)
  }
  # } 
  # else if(option == "Z-score"){
  #   val_z <- assay(data)
  #   res <- c()
  #   for(i in 1:nrow(val_z)){
  #     tmp <- as.numeric(val_z[i,])
  #     res <- rbind(res,(tmp-mean(tmp))/sd(tmp))
  #   }
  #   data_norm <- SummarizedExperiment(assays = list(as.data.frame(res)), rowData = rowData(data), colData = colData(data))
  # }
  else {
    data_norm <- data
  }
  return(data_norm)
}


use_imputation_option <- function(data, case, control, option) {
  switch(option,
         "QRILC" =
           {data_imp <- impute(data, fun="QRILC")},
         "MinProb" =
           {data_imp <- impute(data, fun="MinProb", q=0.01)},
         "man" =
           {data_imp <- impute(data, fun = "man", shift = 1.8, scale = 0.3)},
         "zero" =
           {data_imp <- impute(data, fun="zero")}
  )

  return(data_imp)
}

test <- function(data,design,p_option,q_option){
  total_sample_num <- nrow(design)
  condition <- unique(design$condition)
  case_num <- length(grep(condition[1],design$condition))
  ctrl_num <- length(grep(condition[2],design$condition))
  case_dt <- data[,c(1:case_num)]
  ctrl_dt <- data[,c((case_num+1):(case_num+ctrl_num))]
  pval <- c()
  if(p_option == "T.Test"){
    case_df <- t(case_dt)
    case_df <- as.data.frame(case_df)
    ctrl_df <- t(ctrl_dt)
    ctrl_df <- as.data.frame(ctrl_df)
    pval <- mapply(function(x,y){t.test(x,y,alternative = c("two.sided"))},case_df,ctrl_df)
    pval <- t(pval)
    pval <- as.data.frame(pval)
    pval <- unlist(pval$p.value)
  }
  else if(p_option == "Wilcoxon-Ranksum"){
    case_df <- t(case_dt)
    case_df <- as.data.frame(case_df)
    ctrl_df <- t(ctrl_dt)
    ctrl_df <- as.data.frame(ctrl_df)
    pval <- mapply(function(x,y){wilcox.test(x,y,alternative = c("two.sided"))},case_df,ctrl_df)
    pval <- t(pval)
    pval <- as.data.frame(pval)
    pval <- unlist(pval$p.value)
  }
  else{
    #edgeR
    condition <- design[,c(2,3)]
    condition$condition <- as.factor(condition$condition)
    design<-model.matrix(~condition,data=condition)
    
    edr<-DGEList(counts=data, genes=rownames(data))
    noredr<-calcNormFactors(edr, method="TMM")
    
    disedr<-estimateDisp(noredr, design)
    fitq<-glmQLFit(disedr,design)
    glf<-glmQLFTest(fitq)
    res_edgeR<-glf$table
    pval <- res_edgeR$PValue
  }
  
  qval <- p.adjust(pval,method = q_option)
  df <- data.frame(pval=pval,padj=qval)
  return(df)
}

change_Sig <- function(data_rejection,pvalue,log2fc){
  dep_rowData <- rowData(data_rejection)
  pval_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
  log2fc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
  sig_pos <- grep("significant",colnames(dep_rowData),fixed = T)
  
  pval_df <- as.numeric(dep_rowData[,pval_pos])
  log2fc_df <- abs(as.numeric(dep_rowData[,log2fc_pos]))
  res_pval <- as.logical(pval_df < pvalue)
  res_lfc <-  as.logical(log2fc_df > log2fc)
  TT_pos <- which(res_pval=="TRUE" & res_lfc=="TRUE")
  
  res <- c()
  res[1:nrow(dep_rowData)] <- "FALSE"
  res[TT_pos] <- "TRUE"
  res <- as.logical(res)
  
  dep_rowData[,sig_pos] <- res
  dep_rowData$name <- as.character(dep_rowData$name)
  return(dep_rowData)
}

gsa <- function(data,set,tool){
  data <- data[,c(1,7)]
  colnames(data) <- c("genes","log2fc")
  if(set == "caseup"){
    genes <- data$genes[data$log2fc > 0]
  } else{
    genes <- data$genes[data$log2fc < 0]
  }
  
  if(tool == "enrichR"){
    res_gsa <- enrichR(genes)
    return(res_gsa)
  } else{
    
  }
}

enrichR <- function(genes){
  dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human") #, "WikiPathways_2016", "Reactome_2016"
  res_er<-enrichr(genes,dbs)
  return(res_er)
}