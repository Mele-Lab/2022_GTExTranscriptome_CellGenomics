#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# metadata ----
load("00.metadata.RData")
autosomes <- gene_annotation[!gene_annotation$chr %in% c("chrX", "chrY", "chrX"), "ensembl.id"]
sex_chr <- gene_annotation[gene_annotation$chr %in% c("chrX", "chrY", "chrX"), "ensembl.id"]

# differential expression analysis ----
dea <- readRDS("00.differential_expression_analysis.rds")

# number of differentially expressed genes (DEGs) per tissue and trait ---
number_of_DEGs <- t(sapply(tissues, function(tissue)
  sapply(traits, function(trait)
    ifelse(tissue %in% sex_tissues & trait=="Sex",
           NA,
           sum(dea[[tissue]][[trait]]$adj.P.Val < 0.05)
    )
  )
))

for(tissue in tissues){
  if(tissue %in% sex_tissues){
    for(trait in traits[-2]){
      dea[[tissue]][[trait]]$tissue <- rep(tissue, nrow(dea[[tissue]][[trait]]))
      dea[[tissue]][[trait]]$ensembl_id <- rownames(dea[[tissue]][[trait]])
    }
  }else{
    for(trait in traits){
      dea[[tissue]][[trait]]$tissue <- rep(tissue, nrow(dea[[tissue]][[trait]]))
      dea[[tissue]][[trait]]$ensembl_id <- rownames(dea[[tissue]][[trait]])
    } 
  }
}

data <- list()
trait_tissues <- apply(number_of_DEGs, 2, function(x) tissues[x>=100][!is.na(tissues[x>=100])])
names(trait_tissues) <- traits

for(trait in traits){
  if(trait == "Sex"){
    data[[trait]] <- do.call(rbind.data.frame, lapply(trait_tissues[[trait]], function(tissue) dea[[tissue]][[trait]][!is.na(dea[[tissue]][[trait]]$R2),]))  
    data[["Sex_autosomes"]] <- data[["Sex"]][data[["Sex"]]$ensembl_id %in% autosomes,]
    data[["Sex_sex_chr"]] <- data[["Sex"]][data[["Sex"]]$ensembl_id %in% sex_chr,]
  }else{
    data[[trait]] <- do.call(rbind.data.frame, lapply(trait_tissues[[trait]], function(tissue) dea[[tissue]][[trait]][!is.na(dea[[tissue]][[trait]]$R2),]))
  }
}

for(trait in names(data)){
  print(trait)
  data[[trait]] <- data[[trait]][order(data[[trait]]$R2, decreasing = T), c("tissue", "ensembl_id", "gene_name", "R2")]
  data[[trait]] <- data[[trait]][1:50,]
  data[[trait]]$R2 <- round(data[[trait]]$R2, 2)
}

for(trait in names(data)){
  write.table(data[[trait]],
              paste0("~/GTEx_v8/Raquel/Draft/Cell_submission/data/", trait, "Table_S1B.expression_hierpart_top50.tab"),
              col.names = T, row.names = F,
              quote = F, sep = "\t")
}

# differential splicing analysis ----
dsa <- readRDS("00.differential_splicing_analysis.rds")

# number of differentially spliced events (DSEs) per tissue and trait ---
number_of_DSEs <- t(sapply(tissues, function(tissue)
  sapply(traits, function(trait)
    ifelse(tissue %in% sex_tissues & trait=="Sex",
           NA,
           sum(dsa[[tissue]][[trait]]$adj.P.Val < 0.05)
    )
  )
))

for(tissue in tissues){
  if(tissue %in% sex_tissues){
    for(trait in traits[-2]){
      dsa[[tissue]][[trait]]$tissue <- rep(tissue, nrow(dsa[[tissue]][[trait]]))
      dsa[[tissue]][[trait]]$event_id <- rownames(dsa[[tissue]][[trait]])
    }
  }else{
    for(trait in traits){
      dsa[[tissue]][[trait]]$tissue <- rep(tissue, nrow(dsa[[tissue]][[trait]]))
      dsa[[tissue]][[trait]]$event_id <- rownames(dsa[[tissue]][[trait]])
    } 
  }
}

data <- list()
trait_tissues <- apply(number_of_DEGs, 2, function(x) tissues[x>=100][!is.na(tissues[x>=100])])
names(trait_tissues) <- traits

for(trait in traits){
  if(trait == "Sex"){
    data[[trait]] <- do.call(rbind.data.frame, lapply(trait_tissues[[trait]], function(tissue) dsa[[tissue]][[trait]][!is.na(dsa[[tissue]][[trait]]$R2),]))    
  }else{
    data[[trait]] <- do.call(rbind.data.frame, lapply(trait_tissues[[trait]], function(tissue) dsa[[tissue]][[trait]][!is.na(dsa[[tissue]][[trait]]$R2),]))
  }
}
for(trait in traits){
  data[[trait]] <- data[[trait]][order(data[[trait]]$R2, decreasing = T), c("tissue", "event_id", "gene_name", "R2")]
  data[[trait]] <- data[[trait]][1:50,]
  data[[trait]]$R2 <- round(data[[trait]]$R2, 2)
}

for(trait in traits){
  write.table(data[[trait]],
              paste0("~/GTEx_v8/Raquel/Draft/Cell_submission/data/", trait, "Table_S1B.splicing_hierpart_top50.tab"),
              col.names = T, row.names = F,
              quote = F, sep = "\t")
}
