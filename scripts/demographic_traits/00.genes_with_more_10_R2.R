#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(ggpubr)

# metadata ----
load("00.metadata.RData")

# differential expression analysis ----
dea <- readRDS("differential_expression_analysis.rds")

# DEGs woth more than 10% of residual expression variation explained
dea_10 <- list()
for(tissue in tissues){
  dea_10[[tissue]] <- list()
  for(trait in traits){
    if(tissue %in% sex_tissues & trait=="Sex"){
      dea_10[[tissue]][[trait]] <-  NA
    }else{
      dea_10[[tissue]][[trait]] <-  dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05 & dea[[tissue]][[trait]]$R2 > 10,]  
    }
  }
}

# total unique number of DEGs per tissue
get_DEGs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(dea_10[[tissue]][[trait]][dea_10[[tissue]][[trait]]$adj.P.Val < 0.05,])
  }
}

DEGs <- lapply(traits, function(trait) lapply(tissues, function(tissue) get_DEGs(tissue, trait)))
names(DEGs) <- traits
for(trait in traits){names(DEGs[[trait]]) <- tissues}


genes <- unique(unlist(lapply(traits, function(trait) unlist(DEGs[[trait]]))))
length(genes)
