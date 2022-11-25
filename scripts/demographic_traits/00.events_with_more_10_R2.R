#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)


# metadata ----
load("00.metadata.RData")

# differential splicing analysis ----
dsa <- readRDS("differential_splicing_analysis.rds")

# DSEs with more than 10% of residual splicing variation explained
dsa_10 <- list()
for(tissue in tissues){
  dsa_10[[tissue]] <- list()
  for(trait in traits){
    if(tissue %in% sex_tissues & trait=="Sex"){
      dsa_10[[tissue]][[trait]] <-  NA
    }else{
      dsa_10[[tissue]][[trait]] <-  dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05 & dsa[[tissue]][[trait]]$R2 > 10,]  
    }
  }
}

# total unique number of DEGs per tissue
get_DSEs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(dsa_10[[tissue]][[trait]][dsa_10[[tissue]][[trait]]$adj.P.Val < 0.05,])
  }
}

DSEs <- lapply(traits, function(trait) lapply(tissues, function(tissue) get_DSEs(tissue, trait)))
names(DSEs) <- traits
for(trait in traits){names(DSEs[[trait]]) <- tissues}

events <- unique(unlist(lapply(traits, function(trait) unlist(DSEs[[trait]]))))
length(events)
