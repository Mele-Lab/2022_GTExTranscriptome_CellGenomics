#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
table_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/tables/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}
if(!dir.exists(table_path)){dir.create(table_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(ggtext)

# metadata ----
load("00.metadata.RData")

# DEGs tissue sharing --
tissue_sharing <- readRDS("genes_DE.tissue_sharing.rds")
for(trait in traits){
  tissue_sharing[[trait]]$trait <- rep(trait, nrow(tissue_sharing[[trait]]))
}

quantile(tissue_sharing$Age$n_DE, probs = seq(0,1,0.1))
quantile(tissue_sharing$Ancestry$n_DE, probs = seq(0,1,0.1))
quantile(tissue_sharing$Sex$n_DE, probs = seq(0,1,0.1))
quantile(tissue_sharing$BMI$n_DE, probs = seq(0,1,0.1))

# genes exclusively expressed in the tissues where they are DE --
data <- lapply(traits, function(trait) tissue_sharing[[trait]][tissue_sharing[[trait]]$n_DE==tissue_sharing[[trait]]$n_expressed, ])
names(data) <- traits
sum(sapply(traits, function(trait) nrow(data[[trait]])))
sapply(traits, function(trait) nrow(data[[trait]]))
length(unique(unlist(lapply(traits, function(trait) data[[trait]]$ensembl_id))))

# write table
for(trait in traits){
  write.table(data[[trait]][, c("ensembl_id", "gene_name", "tissues_DE")],
              paste0(table_path, "Table_S1E.", trait, ".tab"),
              col.names = T, row.names = T,
              quote = F,
              sep = "\t")  
}

nrow(tissue_sharing$BMI[tissue_sharing$BMI$n_DE>1,])
nrow(tissue_sharing$BMI[tissue_sharing$BMI$n_DE>3,])
tissue_sharing$BMI[tissue_sharing$BMI$gene_name=="LEP",]
tissue_sharing$BMI[tissue_sharing$BMI$gene_name=="AKAP1",]
