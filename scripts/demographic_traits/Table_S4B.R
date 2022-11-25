#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
table_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/tables/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}
if(!dir.exists(table_path)){dir.create(table_path, recursive = T)}

# libraries ----
library(rcompanion)
library(gtools)

# metadata ----
load("00.metadata.RData")

# differential splicing tables ----
dsa <- readRDS("differential_splicing_analysis.rds")
for(tissue in tissues){
  if(tissue %in% sex_tissues){
    for(trait in traits[-2]){
      dsa[[tissue]][[trait]]$type <- sapply(rownames(dsa[[tissue]][[trait]]), function(e) unlist(strsplit(unlist(strsplit(e, split = ";"))[[2]], split = ":"))[[1]])
      dsa[[tissue]][[trait]]$iso_id <- paste0(dsa[[tissue]][[trait]]$spliced_in, "-", dsa[[tissue]][[trait]]$spliced_out)
      dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][sample(1:nrow(dsa[[tissue]][[trait]])),]
      dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][!duplicated(dsa[[tissue]][[trait]]$iso_id),]
    } 
  }else{
    for(trait in traits){
      dsa[[tissue]][[trait]]$type <- sapply(rownames(dsa[[tissue]][[trait]]), function(e) unlist(strsplit(unlist(strsplit(e, split = ";"))[[2]], split = ":"))[[1]])
      dsa[[tissue]][[trait]]$iso_id <- paste0(dsa[[tissue]][[trait]]$spliced_in, "-", dsa[[tissue]][[trait]]$spliced_out)
      dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][sample(1:nrow(dsa[[tissue]][[trait]])),]
      dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][!duplicated(dsa[[tissue]][[trait]]$iso_id),]
    }
  }
}

# binomial test ----
data <- lapply(splicing_events, function(event)
  t(sapply(tissues, function(tissue) table(dsa[[tissue]][["Age"]][dsa[[tissue]][["Age"]]$type == event, "biotype"])[c("NC-PC", "PC-NC")]))
)
names(data) <- splicing_events
p_value <- matrix(data = unlist(lapply(splicing_events, function(event) 
  apply(data[[event]], 1, function(x) binom.test(x)$p.value)
)), 
nrow = length(tissues),
ncol = length(splicing_events), 
byrow = F)
adj_p_value <- matrix(data = p.adjust(unlist(lapply(splicing_events, function(event) 
  apply(data[[event]], 1, function(x) binom.test(x)$p.value)
)), method = "BH"), 
nrow = length(tissues),
ncol = length(splicing_events), 
byrow = F)
fdr <- adj_p_value
fdr[fdr>=0.05] <- NA
fdr <- -log10(fdr)
colnames(p_value) <- splicing_events
colnames(fdr) <- splicing_events
rownames(p_value) <- tissues
rownames(fdr) <- tissues

write.table(cbind.data.frame(p_value, fdr),
            paste0(table_path, "Table_S4B.tab"),
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")
