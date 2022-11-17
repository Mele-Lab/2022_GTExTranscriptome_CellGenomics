#!/usr/bin/env Rscript

# Tissues ----
first_dir <- "~/Documents/mn4/Jose/"
tissues <- list.dirs(paste0(first_dir, "03_Models/Tissues/"), full.names = F)[-1]
sex_tissues <- c("Uterus", "Vagina", "Ovary", "Prostate", "Testis")

# Differential expression analyses: results tables ----
function_to_read_dea <- function(tissue){
  files <- list.files(paste0(first_dir, "03_Models/Tissues/", tissue, "/"), full.names = T)
  files <- files[grep("voom", files)]
  file <- files[grep("interaction", files, invert = T)]
  readRDS(file)
}
dea_res <- lapply(tissues, function(tissue) function_to_read_dea(tissue))
names(dea_res) <- tissues

# Hier.part ----
function_to_read_hier <- function(tissue){
  files <- list.files(paste0(first_dir, "03_Models/Tissues/", tissue, "/"), full.names = T)
  files <- files[grep("hier_part", files)]
  file <- files[grep("spliced", files, invert = T)]
  readRDS(file)
}
hier.part.exprs <- lapply(tissues, function(tissue) function_to_read_hier(tissue))
names(hier.part.exprs) <- tissues


for(tissue in tissues){
  print(tissue)
  names(dea_res[[tissue]])[!names(dea_res[[tissue]]) %in% c("Age", "BMI")] <- sapply(names(dea_res[[tissue]])[!names(dea_res[[tissue]]) %in% c("Age", "BMI")], function(trait) substr(trait, 1, nchar(trait)-1))
  for(trait in names(dea_res[[tissue]])){
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    }else{
      dea_res[[tissue]][[trait]]$R2 <- NA
      for(gene in rownames(dea_res[[tissue]][[trait]])){
        if(dea_res[[tissue]][[trait]][gene, "adj.P.Val"] < 0.05){
          dea_res[[tissue]][[trait]][gene, "R2"] <- hier.part.exprs[[tissue]][gene, paste0(trait, "_abs")]
        }
      }
    }
  }
}

#Add gene_name
gene_annotation <- read.delim(paste0(first_dir,"00_Data/gencode.v26.GRCh38.genes.bed"), header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

for(tissue in tissues){
  print(tissue)
  for(trait in names(dea_res[[tissue]])){
    if(tissue %in% sex_tissues & trait == "Sex"){
      print(paste0(tissue, ": ", trait))
      next
    }else{
      dea_res[[tissue]][[trait]]$gene_name <- sapply(rownames(dea_res[[tissue]][[trait]]), function(ensembl) gene_annotation$symbol[gene_annotation$gene==ensembl])
      dea_res[[tissue]][[trait]] <- dea_res[[tissue]][[trait]][, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "gene_name", "R2")]
    }
  }
}

saveRDS(dea_res, paste0(first_dir, "Tables/00.differential_expression_analysis.clinical_traits.rds"))
