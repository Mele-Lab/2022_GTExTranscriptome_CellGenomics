#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)

# metadata ----
load("00.metadata.RData")

# differential expression analysis ----
dea <- readRDS("differential_expression_analysis.rds")

# number of differentially expressed genes (DEGs) per tissue and trait ---
number_of_DEGs <- t(sapply(tissues, function(tissue)
  sapply(traits, function(trait)
    ifelse(tissue %in% sex_tissues & trait=="Sex",
           NA,
           sum(dea[[tissue]][[trait]]$adj.P.Val < 0.05)
    )
  )
))

# differential splicing analysis ----
dsa <- readRDS("differential_splicing_analysis.rds")

# number of differentially expressed genes (DEGs) per tissue and trait ---
number_of_DSEs <- t(sapply(tissues, function(tissue)
  sapply(traits, function(trait)
    ifelse(tissue %in% sex_tissues & trait=="Sex",
           NA,
           sum(dsa[[tissue]][[trait]]$adj.P.Val < 0.05)
    )
  )
))

#  subset autosomal genes
autosomal_genes <- gene_annotation[!gene_annotation$chr %in% c("chrX", "chrY", 'chrM'), "ensembl.id"]
for(tissue in tissues){
  for(trait in traits){
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    }else{
      dea[[tissue]][[trait]] <- dea[[tissue]][[trait]][rownames(dea[[tissue]][[trait]]) %in% autosomal_genes,]
    }
  }
}

for(tissue in tissues){
  for(trait in traits){
    if(tissue %in% sex_tissues & trait == "Sex"){
      next
    }else{
      dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$ensembl_id %in% autosomal_genes,]
    }
  }
}

# tissues with at lest 5 DEGs & DSEs per trait
tissue_traits <- lapply(traits, function(trait)
  sapply(tissues, function(tissue)
    ifelse(number_of_DEGs[tissue, trait]>=5 & number_of_DSEs[tissue, trait]>=5, tissue, NA)
  )[!is.na(sapply(tissues, function(tissue)
    ifelse(number_of_DEGs[tissue, trait]>=5 & number_of_DSEs[tissue, trait]>=5, tissue, NA)
  ))]
)
names(tissue_traits) <- traits

# Comparison of the contribution of each demographic trait to the gene expression and alternative splicing variation (S4F) ----
a <- round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2)
b <- round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2)
round(round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2)/round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2),2)
1 - round(round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2)/round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2),2)
mean(1 - round(round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2)/round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) median(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))),2),2))
mean((a-b)/a)

# Comparison of the relative contribution of each demographic trait to the total tissue expression and splicing variation explained (4H) ----
# proportion of total tissue expression variation explained by each trait 
get_tissue_expression_variation_explained <- function(tissue){
  print(tissue)
  if(tissue %in% sex_tissues){
    tissue_expression_variation_explained <- sapply(traits[-2], function(trait) sum(dea[[tissue]][[trait]]$R2[!is.na(dea[[tissue]][[trait]]$R2)]))/sum(sapply(traits[-2], function(trait) sum(dea[[tissue]][[trait]]$R2[!is.na(dea[[tissue]][[trait]]$R2)])))
    tissue_expression_variation_explained <- c(tissue_expression_variation_explained[c(1)], 0, tissue_expression_variation_explained[c(2,3)])
    names(tissue_expression_variation_explained)[2] <- "Sex"
  }else{
    tissue_expression_variation_explained <- sapply(traits, function(trait) sum(dea[[tissue]][[trait]]$R2[!is.na(dea[[tissue]][[trait]]$R2)]))/sum(sapply(traits, function(trait) sum(dea[[tissue]][[trait]]$R2[!is.na(dea[[tissue]][[trait]]$R2)])))
  }
  return(tissue_expression_variation_explained)   
}
tissue_expression_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues, function(tissue) get_tissue_expression_variation_explained(tissue)))
colnames(tissue_expression_variation_explained) <- traits
rownames(tissue_expression_variation_explained) <- tissues    

get_tissue_splicing_variation_explained <- function(tissue){
  print(tissue)
  if(tissue %in% sex_tissues){
    tissue_splicing_variation_explained <- sapply(traits[-2], function(trait) sum(dsa[[tissue]][[trait]]$R2[!is.na(dsa[[tissue]][[trait]]$R2)]))/sum(sapply(traits[-2], function(trait) sum(dsa[[tissue]][[trait]]$R2[!is.na(dsa[[tissue]][[trait]]$R2)])))
    tissue_splicing_variation_explained <- c(tissue_splicing_variation_explained[c(1)], 0, tissue_splicing_variation_explained[c(2,3)])
    names(tissue_splicing_variation_explained)[2] <- "Sex"
  }else{
    tissue_splicing_variation_explained <- sapply(traits, function(trait) sum(dsa[[tissue]][[trait]]$R2[!is.na(dsa[[tissue]][[trait]]$R2)]))/sum(sapply(traits, function(trait) sum(dsa[[tissue]][[trait]]$R2[!is.na(dsa[[tissue]][[trait]]$R2)])))
  }
  return(tissue_splicing_variation_explained)   
}
tissue_splicing_variation_explained <- do.call(rbind.data.frame,
                                                 lapply(tissues, function(tissue) get_tissue_splicing_variation_explained(tissue)))
colnames(tissue_splicing_variation_explained) <- traits
rownames(tissue_splicing_variation_explained) <- tissues    



a <- round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_expression_variation_explained[tissue, trait]))),2)
b <- round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_splicing_variation_explained[tissue, trait]))),2)
a
b
x <- (a-b)/a
mean(x[-1])
(b-a)/b
round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_splicing_variation_explained[tissue, trait]))),2)/round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_expression_variation_explained[tissue, trait]))),2)
round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_expression_variation_explained[tissue, trait]))),2)/round(sapply(traits, function(trait) mean(sapply(tissue_traits[[trait]], function(tissue) tissue_splicing_variation_explained[tissue, trait]))),2)






sd(sapply(tissue_traits[[trait]], function(tissue) median(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05, "R2"])))
