#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggrepel)

# metadata ----
load("00.metadata.RData")

# DEGs tissue sharing ----
tissue_sharing <- readRDS("genes_DE.tissue_sharing.rds")
for(trait in traits){
  tissue_sharing[[trait]]$trait <- rep(trait, nrow(tissue_sharing[[trait]]))
}

# Table S3E and S3F
table_S3E <- read.delim("Table_S3E.tsv", skip = 2)
ancestry_genes <- unique(unlist(lapply(table_S3E$gene_name, function(i) unlist(strsplit(i, split = "/")))))
table_S3F <- read.delim("Table_S3F.tsv", skip = 2)
age_genes <- unique(unlist(lapply(table_S3F$gene_name, function(i) unlist(strsplit(i, split = "/")))))

# XCI escapee genes 
oliva_et_al <- read.csv("Oliva_et_al.table_S3.csv")
XCI_genes <- oliva_et_al[oliva_et_al$Reported.Escapee.==1, "ENSEMBL_gene_id"]
XCI_genes <- XCI_genes[XCI_genes %in% gene_annotation$ensembl.id] # # oliva et al. includie pseudoautosomal genes
xci_genes <- gene_annotation[gene_annotation$ensembl.id %in%  XCI_genes, "gene.name"]
xci_genes_highly_shared <- xci_genes[xci_genes %in% tissue_sharing$Sex[tissue_sharing$Sex$n_DE > 9, "gene_name"]]
Y_genes <- gene_annotation[gene_annotation$chr == "chrY", "gene.name"]
Y_genes_highly_shared <- Y_genes[Y_genes %in% tissue_sharing$Sex[tissue_sharing$Sex$n_DE > 9, "gene_name"]]


traits_cols <- c(traits_cols, "black")
names(traits_cols)[5] <- "gene"

df <- do.call(rbind.data.frame, lapply(traits, function(trait) tissue_sharing[[trait]][, c(1,2,3,8)]))
color_variable <- function(gene, trait){
  if(gene %in% ancestry_genes & trait == "Ancestry"){
    "gene"
  }else if(gene %in% age_genes & trait == "Age"){
    "gene"
  }else if(gene %in% c(xci_genes_highly_shared, Y_genes_highly_shared) & trait == "Sex"){
    "gene"
  }else{
    trait
  }
}
df$col <- sapply(1:nrow(df), function(row) color_variable(df[row, "gene_name"], df[row, "trait"]))
df$trait <- factor(df$trait, levels = c(traits, "gene"), order = T)


df <- rbind.data.frame(df[df$col != "gene",], df[df$col == "gene",])

ggplot(df,
       aes(x = trait,
           y = n_DE,
           col = col)) +
  geom_jitter() +
  scale_color_manual(values = traits_cols) +
  theme_bw() +
  ylab("Number of tissues") + xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = 10, lty = 2, size = 0.25) +
  geom_hline(yintercept = 1.5, lty = 2, size = 0.25) +
  geom_hline(yintercept = 5.5, lty = 2, size = 0.25) +
  geom_text_repel(data = rbind.data.frame(df[df$col == "gene" & df$trait %in% c("Ancestry", "Age"),],
                                          df[df$trait == "Sex" & df$gene_name %in% c("XIST", "UTY"),],
                                          df[df$trait == "BMI",][1:3,]), 
                  aes(x = trait, y = n_DE, label = gene_name), 
                  cex = 3.5, col = "black")

