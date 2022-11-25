#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggtext)
library(ggpubr)

# metadata ----
load("00.metadata.RData")
autosomal_genes <- gene_annotation[!gene_annotation$chr %in% c("chrX", "chrY"), "ensembl.id"]
X_genes <- gene_annotation[gene_annotation$chr == "chrX", "ensembl.id"]
Y_genes <- gene_annotation[gene_annotation$chr == "chrY", "ensembl.id"]

# differential expression analysis ----
dea <- readRDS("00.differential_expression_analysis.rds")

# genes DE with each trait in each tissue
get_deg <- function(tissue,trait){
  if(length(dea[[tissue]][[trait]])==1){
    return(NA)
  }else{
    de_genes <- rownames(dea[[tissue]][[trait]])[ dea[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(de_genes)
  }
}

DEGs <- lapply(tissues, function(tissue)
  lapply(traits, function(trait) 
    get_deg(tissue,trait)
  )
)
names(DEGs) <- tissues
for(tissue in tissues){names(DEGs[[tissue]]) <- traits}

# expression variation explained by each trait in each tissue for each DEG
get_value <- function(tissue){
  genes <- unique(unlist(DEGs[[tissue]]))
  if(tissue %in% sex_tissues){
    data <- t(sapply(genes, function(gene) sapply(traits[-2], function(trait) dea[[tissue]][[trait]][gene, "R2"]))) 
  }else{
    data <- t(sapply(genes, function(gene) sapply(traits, function(trait) dea[[tissue]][[trait]][gene, "R2"])))  
  }
  d <- melt(data)
  d$tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"])
  d <- d[,c(4,1,2,3)]
  colnames(d) <- c("tissue", "gene", "variable", "value")
  d$chr <- sapply(d$gene, function(gene)  ifelse(gene %in% autosomal_genes,
                                                 "autosomes",
                                                 ifelse(gene %in% X_genes,
                                                        "chrX",
                                                        ifelse(gene %in% Y_genes,
                                                               "chrY",
                                                               NA))) )
  return(d)
}

data <- do.call(rbind.data.frame, lapply(tissues, function(tissue) get_value(tissue)))
data$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$chr <- factor(data$chr, levels = c("autosomes", "chrX", "chrY"), order = T)
data$variable <- gsub("_abs", "", data$variable)
data$variable <- factor(data$variable, levels = rev(traits[c(2,3,1,4)]), order = T)
data <- data[!is.na(data$value),] 

p1 <- ggplot(data[data$tissue %in% tissue_info[1:21, "tissue_abbrv"] &
                    data$chr == "autosomes",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

p2 <- ggplot(data[data$tissue %in% tissue_info[23:46, "tissue_abbrv"] &
                    data$chr == "autosomes",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,10))

p3 <- ggplot(data[data$tissue %in% tissue_info[1:21, "tissue_abbrv"] &
                    data$chr == "chrX",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

p4 <- ggplot(data[data$tissue %in% tissue_info[22:46, "tissue_abbrv"] &
                    data$chr == "chrX",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

p5 <- ggplot(data[data$tissue %in% tissue_info[1:21, "tissue_abbrv"] &
                    data$chr == "chrY",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats_chrY, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

p6 <- ggplot(data[data$tissue %in% tissue_info[23:46, "tissue_abbrv"] &
                    data$chr == "chrY",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~chr, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats_chrY, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))


pdf(paste0(plot_path, "Figure_S1D.gene_expression_variation_explained.pdf"),
    width = 12, height = 12)
ggarrange(p1, p3, p5, p2, p4, p6, nrow = 2, ncol = 3, heights = c(0.45, 0.55))
dev.off()
