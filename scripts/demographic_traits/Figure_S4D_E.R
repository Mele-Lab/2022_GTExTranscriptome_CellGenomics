#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(reshape2)
library(ggplot2)
library(ggtext)
library(ggpubr)

# metadata ----
load("00.metadata.RData")

# differential splicing analysis ----
dsa <- readRDS("00.differential_splicing_analysis.rds")

# Figure S4D ----
# total unique number of DSEs per tissue
get_DSEs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05,])
  }
}
DSEs <- lapply(traits, function(trait) lapply(tissues, function(tissue) get_DSEs(tissue, trait)))
names(DSEs) <- traits
for(trait in traits){names(DSEs[[trait]]) <- tissues}

DSEs_per_trait_and_tissue <- unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, sum(dsa[[tissue]][[trait]]$adj.P.Val < 0.05)))))
DSEs_per_tissue <- sapply(tissues, function(tissue) length(unique(unlist(lapply(traits, function(trait) DSEs[[trait]][[tissue]])))[!is.na(unique(unlist(lapply(traits, function(trait) DSEs[[trait]][[tissue]]))))]))
proportion_of_DSEs_per_trait <- 100*DSEs_per_trait_and_tissue/rep(DSEs_per_tissue,4)

# average gene expression variation explained per trait and tissue --
data <- cbind.data.frame("Tissue" = rep(tissue_info$tissue_abbrv, 4),
                         "Trait" = unlist(lapply(traits, function(trait) rep(trait, length(tissues)))),
                         "DSEs" = unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, sum(dsa[[tissue]][[trait]]$adj.P.Val < 0.05))))),
                         "proportion_of_DSEs" <- proportion_of_DSEs_per_trait,
                         "R2" = unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, mean(dsa[[tissue]][[trait]]$R2, na.rm = T)))))
)
data$Tissue <- factor(data$Tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$Trait <- factor(data$Trait, levels = traits, order = T)

pdf(paste0(plot_path, "Figure_S4D.alternative_splicing_variation_explained.pdf"),
    width = 4, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DSEs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean alternative splicing variation explained (%)") +
  xlab("") +
  theme(#axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

pdf(paste0(plot_path, "Figure_S4D.alternative_splicing_variation_explained.legend.pdf"),
    width = 5, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DSEs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean alternative splicing variation explained (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

# Figure S4E ----
# events DS with each trait in each tissue
get_dse <- function(tissue,trait){
  if(length(dsa[[tissue]][[trait]])==1){
    return(NA)
  }else{
    ds_events <- rownames(dsa[[tissue]][[trait]])[ dsa[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(ds_events)
  }
}

DSEs <- lapply(tissues, function(tissue)
  lapply(traits, function(trait) 
    get_dse(tissue,trait)
  )
)
names(DSEs) <- tissues
for(tissue in tissues){names(DSEs[[tissue]]) <- traits}

# splicing variation explained by each trait in each tissue for each DEG
get_value <- function(tissue){
  events <- unique(unlist(DSEs[[tissue]]))
  if(tissue %in% sex_tissues){
    data <- t(sapply(events, function(event) sapply(traits[-2], function(trait) dsa[[tissue]][[trait]][event, "R2"]))) 
  }else{
    data <- t(sapply(events, function(event) sapply(traits, function(trait) dsa[[tissue]][[trait]][event, "R2"])))  
  }
  d <- melt(data)
  d$tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"])
  d <- d[,c(4,1,2,3)]
  colnames(d) <- c("tissue", "event", "variable", "value")
  return(d)
}

data <- do.call(rbind.data.frame, lapply(tissues, function(tissue) get_value(tissue)))
data$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$chr <- factor(data$chr, levels = c("autosomes", "chrX", "chrY"), order = T)
data$variable <- gsub("_abs", "", data$variable)
data$variable <- factor(data$variable, levels = rev(traits[c(2,3,1,4)]), order = T)
data <- data[!is.na(data$value),] 
data$ensembl_id <- sapply(as.character(data$event), function(e) unlist(strsplit(e, split = ";"))[[1]])
data$gene_name <- sapply(data$ensembl_id, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name"])

pdf(paste0(plot_path, "Figure_S4E.alternative_splicing_variation_explained.box_plots.pdf"),
    width = 4, height = 9)
ggplot(data,
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Alternative splicing variation explained (%)") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) 
dev.off()


