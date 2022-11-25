#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(reshape2)
library(ggplot2)

# metadata ----
load("00.metadata.RData")
autosomal_genes <- gene_annotation[!gene_annotation$chr %in% c("chrX", "chrY"), "ensembl.id"]

# differential expression analysis ----
dea <- readRDS("00.differential_expression_analysis.rds")

# Figure S3C ----
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
    data <- as.data.frame(t(sapply(genes, function(gene) sapply(traits[-2], function(trait) dea[[tissue]][[trait]][gene, "R2"]))) )
  }else{
    data <- as.data.frame(t(sapply(genes, function(gene) sapply(traits, function(trait) dea[[tissue]][[trait]][gene, "R2"])))  )
  }
  data$n <- as.numeric(apply(data, 1, function(x) sum(!is.na(x))))
  if(tissue %in% sex_tissues){
    data$value <- as.numeric(apply(data[,1:3], 1, function(x) sum(x[!is.na(x)])))
  }else{
    data$value <- as.numeric(apply(data[,1:4], 1, function(x) sum(x[!is.na(x)])))
  }
  data$tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"])
  if(tissue %in% sex_tissues){
    data$Sex <- rep(NA, nrow(data))
    data <- data[,c(traits, "n", "tissue", "value")]
  }
  data$ensembl_id <- rownames(data)
  data <- data[data$ensembl_id %in% autosomal_genes,]
  rownames(data) <- NULL
  return(data)
}

data <- do.call(rbind.data.frame, lapply(tissues, function(tissue) get_value(tissue)))
data$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$n <- factor(data$n, levels = c(4:1), order = T)


n_cols <- c("#EA6B70",
            "#CA3C70",
            "#8F2390",
            "#00299C")
names(n_cols) <- c(1:4)

pdf(paste0(plot_path, "Figure_S3C.pdf"),
    width = 4, height = 9)
ggplot(data,
       aes(x = tissue,
           y = value,
           fill = n,
           col = n)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = n_cols) +
  scale_fill_manual(values = sapply(n_cols, function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) 
dev.off()

p2 <- ggplot(data,
             aes(x = tissue,
                 y = value,
                 fill = n,
                 col = n)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = n_cols) +
  scale_fill_manual(values = sapply(n_cols, function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

ggarrange(p1, p2)
ggplot(data,
       aes(x = tissue,
           y = value,
           fill = n,
           col = n)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = n_cols) +
  scale_fill_manual(values = sapply(n_cols, function(i) alpha(i, 0.5))) +
  ylab("Expression variation explained (%)") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(limits=c(60,80))
