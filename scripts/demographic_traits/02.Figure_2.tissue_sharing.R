#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(ggtext)

# metadata ----
load("00.metadata.RData")


# Figure 2D ----
# DEGs tissue sharing --
tissue_sharing <- readRDS("genes_DE.tissue_sharing.rds")
for(trait in traits){
  tissue_sharing[[trait]]$trait <- rep(trait, nrow(tissue_sharing[[trait]]))
}

# Table S3E and S3F --
table_S3E <- read.delim("Table_S3E.tsv", skip = 2)
ancestry_genes <- unique(unlist(lapply(table_S3E$gene_name, function(i) unlist(strsplit(i, split = "/")))))
table_S3F <- read.delim("Table_S3F.tsv", skip = 2)
age_genes <- unique(unlist(lapply(table_S3F$gene_name, function(i) unlist(strsplit(i, split = "/")))))

# XCI escapee genes --
oliva_et_al <- read.csv("Oliva_et_al.table_S3.csv")
XCI_genes <- oliva_et_al[oliva_et_al$Reported.Escapee.==1, "ENSEMBL_gene_id"]
XCI_genes <- XCI_genes[XCI_genes %in% gene_annotation$ensembl.id] # # oliva et al. includie pseudoautosomal genes
xci_genes <- gene_annotation[gene_annotation$ensembl.id %in%  XCI_genes, "gene.name"]
xci_genes_highly_shared <- xci_genes[xci_genes %in% tissue_sharing$Sex[tissue_sharing$Sex$n_DE > 9, "gene_name"]]
Y_genes <- gene_annotation[gene_annotation$chr == "chrY", "gene.name"]
Y_genes_highly_shared <- Y_genes[Y_genes %in% tissue_sharing$Sex[tissue_sharing$Sex$n_DE > 9, "gene_name"]]

# tissue-shared genes distribution --
traits_cols <- c(traits_cols, "black")
names(traits_cols)[5] <- "gene"

d0 <- do.call(rbind.data.frame, lapply(traits, function(trait) tissue_sharing[[trait]][, c(1,2,3,8)]))
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
d0$col <- sapply(1:nrow(d0), function(row) color_variable(d0[row, "gene_name"], d0[row, "trait"]))
d0$trait <- factor(d0$trait, levels = c(traits, "gene"), order = T)
d0 <- rbind.data.frame(d0[d0$col != "gene",], d0[d0$col == "gene",])
# so the plot is not so heavy
df <- rbind.data.frame(do.call(rbind.data.frame, 
                               lapply(traits[-4], function(trait) 
  do.call(rbind.data.frame, lapply(1:4, function(i) d0[d0$trait == trait & d0$n_DE == i,][sample(nrow(d0[d0$trait == trait & d0$n_DE == i,]), 500),]))
  )),
  d0[d0$n_DE > 4 & d0$trait %in% traits[c(1:3)],],
  d0[d0$trait == "BMI",])

# plot
p1 <- ggplot(df,
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

pdf(paste0(plot_path, "Figure_2D.tissue_sharing.pdf"),
    width = 6, height = 6)
p1
dev.off()

# proportion of genes with a given tissue-sharing per trait --
total_DEGs <- sapply(traits, function(trait) nrow(tissue_sharing[[trait]]))
d <- rbind.data.frame(sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DE==1)/total_DEGs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DE>=2 & tissue_sharing[[trait]]$n_DE<=5)/total_DEGs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DE>=6 & tissue_sharing[[trait]]$n_DE<=9)/total_DEGs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DE>=10 )/total_DEGs[trait])
))
colnames(d) <- traits
rownames(d) <- c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")


df2 <- melt(d)
df2$type <- rep(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing"), 4)
df2$type <- factor(df2$type, levels = rev(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")), order = T) 
cols <- brewer.pal(4, "Greys")
names(cols) <- c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")
p2 <- ggplot(df2, aes(x = 1,
               y = value,
               fill = type)) +
  geom_bar(stat= "identity") + 
  coord_flip() +
  theme_bw() +
  facet_grid(~variable) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("") + ylab("DEGs (%)")


pdf(paste0(plot_path, "Figure_2D.tissue_sharing.bar_plot.pdf"),
    width = 12, height = 1.25)
p2
dev.off()


# Figure S2E ----
# The expression data  corresponds to gene-level TPM quantifications from the GTEx v8 main paper (GTEx Consortium. 2020), which are available on the GTEx portal
# From each tissue, we used a subset of samples with available metadata for the covariates included in the linear regression (STAR methods)
# gene examples of highly shared genes
tpm <- lapply(tissues, function(tissue) readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds")))
names(tpm) <- tissues

# Ancestry --
gene_name <- "GSTM3"
EA_value <- sapply(unlist(strsplit(tissue_sharing$Ancestry[tissue_sharing$Ancestry$gene_name==gene_name,"tissues_DE"], split = ";")), function(tissue)
  median(as.numeric(tpm[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"],
                                  mdata[[tissue]][mdata[[tissue]]$Ancestry=="EUR", "Sample"]]))
)
AA_value <- sapply(unlist(strsplit(tissue_sharing$Ancestry[tissue_sharing$Ancestry$gene_name==gene_name,"tissues_DE"], split = ";")), function(tissue)
  median(as.numeric(tpm[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"],
                                  mdata[[tissue]][mdata[[tissue]]$Ancestry=="AFR", "Sample"]]))
)

data <- cbind.data.frame("Ancestry" = c(rep("EA",length(EA_value)),
                                        rep("AA",length(AA_value))),
                         "value" = c(EA_value, AA_value),
                         "Tissue" = rep(unlist(strsplit(tissue_sharing$Ancestry[tissue_sharing$Ancestry$gene_name==gene_name,"tissues_DE"], split = ";")), 2))
data$Ancestry <- factor(data$Ancestry, levels = c("EA", "AA"), order = T)
data$Tissue <- factor(data$Tissue, levels = unlist(strsplit(tissue_sharing$Ancestry[tissue_sharing$Ancestry$gene_name==gene_name,"tissues_DE"], split = ";")),
                      order = T)
tissue_cols <- tissue_info[tissue_info$tissue_ID %in% as.character(levels(data$Tissue)), "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$gene <- rep(gene_name, nrow(data))
data$gene <- factor(data$gene)

p3 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = log2(value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("TPM (median log<sub>2</sub>)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  facet_grid(~gene) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") 

# Age --
gene_name <- "ZMAT3"
young_value <- sapply(unlist(strsplit(tissue_sharing$Age[tissue_sharing$Age$gene_name==gene_name,"tissues_DE"], split = ";")), function(tissue)
  median(as.numeric(tpm[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"],
                                  mdata[[tissue]][mdata[[tissue]]$Age < 45, "Sample"]]))
)
old_value <- sapply(unlist(strsplit(tissue_sharing$Age[tissue_sharing$Age$gene_name==gene_name,"tissues_DE"], split = ";")), function(tissue)
  median(as.numeric(tpm[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"],
                                  mdata[[tissue]][mdata[[tissue]]$Age >= 45, "Sample"]]))
)

data <- cbind.data.frame("Age" = c(rep("[20-45)",length(young_value)),
                                   rep("[45-70]",length(old_value))),
                         "value" = c(young_value, old_value),
                         "Tissue" = rep(unlist(strsplit(tissue_sharing$Age[tissue_sharing$Age$gene_name==gene_name,"tissues_DE"], split = ";")), 2))
data$Age <- factor(data$Age, levels = c("[20-45)", "[45-70]"), order = T)
data$Tissue <- factor(data$Tissue, levels = unlist(strsplit(tissue_sharing$Age[tissue_sharing$Age$gene_name==gene_name,"tissues_DE"], split = ";")),
                      order = T)
tissue_cols <- tissue_info[tissue_info$tissue_ID %in% as.character(levels(data$Tissue)), "colcodes"]
names(tissue_cols) <- as.character(levels(data$Tissue))
data$gene <- rep(gene_name, nrow(data))
data$gene <- factor(data$gene)

p4 <- ggplot(data = data,
             aes(x = Age,
                 y = log2(value),
                 col = Tissue),
) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_point(aes(col = Tissue),
             size = 0.8) +
  geom_line(aes(group=Tissue)) +
  xlab("") +
  ylab("TPM (median log<sub>2</sub>)") +
  scale_color_manual(values = tissue_cols) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Age"]),
        legend.position = "none") +
  facet_grid(~gene)

pdf(paste0(plot_path, "Figure_2E.tissue_sharing_examples.pdf"),
    width = 2, height = 8)
ggarrange(p3, p4, ncol = 1)
dev.off()
