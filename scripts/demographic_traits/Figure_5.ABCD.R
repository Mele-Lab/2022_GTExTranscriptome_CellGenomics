#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

# metadata ----
load("00.metadata.RData")

# ancestry DSEs in sGenes 
ancestry_DSE_sGenes <- readRDS("ancestry_DSE_sGenes.rds")

# Figure 5A ----
d <- cbind.data.frame("tissue" = tissues,
                      "cis-driven" = 100*sapply(tissues, function(tissue) nrow(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "cis-driven",])/nrow(ancestry_DSE_sGenes[[tissue]])))
d$tissue <- factor(d$tissue, levels = tissues, order = T)
d$dummy <- rep("ancestry-DSEs (sGenes)", nrow(d))

# spearman correlation between sample size and number/proportion of cis-driven DSEs
cor.test(sapply(tissues, function(tissue) nrow(mdata[[tissue]])),
  sapply(tissues, function(tissue) nrow(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "cis-driven",])),
  method = "spearman"
)
cor.test(sapply(tissues, function(tissue) nrow(mdata[[tissue]])),
         sapply(tissues, function(tissue) nrow(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "cis-driven",])/nrow(ancestry_DSE_sGenes[[tissue]])),
         method = "spearman"
)

# plot
p1 <- ggplot(d,
             aes(x = 1, y = `cis-driven`, col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Cis-driven events (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") +
  facet_grid(~dummy)

# Figure 5B ---
# update
Fst <- lapply(tissues, function(tissue) readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Fst/Tissues/",tissue,"/", tissue, ".sVariants.Fst.rds")))
names(Fst) <- tissues
for(tissue in tissues){
  ancestry_DSE_sGenes[[tissue]]$ensembl_id <- sapply(rownames(ancestry_DSE_sGenes[[tissue]]), function(gene) unlist(strsplit(gene, split = ";"))[[1]])
}
Fst_data <- cbind.data.frame("tissue" = tissues, 
                             "cis-driven" = sapply(tissues, function(tissue) 
                               median(Fst[[tissue]][unique(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "cis-driven", "ensembl_id"]), "avg.Fst"], na.rm = T)),
                             "not cis-driven" = sapply(tissues, function(tissue) 
                               median(Fst[[tissue]][unique(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "not cis-driven", "ensembl_id"]), "avg.Fst"], na.rm = T))
)
Fst_data <- melt(Fst_data)
Fst_data$variable <- factor(Fst_data$variable, levels = c("cis-driven", "not cis-driven"), order = T)

p2 <- ggplot(Fst_data,
             aes(x = variable, 
                 y = value, 
                 col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  #geom_line(aes(group=tissue)) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Fst") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=brewer.pal(11, "PRGn")[c(3,9)]),
        legend.position = "none") +
  facet_grid(~variable, scale = "free")

wilcox.test(Fst_data[Fst_data$variable == "cis-driven", "value"],
            Fst_data[Fst_data$variable == "not cis-driven", "value"],
            paired = T)

# Figure 5C ----
tissue_sharing <- readRDS("events_DS.tissue_sharing.rds")
tissue_sharing_data <- cbind.data.frame("tissue" = tissues, 
                                        "cis-driven" = sapply(tissues, function(tissue) 
                                          median(tissue_sharing$Ancestry[tissue_sharing$Ancestry$event_id %in% rownames(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "cis-driven",]), "n_DS"])),
                                        "not cis-driven" = sapply(tissues, function(tissue) 
                                          median(tissue_sharing$Ancestry[tissue_sharing$Ancestry$event_id %in% rownames(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type == "not cis-driven",]), "n_DS"]))
)
tissue_sharing_data <- melt(tissue_sharing_data)
tissue_sharing_data$variable <- factor(tissue_sharing_data$variable, levels = c("cis-driven", "not cis-driven"), order = T)

p3 <- ggplot(tissue_sharing_data,
             aes(x = variable, 
                 y = value, 
                 col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  #geom_line(aes(group=tissue)) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Number of tissues") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=brewer.pal(11, "PRGn")[c(3,9)]),
        legend.position = "none") +
  facet_grid(~variable, scale = "free")

wilcox.test(tissue_sharing_data[tissue_sharing_data$variable == "cis-driven", "value"],
            tissue_sharing_data[tissue_sharing_data$variable == "not cis-driven", "value"],
            paired = T)

# Figure 5D ----
data <- cbind.data.frame("tissue" = rep(tissues, 3),
                         "value" = c(sapply(tissues, function(tissue) summary(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type=="cis-driven", "isQTLs_R2"])[3]),
                                     sapply(tissues, function(tissue) summary(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type=="not cis-driven", "isQTLs_R2"])[3]),
                                     sapply(tissues, function(tissue) summary(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type=="not cis-driven", "Ancestry_R2"])[3])),
                         "feature" = c(rep("sQTLs", length(tissues)*2), rep("Ancestry", length(tissues))),
                         "type" = c(rep("cis-driven", length(tissues)), rep("cis-independent", length(tissues)*2)))
data$tissue <- factor(data$tissue, levels = tissues, order = T)
data$feature <- factor(data$feature, levels = c("sQTLs", "Ancestry"), order = T)
data$type <- factor(data$type, levels = c("cis-driven", "cis-independent"), order = T)
data$class <- paste0(data$feature, "_", data$type)
data$class <- factor(data$class, levels = c("sQTLs_cis-driven", "sQTLs_cis-independent", "Ancestry_cis-independent"), order = T)
unique(data$class)

p4 <- ggplot(data,
             aes(x= feature, y =  value, col = feature, fill = feature)) +
  geom_violin(col = "black", fill = NA) +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  scale_color_manual(values = tissue_info$colcodes) +
  facet_grid(~class, scale = "free") +
  ylab("Splicing variation explained (%)") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") 
p4
pdf(paste0(plot_path, "Figure_5ABCD.cis_driven.NatureGenetics.pdf"),
    width = 9, height = 3)
ggarrange(p1, p2, p3, p4, widths = c(1,2,2,3), ncol = 4)
dev.off()

# data <- cbind.data.frame("tissue" = rep(tissues, 2),
#                          "value" = c(sapply(tissues, function(tissue) summary(ancestry_DSE_sGenes[[tissue]][, "isQTLs_R2"])[3]),
#                                      sapply(tissues, function(tissue) summary(ancestry_DSE_sGenes[[tissue]][ancestry_DSE_sGenes[[tissue]]$type=="not cis-driven", "Ancestry_R2"])[3])),
#                          "feature" = c(rep("sQTLs", length(tissues)), rep("Ancestry", length(tissues))))
# data$tissue <- factor(data$tissue, levels = tissues, order = T)
# data$feature <- factor(data$feature, levels = c("sQTLs", "Ancestry"), order = T)
# 
# p4 <- ggplot(data,
#              aes(x= feature, y =  value, col = feature, fill = feature)) +
#   geom_violin(col = "black", fill = NA) +
#   geom_boxplot(col = "black",
#                fill = "white",
#                outlier.shape = NA,
#                notch = T,
#                width = 0.25) + 
#   geom_jitter(aes(col = tissue), 
#               size = 2) +
#   scale_color_manual(values = tissue_info$colcodes) +
#   facet_grid(~feature, scale = "free") +
#   ylab("Gene expression variation explained (%)") +
#   xlab("") +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(hjust = 0.5,
#                                   size = 15),
#         strip.background = element_rect(fill=traits_cols["Ancestry"]),
#         legend.position = "none") 

wilcox.test(data[data$feature == "sQTLs", "value"],
            data[data$feature == "Ancestry", "value"],
            paired = T)

pdf(paste0(plot_path, "Figure_4IJKL.cis_driven.pdf"),
    width = 9, height = 3)
ggarrange(p1, p2, p3, p4, widths = c(1,2,2,2), ncol = 4)
dev.off()

