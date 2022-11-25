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
library(dplyr)

# metadata ----
load("00.metadata.RData")

# differential splicing analysis ----
dsa <- readRDS("differential_splicing_analysis.rds")

# ancestry DSEs sGenes ----
ancestry_DSEs_sGenes <- readRDS("ancestry_DSE_sGenes.rds")

# sGene data ----
inpath <- "~/GTEx_v8_data/cisQTLs/"
sgene_data <- lapply(tissues, function(tissue)
  as.data.frame(data.table::fread(paste0(inpath,"sGenes/GTEx_Analysis_v8_sQTL/",tissue_info[tissue_info$tissue_ID==tissue,"tissue_id"],".v8.sgenes.txt.gz"))))
names(sgene_data) <- tissues

# select sGenes (q.val < 0.05) --
sGenes <- lapply(tissues, function(tissue)
  sgene_data[[tissue]][sgene_data[[tissue]]$qval <= 0.05,"gene_id"]
)
names(sGenes) <- tissues

# Figure S5A ----
ancestry_DSEs <- lapply(tissues, function(tissue) rownames(dsa[[tissue]][["Ancestry"]][dsa[[tissue]][["Ancestry"]]$adj.P.Val < 0.05,]) )
names(ancestry_DSEs) <- tissues
for(tissue in tissues){
  names(ancestry_DSEs[[tissue]]) <- sapply(ancestry_DSEs[[tissue]], function(i) unlist(strsplit(i, split = ";"))[[1]])
}

# number of ancestry DSEs
data <- melt(cbind.data.frame(
  "tissue" = tissue_info$tissue_abbrv,
  "not sGene" = sapply(tissues, function(tissue) sum(!names(ancestry_DSEs[[tissue]]) %in% sGenes[[tissue]])),
  "not cis-driven" = sapply(tissues, function(tissue) sum(ancestry_DSEs_sGenes[[tissue]]$type == "not cis-driven")),
  "cis-driven" = sapply(tissues, function(tissue) sum(ancestry_DSEs_sGenes[[tissue]]$type == "cis-driven")),
  "not classified" = sapply(tissues, function(tissue) sum(!ancestry_DSEs[[tissue]][names(ancestry_DSEs[[tissue]]) %in% sGenes[[tissue]]] %in% rownames(ancestry_DSEs_sGenes[[tissue]])) )))
data$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$variable <- factor(data$variable, levels = rev(c("not sGene", "not cis-driven", "cis-driven", "not classified")), order = T)

cols <- c(brewer.pal(11,"PRGn")[c(2,3,9)], "light grey")
names(cols) <- c("not sGene", "not cis-driven", "cis-driven", "not classified")

p1 <- ggplot(data, mapping = aes(x = value, 
                                 fill = variable, #actor(var_2, levels = rev(x_levs)), 
                                 y = tissue, #factor(var_1, levels = rev(y_levs)),
                                 label = value)) +
  geom_bar( stat = "identity") + 
  #geom_text(size = 2,
  #          angle = 45,
  #          position = position_stack(vjust = 0.5)) + # Add number labelsfacet_grid(~trait) +
  ylab("") + xlab("Number of DSEs") +
  labs(fill = "") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  theme(legend.position = "none")

# proportion of ancestry DEGs 
data2 <- data <- melt(cbind.data.frame(
  "tissue" = tissue_info$tissue_abbrv,
  "not sGene" = sapply(tissues, function(tissue) sum(!names(ancestry_DSEs[[tissue]]) %in% sGenes[[tissue]]))/sapply(tissues, function(tissue) length(ancestry_DSEs[[tissue]])),
  "not cis-driven" = sapply(tissues, function(tissue) sum(ancestry_DSEs_sGenes[[tissue]]$type == "not cis-driven"))/sapply(tissues, function(tissue) length(ancestry_DSEs[[tissue]])),
  "cis-driven" = sapply(tissues, function(tissue) sum(ancestry_DSEs_sGenes[[tissue]]$type == "cis-driven"))/sapply(tissues, function(tissue) length(ancestry_DSEs[[tissue]])),
  "not classified" = sapply(tissues, function(tissue) sum(!ancestry_DSEs[[tissue]][names(ancestry_DSEs[[tissue]]) %in% sGenes[[tissue]]] %in% rownames(ancestry_DSEs_sGenes[[tissue]])) )/sapply(tissues, function(tissue) length(ancestry_DSEs[[tissue]]))
))

data2$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data2$variable <- factor(data$variable, levels = rev(c("not sGene", "not cis-driven", "cis-driven", "not classified")), order = T)
#data2$value <- round(100*data2$value)

p2 <- ggplot(data2, mapping = aes(x = 100*value, 
                                  fill = variable, #actor(var_2, levels = rev(x_levs)), 
                                  y = tissue, #factor(var_1, levels = rev(y_levs)),
                                  label = value)) +
  geom_bar( stat = "identity") + 
  #geom_text(size = 2
  #          , position = position_stack(vjust = 0.5)) + # Add number labelsfacet_grid(~trait) +
  ylab("") + xlab("DSEs (%)") +
  labs(fill = "") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  theme(legend.position = "none")

# plot
pdf(paste0(plot_path, "Figure_S5A.pdf"),
    width = 12, height = 9)
ggarrange(p1, p2, widths = c(1,0.5))
dev.off()

# Figure S5 B-C ----
# S5C Fst
Fst <- readRDS("Fst_sGenes.rds")

# S5D tissue-sharing
tissue_sharing <- readRDS("events_DS.tissue_sharing.rds")[["Ancestry"]]

for(tissue in tissues){
  ancestry_DSEs_sGenes[[tissue]]$ensembl_id <- sapply(rownames(ancestry_DSEs_sGenes[[tissue]]), function(e) unlist(strsplit(e, split = ";"))[[1]])
  ancestry_DSEs_sGenes[[tissue]]$tissue <- rep(tissue_info[tissue_info$tissue_ID == tissue, "tissue_abbrv"])
}
parse_data <- function(tissue){
  d1 <- merge(ancestry_DSEs_sGenes[[tissue]], Fst[[tissue]], by = "ensembl_id")
  d2 <- merge(d1, tissue_sharing[, c("ensembl_id", "n_DS")])
  d3 <- d2[,c(1,2,5,6,7)]
  colnames(d3)[4] <- "Fst"
  return(d3)
}

data <- do.call(rbind.data.frame, lapply(tissues, function(tissue) parse_data(tissue)))
data$tissue <- factor(data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
data$type <- factor(data$type, levels = c("cis-driven", "not cis-driven"), order = T)

# number of tissues with significantly different distributions between cis-driven and not cis-driven
sum(p.adjust(sapply(tissue_info$tissue_abbrv, function(tissue)
  wilcox.test(data[data$tissue == tissue & data$type == "cis-driven", "Fst"],
              data[data$tissue == tissue & data$type == "not cis-driven", "Fst"],
              alternative  = "greater")$p.value
)) < 0.05)
sum(p.adjust(sapply(tissue_info$tissue_abbrv, function(tissue)
  wilcox.test(data[data$tissue == tissue & data$type == "cis-driven", "n_DS"],
              data[data$tissue == tissue & data$type == "not cis-driven", "n_DS"],
              alternative  = "greater")$p.value
)) < 0.05)


# plot 
data <-
  data %>%
  group_by(tissue) %>%
  mutate(outlier = Fst > median(Fst) + IQR(Fst) * 1.5) %>%
  ungroup

p1 <- ggplot(data, aes(x = tissue,
                       y = Fst,
                       col = type,
                       fill = type)) +
  geom_boxplot(outlier.size = 0.1) +
  #geom_jitter(data = function(x) dplyr::filter_(x, ~ outlier), position = 'jitter', size = 0.1) + # Outliers
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values =  c(brewer.pal(11,"PRGn")[9], brewer.pal(11,"PRGn")[2])) +
  scale_fill_manual(values = c(alpha(brewer.pal(11,"PRGn")[9], 0.5), alpha(brewer.pal(11,"PRGn")[2], 0.5))) +
  ylab("Fst") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))



data2 <-
  data %>%
  group_by(tissue) %>%
  mutate(outlier = n_DS > median(n_DS) + IQR(n_DS) * 1.5) %>%
  ungroup

p2 <- ggplot(data2, aes(x = tissue,
                        y = n_DS,
                        col = type,
                        fill = type)) +
  geom_boxplot(outlier.shape = NA) +  # no outliers
  geom_point(data = function(x) dplyr::filter_(x, ~ outlier), 
             position = 'jitter', 
             size = 0.1) + # Outliers 
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values =  c(brewer.pal(11,"PRGn")[9], brewer.pal(11,"PRGn")[2])) +
  scale_fill_manual(values = c(alpha(brewer.pal(11,"PRGn")[9], 0.5), alpha(brewer.pal(11,"PRGn")[2], 0.5))) +
  ylab("Number of tissues") + xlab("") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))

pdf(paste0(plot_path, "Figure_S5_BC.pdf"),
    width = 9, height = 10)
ggarrange(p1, p2)
dev.off()


# Figure S5D ----
for(tissue in tissues){
  ancestry_DSEs_sGenes[[tissue]]$tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(ancestry_DSEs_sGenes[[tissue]]))
}
R2_data <- melt(do.call(rbind.data.frame, ancestry_DSEs_sGenes))
R2_data$tissue <- factor(R2_data$tissue, levels = rev(tissue_info$tissue_abbrv), order = T)
R2_data$variable <- factor(R2_data$variable, levels = c("isQTLs_R2", "Ancestry_R2"))

# p1 <- ggplot(R2_data,
#              aes(x = tissue,
#                  y = value,
#                  fill = variable,
#                  col = variable)) +
#   geom_boxplot(outlier.size = 0.1) +
#   facet_grid(~variable, scale = "free_x") +
#   coord_flip() +
#   theme_bw() +
#   scale_colour_manual(values = c("grey", "#E69F00") ) +
#   scale_fill_manual(values = c(alpha("grey", 0.5), alpha("#E69F00", 0.5))) +
#   ylab("Splicing variation explained (%)") + xlab("") +
#   #stat_summary(fun.data = get_box_stats, geom = "text",
#   #             hjust = 0.5, vjust = 0.9, size = 2) +
#   theme(legend.position = 'none',
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 10)) +
#   scale_y_continuous(limits = c(0,30))
# 
# pdf(paste0(plot_path, "Figure_S5D.pdf"),
#     width = 6, height = 10)
# ggarrange(p1)
# dev.off()

# number of tissues with significantly different amount of splicing variation explained by sQTLs and ancestry
sum(p.adjust(sapply(tissue_info$tissue_abbrv, function(tissue)
  wilcox.test(R2_data[R2_data$tissue == tissue & R2_data$variable == "isQTLs_R2", "value"],
              R2_data[R2_data$tissue == tissue & R2_data$variable == "Ancestry_R2", "value"],
              alternative  = "greater")$p.value
)) < 0.05)

p1 <- ggplot(R2_data[R2_data$variable == "isQTLs_R2",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~type, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = c("grey", "#E69F00") ) +
  scale_fill_manual(values = c(alpha("grey", 0.5), alpha("#E69F00", 0.5))) +
  ylab("Splicing variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10))
p2 <- ggplot(R2_data[R2_data$variable == "Ancestry_R2" &
                       R2_data$type == "not cis-driven",],
             aes(x = tissue,
                 y = value,
                 fill = variable,
                 col = variable)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_grid(~type, scale = "free_x") +
  coord_flip() +
  theme_bw() +
  scale_colour_manual(values = c("#E69F00", "#E69F00") ) +
  scale_fill_manual(values = c(alpha("#E69F00", 0.5), alpha("#E69F00", 0.5))) +
  ylab("Splicing variation explained (%)") + xlab("") +
  #stat_summary(fun.data = get_box_stats, geom = "text",
  #             hjust = 0.5, vjust = 0.9, size = 2) +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  scale_y_continuous(limits = c(0,30))

pdf(paste0(plot_path, "Figure_S5D.NatureGenetics.pdf"),
    width = 10, height = 9)
ggarrange(p1, p2, widths = c(2,1.25))
dev.off()

round(mean(sapply(tissue_info$tissue_abbrv, function(tissue) median(R2_data[R2_data$tissue == tissue & 
                                                                              R2_data$type == "cis-driven" &
                                                                              R2_data$variable == "isQTLs_R2", "value"], na.rm = T))))
round(mean(sapply(tissue_info$tissue_abbrv, function(tissue) median(R2_data[R2_data$tissue == tissue & 
                                                                              R2_data$type == "not cis-driven" &
                                                                              R2_data$variable == "isQTLs_R2", "value"], na.rm = T))))
round(mean(sapply(tissue_info$tissue_abbrv, function(tissue) median(R2_data[R2_data$tissue == tissue & 
                                                                              R2_data$type == "not cis-driven" &
                                                                              R2_data$variable == "Ancestry_R2", "value"], na.rm = T))))
