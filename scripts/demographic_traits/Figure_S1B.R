#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)

# metadata ----
load("00.metadata.RData")

# differential expression analysis ----
dea <- readRDS("00.differential_expression_analysis.rds")

# differential splicing analysis ----
dsa <- readRDS("00.differential_splicing_analysis.rds")
for(tissue in tissues){
  dsa[[tissue]][["Age"]]$type <- sapply(rownames(dsa[[tissue]][["Age"]]), function(i) unlist(strsplit(unlist(strsplit(i,split = ";"))[[2]], split = ":"))[[1]])
}

# number of samples per tissue
n_samples <- sapply(tissues, function(tissue) nrow(mdata[[tissue]]))
names(n_samples) <- tissues
# sum(n_samples) 13,684 total number of samples analyzed

# number of unique donor
d <- do.call(rbind.data.frame, lapply(tissues, function(tissue) mdata[[tissue]][,c("Sample","Donor")]))
# length(unique(d$Donor)) 781

# Figure S1B and C ----
# ancestry 
t <- sapply(tissues, function(tissue) sapply(c("EUR","AFR"), function(i) sum(mdata[[tissue]]$Ancestry==i)))
colnames(t) <- sapply(tissues, function(tissue) tissue_info[tissue_info$tissue_ID==tissue,]$tissue_abbrv)
ancestry_data <- as.data.frame(t(t))

# sex 
t <- sapply(tissues, function(tissue) sapply(c("1","2"), function(i) sum(mdata[[tissue]]$Sex==i)))
colnames(t) <- sapply(tissues, function(tissue) tissue_info[tissue_info$tissue_ID==tissue,]$tissue_abbrv)
sex_data <- as.data.frame(t(t)[,1:2])
colnames(sex_data) <- c("Male","Female")

# number of genes expressed per tissue
exprs_genes <- sapply(tissues, function(tissue) nrow(dea[[tissue]][["Age"]]))
names(exprs_genes) <- tissue_info$tissue_abbrv

# biotype of genes expressed per tissue
exprs_genes_biotype <- sapply(tissues, function(tissue) table(gene_annotation[gene_annotation$ensembl.id %in% rownames(dea[[tissue]][["Age"]]), "biotype"]))
colnames(exprs_genes_biotype) <- tissue_info$tissue_abbrv

# alternatively spliced events (ASEs)
ASEs <- sapply(tissues, function(tissue) nrow(dsa[[tissue]][["Age"]]))
names(ASEs) <- tissue_info$tissue_abbrv

# type of ASEs
ASEs_biotype <- sapply(tissues, function(tissue) table(dsa[[tissue]][["Age"]]$type)[c("SE","MX","A5","A3","RI","AF","AL")])
colnames(ASEs_biotype) <- tissue_info$tissue_abbrv

# total number of ASEs analyzed
unique_ASEs <- unique(unlist(lapply(tissues, function(tissue) rownames(dsa[[tissue]][["Age"]]))))
# sum(table(sapply(unique_ASEs, function(i) unlist(strsplit(unlist(strsplit(i, split = ";"))[[2]], split = "\\:"))[[1]]))) 62,269

# heatmap 1 annotation --
df <- t(apply(sex_data, 1, function(x) x/sum(x)))
df[tissue_info[tissue_info$tissue_ID %in% sex_tissues, "tissue_abbrv"],] <- 0
row_ha_left <- HeatmapAnnotation(
                                  "Samples" = anno_barplot(n_samples,
                                                          gp = gpar(fill = tissue_info$colcodes,
                                                                    col = tissue_info$colcodes),
                                                          border=F),
                                 "Ancestry" = anno_barplot(t(apply(ancestry_data, 1, function(x) x/sum(x))),
                                                           gp = gpar(fill = c("#f5d899","#E69F00","light grey"),
                                                                     col = c("#f5d899","#E69F00","light grey")),
                                                           border=F),
                                 "Sex" = anno_barplot(df,
                                                      gp = gpar(fill = c("#66c4ab","#009E73"),
                                                                col = c("#66c4ab","#009E73")),
                                                      border=F),
                                 "Age" = anno_boxplot(lapply(tissues, function(tissue) mdata[[tissue]]$Age),
                                                      height = unit(4, "cm"),
                                                      gp = gpar(fill = tissue_info$colcodes),
                                                      border=F,axis=T,
                                                      size = unit(0.5, "mm")),
                                 "BMI" = anno_boxplot(lapply(tissues, function(tissue) mdata[[tissue]]$BMI),
                                                      height = unit(4, "cm"),
                                                      gp = gpar(fill = tissue_info$colcodes),
                                                      border=F,axis=T,
                                                      size = unit(0.5, "mm")),
                                 # "Expressed genes" = anno_barplot(exprs_genes,
                                 #                                  gp = gpar(fill = "gray",
                                 #                                            col = "gray"),
                                 #                                  border=F),
                                 # "Gene biotype" = anno_barplot(t(apply(exprs_genes_biotype[2:1,], 2, function(x) x/sum(x))),
                                 #                                  gp = gpar(fill = c("#2288dd", "#6dc781"),
                                 #                                            col = c("#2288dd", "#6dc781")),
                                 #                                  border=F),
                                 # "Alternatively spliced events" = anno_barplot(ASEs,
                                 #                                               gp = gpar(fill = "gray",
                                 #                                                              col = "gray"),
                                 #                                               border=F),
                                 # "Type of splicing event" = anno_barplot(t(apply(ASEs_biotype[1:7,], 2, function(x) x/sum(x))),
                                 #                                         gp = gpar(fill = splicing_events_cols,
                                 #                                                   col = splicing_events_cols),
                                 #                                         border = F),
                                 gap = unit(0.25,"cm"),
                                 #width = unit(12, "cm"),
                                 show_legend = T, 
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")


mock_matrix <- matrix(1, nrow = length(tissues), ncol = 1)

ht1 <- Heatmap(mock_matrix,
          na_col = "white",
          cluster_rows = F,
          cluster_columns = F,
          name = "DE signal",
          row_names_side = "left",
          column_names_side = "bottom",
          column_names_rot =  90,
          left_annotation = row_ha_left,
          row_names_gp = gpar(fontsize = 9),
          heatmap_legend_param = list(direction = "horizontal")
)

pdf(paste0(plot_path, "Figure_S1B.annotation_demographic_traits.pdf"),
    width = 3, height = 9)
draw(ht1)
dev.off()

# heatmap 2 annotation --
row_ha_left <- HeatmapAnnotation(
  "Expressed genes" = anno_barplot(exprs_genes,
                                   gp = gpar(fill = "gray",
                                             col = "gray"),
                                   border=F),
  "Gene biotype" = anno_barplot(t(apply(exprs_genes_biotype[2:1,], 2, function(x) x/sum(x))),
                                   gp = gpar(fill = c("#2288dd", "#6dc781"),
                                             col = c("#2288dd", "#6dc781")),
                                   border=F),
  "Alternatively spliced events" = anno_barplot(ASEs,
                                                gp = gpar(fill = "gray",
                                                               col = "gray"),
                                                border=F),
  "Type of splicing event" = anno_barplot(t(apply(ASEs_biotype[1:7,], 2, function(x) x/sum(x))),
                                          gp = gpar(fill = splicing_events_cols,
                                                    col = splicing_events_cols),
                                          border = F),
  gap = unit(0.25,"cm"),
  width = unit(8, "cm"),
  show_legend = T, 
  show_annotation_name = T,
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 10),
  which = "row")


mock_matrix <- matrix(1, nrow = length(tissues), ncol = 1)
ht2 <- Heatmap(mock_matrix,
               na_col = "white",
               cluster_rows = F,
               cluster_columns = F,
               name = "DE signal",
               row_names_side = "left",
               column_names_side = "bottom",
               column_names_rot =  90,
               left_annotation = row_ha_left,
               row_names_gp = gpar(fontsize = 9),
               heatmap_legend_param = list(direction = "horizontal")
)

pdf(paste0(plot_path, "Figure_S1B.genes_and_ASEs.pdf"),
    width = 4, height = 9)
draw(ht2)
dev.off()

pdf(paste0(plot_path, "Figure_S1B.legend.pdf"),
    width = 12, height = 6)
par(mfrow = c(4,1))
plot.new()
legend("center",
       c("EA", "AA"),
       col = c("grey", "#f5d899"),
       bty = 'n',
       pch = 15, ncol =2, cex = 3)
plot.new()
legend("center",
       c("Male", "Female"),
       col = c("#66c4ab","#009E73"),
       bty = 'n',
       pch = 15, ncol = 2, cex = 3)
plot.new()
legend("center",
       c("protein-coding", "lincRNA"),
       col = c("#2288dd", "#6dc781"),
       bty = 'n',
       pch = 15, ncol = 2, cex = 3)
plot.new()
legend("center",
       splicing_events,
       col = splicing_events_cols,
       bty = 'n',
       pch = 15, ncol = 7, cex = 3)
dev.off()
