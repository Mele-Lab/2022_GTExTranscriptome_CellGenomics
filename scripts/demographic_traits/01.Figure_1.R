#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(ggpubr)

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

# total unique number of DEGs per tissue
get_DEGs <- function(tissue, trait){
  if(tissue %in% sex_tissues & trait == "Sex"){
    NA
  }else{
    rownames(dea[[tissue]][[trait]][dea[[tissue]][[trait]]$adj.P.Val < 0.05,])
  }
}
DEGs <- lapply(traits, function(trait) lapply(tissues, function(tissue) get_DEGs(tissue, trait)))
names(DEGs) <- traits
for(trait in traits){names(DEGs[[trait]]) <- tissues}
total_DEGs <- sapply(traits, function(trait) length(unique(unlist(DEGs[[trait]]))[!is.na(unique(unlist(DEGs[[trait]])))]))

# Figure 1A ----
row_ha_left <- HeatmapAnnotation("Samples" = anno_barplot(sapply(tissues, function(tissue) nrow(mdata[[tissue]])),
                                                          gp = gpar(fill = tissue_info$colcodes,
                                                                    col = tissue_info$colcodes),
                                                          border=F),
                                 gap = unit(0.25,"cm"),
                                 show_legend = T, 
                                 show_annotation_name = T,
                                 annotation_name_rot = 90,
                                 annotation_name_gp = gpar(fontsize = 10),
                                 which = "row")
column_ha_top <- HeatmapAnnotation("total unique number of DEGs" = anno_barplot(total_DEGs,
                                                border = F,
                                                gp = gpar(fill = traits_cols,
                                                          col = traits_cols)),
                                   show_annotation_name = T,
                                   annotation_name_side = "left",
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = 6.5),
                                   height = unit(3, "cm"))

ht_DEGs <- Heatmap(apply(number_of_DEGs, 2,function(x) x/max(x,na.rm=T)),
        col = brewer.pal(9,"BuPu"),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "DE signal",
        row_names_side = "left",
        row_labels = tissue_info$tissue_abbrv,
        column_names_side = "bottom",
        column_names_rot =  90,
        top_annotation = column_ha_top,
        left_annotation = row_ha_left,
        row_names_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 9),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(prettyNum(number_of_DEGs[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
        heatmap_legend_param = list(direction = "horizontal")
)


pdf(paste0(plot_path, "Figure_1A.number_of_DEGs.pdf"),
    width = 3, height = 9)
draw(ht_DEGs,
     heatmap_legend_side = "bottom")
dev.off()

# Figure 1B ----
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

# proportion of total tissue expression variation explained by each trait --
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

# bar plot
pdf(paste0(plot_path, "Figure_1B.tissue_expression_variation_explained.pdf"),
    width = 3, height = 9)
barplot(t(tissue_expression_variation_explained[length(tissues):1,]),
        horiz = T,
        border = NA,
        col = traits_cols,
        xlab = "Tissue expression variation explained (%)",
        las = 2, 
        cex.names = 0.9,
        xaxt = 'n',
        yaxt = 'n')
axis(1, at = axTicks(1))
dev.off()

# number of tissues in which each trait is the main driver of tissue expression variation --
driver_traits <- as.vector(table(apply(tissue_expression_variation_explained, 1, function(x) traits[which.max(x)]))[traits])
names(driver_traits) <- traits

# top annotation bar plot
column_ha_top <- HeatmapAnnotation("number of tissues" = anno_barplot(driver_traits,
                                                                                border = F,
                                                                                gp = gpar(fill = traits_cols,
                                                                                          col = traits_cols)),
                                   show_annotation_name = T,
                                   annotation_name_side = "left",
                                   annotation_name_rot = 90,
                                   annotation_name_gp = gpar(fontsize = 6.5),
                                   height = unit(3, "cm"))

ht_DEGs_2 <- Heatmap(apply(number_of_DEGs, 2,function(x) x/max(x,na.rm=T)),
                   col = brewer.pal(9,"BuPu"),
                   na_col = "white",
                   cluster_rows = F,
                   cluster_columns = F,
                   name = "DE signal",
                   row_names_side = "left",
                   row_labels = tissue_info$tissue_abbrv,
                   column_names_side = "bottom",
                   column_names_rot =  90,
                   top_annotation = column_ha_top,
                   left_annotation = row_ha_left,
                   row_names_gp = gpar(fontsize = 9),
                   column_names_gp = gpar(fontsize = 9),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(prettyNum(number_of_DEGs[i, j], big.mark = ","), x, y, gp = gpar(fontsize = 10))},
                   heatmap_legend_param = list(direction = "horizontal")
)
pdf(paste0(plot_path, "Figure_1B.tissue_expression_variation_explained.bar_plot_top.pdf"),
    width = 3, height = 9)
draw(ht_DEGs_2,
     heatmap_legend_side = "bottom")
dev.off()

# Figure 1C ----
DEGs_per_trait_and_tissue <- unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, sum(dea[[tissue]][[trait]]$adj.P.Val < 0.05)))))
DEGs_per_tissue <- sapply(tissues, function(tissue) length(unique(unlist(lapply(traits, function(trait) DEGs[[trait]][[tissue]])))[!is.na(unique(unlist(lapply(traits, function(trait) DEGs[[trait]][[tissue]]))))]))
proporton_of_DEGs_per_trait <- 100*DEGs_per_trait_and_tissue/rep(DEGs_per_tissue,4)

# average gene expression variation explained per trait and tissue --
data <- cbind.data.frame("Tissue" = rep(tissues, 4),
                  "Trait" = unlist(lapply(traits, function(trait) rep(trait, length(tissues)))),
                  "DEGs" = unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, sum(dea[[tissue]][[trait]]$adj.P.Val < 0.05))))),
                  "proportion_of_DEGs" <- proporton_of_DEGs_per_trait,
                  "R2" = unlist(lapply(traits, function(trait) sapply(tissues, function(tissue)  ifelse(tissue %in% sex_tissues & trait == "Sex", NA, mean(dea[[tissue]][[trait]]$R2, na.rm = T)))))
                 )
data$Tissue <- factor(data$Tissue, levels = rev(tissues), order = T)
data$Trait <- factor(data$Trait, levels = traits, order = T)

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained.pdf"),
    width = 3, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean gene expression variation explained (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

pdf(paste0(plot_path, "Figure_1C.gene_expression_variation_explained.legend.pdf"),
    width = 5, height = 9)
ggplot(data,
       aes(x = Tissue,
           y = R2,
           col = Trait,
           size = proportion_of_DEGs)) +
  geom_point(alpha = 0.5) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = traits_cols) +
  ylab("Mean gene expression variation explained (%)") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
dev.off()

# Figure 1D ----
# The expression data  corresponds to gene-level TPM quantifications from the GTEx v8 main paper (GTEx Consortium. 2020), which are available on the GTEx portal
# From each tissue, we used a subset of samples with available metadata for the covariates included in the linear regression (STAR methods)
# function to plot N in box plots ---
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}

get_gene_data <- function(tissue, tpm, gene_name){
  ensembl_id <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
  df <- cbind.data.frame(tpm[ensembl_id,])
  colnames(df) <- "TPM"
  df$Ancestry <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$Age_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Age"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  df$BMI_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"BMI"])
  df$BMI <- sapply(df$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
  df$BMI <- factor(df$BMI,
                   levels = c("Normal", "Overweight", "Obese"),
                   order = T)
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$ensembl_id <- rep(ensembl_id, nrow(df))
  return(df)
}

get_gene_plot <- function(trait, tissue, gene_name){
  print(paste0(trait, " - ", tissue, " - ", gene_name))
  # Gene TPM --
  tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
  # tpm: matrix with gene TPM expression data for spleen samples used in this study. T
  
  # get gene data --
  data <- get_gene_data(tissue, tpm, gene_name)
  cols <- rep(traits_cols[trait], length(levels(data[,trait])))
  names(cols) <- as.character(levels(data[,trait]))
  
  # plot 1 --
  p1 <- ggplot(data = data,
               aes(x = eval(parse(text=trait)),
                   y = log2(1+TPM),
                   fill =  eval(parse(text=trait))),
  ) +
    geom_violin(col = "black") +
    geom_boxplot(col = "black",
                 fill = "white",
                 outlier.shape = NA,
                 notch = T,
                 width = 0.25) +
    geom_jitter(col = "black", 
                alpha = 0.1,
                size = 0.8) +
    xlab("") +
    ylab("TPM (log<sub>2</sub>)") +
    scale_fill_manual(values = cols) +
    stat_summary(fun.data = get_box_stats, geom = "text",
                 hjust = 0.5, vjust = 0.9, size = 3) +
    labs(title=paste0(gene_name, " (", tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], ")")) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.title.y = ggtext::element_markdown(size=14),
          plot.title = element_text(hjust = 0.5,
                                    size = 15),
          strip.background = element_rect(fill="#B3B3B3"),
          legend.position = "none") 
  
  # plot 2 --
  d <- data.frame("variable" = trait,
                  "value"  = dea[[tissue]][[trait]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], "R2"])
  variable_col <- traits_cols[trait]
  names(variable_col) <- trait
  
  p2 <- ggplot(d, aes(x=variable, y=value)) + 
    geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
    scale_fill_manual(values = variable_col) +
    coord_flip() +
    xlab("") + ylab("Gene expression variation explained (%)") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5,
                                    size = 15),
          strip.background = element_rect(fill="#B3B3B3"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y =  element_blank()) 
  p <- ggarrange(p1, p2, nrow = 2, heights =  c(9,2))
  return(p)
}

gene_examples <- list(c("Ancestry", "Spleen", "ACKR1"),
                      c("Sex", "ColonSigmoid", "LINC01597"),
                      c("Age", "ArteryAorta", "ROBO2"),
                      c("BMI", "AdiposeSubcutaneous", "SLC27A2"))

for(i in 1:length(gene_examples)){
  p <- get_gene_plot(gene_examples[[i]][1], gene_examples[[i]][2], gene_examples[[i]][3])  
  ggexport(p,filename = paste0(plot_path, "Figure_1D.", gene_examples[[i]][1], "_", gene_examples[[i]][2], "_", gene_examples[[i]][3], ".pdf"),
               width = 3, height = 4)
}

pdf(paste0(plot_path, "Demographic_traits.legend.pdf"), width = 8, height = 3)
plot.new()
legend("center",
       traits,
       col = traits_cols,
       pch = 15,
       bty = 'n', ncol = 4)
dev.off()
