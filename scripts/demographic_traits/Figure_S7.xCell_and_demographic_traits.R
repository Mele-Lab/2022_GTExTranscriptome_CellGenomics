#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)
library(circlize)
library(ggplot2)
library(ggarrange)

# metadata ----
load("00.metadata.RData")
for(tissue in c("Uterus","Vagina","Ovary")){
  mdata[[tissue]]$Sex <- rep("2", nrow(mdata[[tissue]]))
}
for(tissue in c("Testis","Prostate")){
  mdata[[tissue]]$Sex <- rep("1", nrow(mdata[[tissue]]))
}
for(tissue in tissues){
  print(tissue)
  mdata[[tissue]]$Tissue <- rep(tissue, nrow(mdata[[tissue]]))
  mdata[[tissue]] <- mdata[[tissue]][,c("Donor","Sample",
                                              "Age","Ancestry","Sex","BMI",
                                              "IschemicTime",  "RIN",
                                              "Tissue")]
}
metadata <- do.call(rbind.data.frame,
                 mdata)
rownames(mdata) <- NULL

# xCell enrichment scores for 7 cell types (adipocytes, epithelial cells, hepatocytes, keratinocytes, myocytes, neurons, neutrophils) across all tissue samples ----
# file: https://gtexportal.org/home/datasets/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz
file <- "~/Downloads/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz"
xCell <- read.table(gzfile(file), row.names = 1, header = T)
colnames(xCell) <- gsub("\\.", "-", colnames(xCell))

# benchmarked cell types 
cell_types <- rownames(xCell)

# select only tissues with benchmarked cell types 
xCell.tissues <- lapply(tissues, function(tissue) xCell[, mdata[[tissue]]$Sample])
names(xCell.tissues) <- tissues
tissue.cellTypes <- lapply(tissues, function(tissue) cell_types[apply(xCell.tissues[[tissue]][cell_types,], 1, function(x) median(x)) > 0.1])
names(tissue.cellTypes) <- tissues
tissue.cellTypes <- tissue.cellTypes[sapply(tissues, function(tissue) length(tissue.cellTypes[[tissue]]) > 0)]
tissues <- names(tissue.cellTypes)

# inverse normal transformed xCell scores 
i.xCell <-  t(apply(xCell,1, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))) )
i.xCell.tissues <- lapply(tissues, function(tissue) i.xCell[, mdata[[tissue]]$Sample])
names(i.xCell.tissues) <- tissues

# prepare data for analysis
parse.data <- function(tissue){
  print(tissue)
  d <- as.data.frame(i.xCell.tissues[[tissue]])
  d$CellType <- rownames(d)
  m <- melt(d)
  colnames(m) <- c("CellType","Sample","ixCell")
  m <- m[m$CellType %in% tissue.cellTypes[[tissue]],]
  m$Sample <- as.character(m$Sample)
  df <- as.data.frame(t(sapply(mdata[[tissue]]$Sample, function(sample) 
    sapply(tissue.cellTypes[[tissue]], function(i)
      m[m$Sample==sample & m$CellType==i,"ixCell"]
    )
  )
  ))
  
  if(length(tissue.cellTypes[[tissue]])==1){
    df <- cbind.data.frame(unlist(df), gsub(paste0(".",tissue.cellTypes[[tissue]]) , "", names(df)))
    colnames(df) <- c(tissue.cellTypes[[tissue]], "Sample")
    colnames(df) <- gsub(" ","_", colnames(df))
  }else{
    df$Sample <- rownames(df)
    colnames(df) <- gsub(" ","_", colnames(df))
  }
  
  return(df)
}

ixCell <- lapply(tissues, function(tissue) parse.data(tissue))
names(ixCell) <- tissues

# merge with tissue metadata 
data <- lapply(tissues, function(tissue) 
  merge(mdata[[tissue]],
        ixCell[[tissue]],
        by = "Sample")
)
names(data) <- tissues

# changes in cell type composition with demographic traits 
compositional_changes_fun <- function(tissue){
  print(tissue)
  
  data <- data[[tissue]]
  
  # cell type(s) in tissue
  covariates <-  c("IschemicTime",  "RIN")
  cell_types <- tissue.cellTypes[[tissue]]
  
  # traits
  if(!tissue %in% sex_tissues){
    traits <- c("Age","Ancestry", "Sex","BMI")
    traits.names <- c("Age","AncestryAFR", "Sex2","BMI")
  }else{
    traits <- c("Age","Ancestry", "BMI")
    traits.names <- c("Age","AncestryAFR", "BMI")
  }
  
  # cell_type ~ ischemic time + RIN + Age + Ancestry + Sex + BMI
  lm.obj <- lapply(cell_types, function(cell_type) 
    lm(as.formula(paste(cell_type," ~  ", paste(paste(c(covariates,traits), collapse = " + "), collapse = " "))),
       data)
  )
  names(lm.obj) <- cell_types
  
  # P.value matrix: rows=traits; columns=cell_types
  m <- sapply(cell_types, function(cell_type) summary(lm.obj[[cell_type]])$coefficients[,4])
  m <- m[-1,]
  if(length(cell_types) > 1){
    m <- m[traits.names,]
    M <- t(m)#t(apply(m, 1, function(x) p.adjust(x, method = "BH")))  
  }else{
    M <- m
    M <- m[traits.names]
  }
  
  return(list("p.value" = M))
}

# results 
results <- lapply(tissues, function(tissue) compositional_changes_fun(tissue))
names(results) <- tissues


# parse and organize data
get.tissue.table <- function(tissue){
  print(tissue)
  if(length(tissue.cellTypes[[tissue]])==1){
    df <- cbind.data.frame(rep(tissue, length(results[[tissue]]$p.value)),
                           gsub("2", "", gsub("AFR","", names(results[[tissue]]$p.value))),
                           rep(gsub(" ", "_", tissue.cellTypes[[tissue]]),length(results[[tissue]]$p.value)),
                           results[[tissue]]$p.value)
    colnames(df) <- c("Tissue","Trait","CellType","p.value")
  }else{
    df <- cbind.data.frame(melt(t(results[[tissue]]$p.value)))
    df$Tissue <- rep(tissue, nrow(df))
    df <- df[,c(4,1,2,3)]
    colnames(df) <- c("Tissue","Trait","CellType","p.value")
    df$Trait <- gsub("AncestryAFR", "Ancestry", df$Trait)
    df$Trait <- gsub("Sex2", "Sex", df$Trait)
  }
  return(df)
}

df <- do.call(rbind.data.frame, 
              lapply(tissues, function(tissue) get.tissue.table(tissue)))
df.traits <- lapply(traits, function(trait)
  df[df$Trait==trait,]
)
names(df.traits) <- traits
for(trait in traits){
  df.traits[[trait]]$adj.P.Val <- p.adjust(df.traits[[trait]]$p.value, method = "BH")
}
df <- do.call(rbind.data.frame, df.traits)
df$FDR <- -log10(df$adj.P.Val)
df$FDR[df$FDR <= -log10(0.05)] <- NA
df <- rbind.data.frame(df, c("Vagina", "Sex", "Epithelial_cells", NA, NA, NA))
df <- rbind.data.frame(df, c("Vagina", "Sex", "Keratinocytes", NA, NA, NA))
df <- rbind.data.frame(df, c("Prostate", "Sex", "Epithelial_cells", NA, NA, NA))
rownames(df) <- NULL
df$tissue_abbrv <- sapply(df$Tissue, function(i) tissue_info[tissue_info$tissue_ID==i, "tissue_abbrv"])
df$tissue_abbrv <- factor(df$tissue_abbrv, levels = tissue_info$tissue_abbrv, order = T)
df <- df[order(df$tissue_abbrv),]
df$CellType <- factor(df$CellType, levels = cell_types, order = T)
df <- df[order(df$CellType),]
d <- do.call(cbind.data.frame, lapply(traits, function(trait) df[df$Trait==trait,c("tissue_abbrv", "CellType","FDR")]))
rownames(d) <- paste0(d[,1], "-", d[,2])
D <- apply(d[,c(3,6,9,12)], 2, as.numeric)
colnames(D) <- traits
rownames(D) <- rownames(d)

# Figure S7H (left) ----
column_bottom_ha <- HeatmapAnnotation("Demographic traits" =  traits,
                                      col = list("Demographic traits" = traits_cols),
                                      show_legend = F, 
                                      show_annotation_name = F,
                                      annotation_name_rot = 90,
                                      annotation_name_gp = gpar(fontsize = 10),
                                      simple_anno_size = unit(0.3, "cm"),
                                      which = "column")
col_fun <- colorRamp2(seq(1.3, 7, length.out=6),
                                brewer.pal(9, "Reds")[4:9])
pdf(paste0(plot_path, "Figure_S7H.xCell_and_demographic_traits.pdf"), 
    width = 3, height = 4)
Heatmap(D[apply(D, 1, function(x) sum(!is.na(x)))>0,],
        na_col = "white",
        name = "FDR (-log10)",
        cluster_rows = F, 
        cluster_columns = F,
        row_names_side = "left",
        col = col_fun,#brewer.pal(9, "Reds")[c(4:9)],
        bottom_annotation = column_bottom_ha,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))
dev.off()


# Figure S7H (right) ----
# function to plot N in box plots ---
get_box_stats <- function(y, upper_limit = max(y) * 1) {
  return(data.frame(
    y = 1.05 * upper_limit,
    label = paste(
      "n =", length(y), "\n"
    )
  ))
}

# examples
# SKINS - Epithelial cells
data$SkinNotSunExposedSuprapubic$Ancestry <- gsub("EUR", "EA", data$SkinNotSunExposedSuprapubic$Ancestry)
data$SkinNotSunExposedSuprapubic$Ancestry <- gsub("AFR", "AA", data$SkinNotSunExposedSuprapubic$Ancestry)
data$SkinNotSunExposedSuprapubic$Ancestry <- factor(data$SkinNotSunExposedSuprapubic$Ancestry, levels = c("EA", "AA"), order = T)

p1 <- ggplot(data = data$SkinNotSunExposedSuprapubic,
       aes(x = Ancestry,
           y = Epithelial_cells,
           fill = Ancestry)) +
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
  ylab("xCell score") +
  scale_fill_manual(values = c("#E69F00", "#E69F00")) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  labs(title="SKINS - Epithelial cells") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")

# ADPSBQ - Adipocytes
data$AdiposeSubcutaneous$Age_int <- data$AdiposeSubcutaneous$Age
data$AdiposeSubcutaneous$Age <- sapply(data$AdiposeSubcutaneous$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
data$AdiposeSubcutaneous$Age <- factor(data$AdiposeSubcutaneous$Age, 
                                       levels = c("[20-45)", "[45-70]"),
                                       order = T)

p2 <- ggplot(data = data$AdiposeSubcutaneous,
       aes(x = Age,
           y = Adipocytes,
           fill = Age)) +
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
  ylab("xCell score") +
  scale_fill_manual(values = c("#56B4E9", "#56B4E9")) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  labs(title="ADPSBQ - Adipocytes") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")

#   LIVER - Adipocytes
pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}
data$Liver$BMI_int <- data$Liver$BMI
data$Liver$BMI <- sapply(data$Liver$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
data$Liver$BMI <- factor(data$Liver$BMI,
                 levels = c("Normal", "Overweight", "Obese"),
                 order = T)

p3 <- ggplot(data = data$Liver,
       aes(x = BMI,
           y = Adipocytes,
           fill = BMI)) +
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
  ylab("xCell score") +
  scale_fill_manual(values = c("#CC79A7", "#CC79A7", "#CC79A7")) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  labs(title="ADPSBQ - Adipocytes") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")

pdf(paste0(plot_path, "Figure_S7H.xCell_and_demographic_traits.examples.pdf"), 
    width = 3, height = 9)
ggarrange(p1, p2, p3, ncol = 1)
dev.off()
