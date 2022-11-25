#!/usr/bin/env Rscript

# Libraries ----
library(ggplot2)
library(ggtext)
library(ggpubr)

# function to plot N in box plots ---
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"#,
      #"Mean =", round(mean(y), 2), "\n",
      #"Median =", round(median(y), 2), "\n"
    )
  ))
}

# ---- Data ---- ####

# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" = "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
tissue_info <- readRDS("~/GTEx_v8/Raquel/Draft/Data/Tissue_info.46_tissues.rds") # tissues ordered by sample size
tissues <- tissue_info$tissue_ID

# Tissue metadata ----
# Number of samples
mdata <- lapply(tissues, function(tissue)
  readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds")))
names(mdata) <- tissues

# Gene Annotation ----
gene_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")

# 1. Hier.part:expression ----
hier_part <- readRDS("~/GTEx_v8/Raquel/Draft/Analysis/Expression.Hier.part/Data/Exprs.Hier.part.rds")
head(hier_part$Spleen)
# All genes expressed in tissue 
hier.part.exprs <- lapply(tissues, function(tissue)
  readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/",tissue,"/", tissue,".hier_part.rds")))
names(hier.part.exprs) <- tissues

# 2. Differential expression results ----
dea_res <- lapply(tissues, function(tissue)
  readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/",tissue,"/", tissue,".voom_limma.covariates_and_traits.results.rds")))
names(dea_res) <- tissues
degs <- readRDS("~/GTEx_v8/Raquel/Draft/01.DiffExprs/Data/Genes_DE.rds")

# AdiposeSubcutaneous: ST5 ----
100*hier.part.exprs$AdiposeSubcutaneous[Reduce(intersect,degs$AdiposeSubcutaneous),c(1,6,7,8,9)]
gene_annotation[gene_annotation$ensembl.id=="ENSG00000166444.18",]
tissue <- "AdiposeSubcutaneous"
gene_name <- "ST5"
gene <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]


# Gene TPM ----
tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))

# Gene residuals ----
residuals <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residuals.rds"))

pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}
get.gene.data <- function(gene_name, tpm, res = residuals){
  #  Expresidualssion data ----
  #tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
  #residuals <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residualsiduals.rds"))
  identical(colnames(tpm), colnames(residuals))
  
  gene <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
  df <- cbind.data.frame(tpm[gene,])
  colnames(df) <- "TPM"
  df$residuals<- residuals[gene,]
  df$Age_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Age"])
  df$BMI_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"BMI"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  df$Ancestry <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$BMI <- sapply(df$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
  df$BMI <- factor(df$BMI,
                   levels = c("Normal", "Overweight", "Obese"),
                   order = T)
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$Gene <- rep(gene, nrow(df))
  head(df)
  return(df)
}

# get gene data ----
data <- get.gene.data(gene_name, tpm, residuals)
cols <- c(traits_cols["Ancestry"], traits_cols["Ancestry"])
names(cols) <- c("EA", "AA")

# plot 1 ----
p1 <- ggplot(data = data,
             aes(x = Ancestry,
                 y = log2(1+TPM),
                 fill = Ancestry),
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
  labs(title=gene_name) +
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

# plot 2 ----
d <- data.frame("variable" = "Ancestry",
                "value"  = 100*hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], "Ancestry_abs"])
variable_col <- traits_cols["Ancestry"]
names(variable_col) <- "Ancestry"

p2 <- ggplot(d, aes(x=variable, y=value)) + 
  geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
  scale_fill_manual(values = variable_col) +
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
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

# pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/demographic_traits_and_cis_driven/figures/Figure_2D.Ancestry_Spleen_ACKR1.pdf",
#     width = 3, height = 4)
plot1 <- ggarrange(p1, p2, nrow = 2, heights =  c(9,2))
# dev.off()


# # ColonSigmoid:  LINC01597 ----
# tissue <- "ColonSigmoid"
# gene_name <- "LINC01597"
# gene <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
# sapply(traits, function(trait) gene %in% degs[[tissue]][[trait]])
# 
# # Gene TPM ----
# tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
# 
# # Gene residuals ----
# residuals <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residuals.rds"))

# get gene data ----
data <- get.gene.data(gene_name, tpm, residuals)
cols <- c(traits_cols["Sex"], traits_cols["Sex"])
names(cols) <- c("Male", "Female")

p3 <- ggplot(data = data,
             aes(x = Sex,
                 y = log2(1+TPM),
                 fill = Sex),
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
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
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

#hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], ]
d <- data.frame("variable" = "Sex",
                "value"  = 100*hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], "Sex_abs"])
variable_col <- traits_cols["Sex"]
names(variable_col) <- "Sex"

p4 <- ggplot(d, aes(x=variable, y=value)) + 
  geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
  scale_fill_manual(values = variable_col) +
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
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
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank()) 

#pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/demographic_traits_and_cis_driven/figures/Figure_2D.Sex_ColonSigmoid_LINC01597.pdf",
#    width = 3, height = 4)
plot2 <- ggarrange(p3, p4, nrow = 2, heights =  c(9,2))

#dev.off()

# # ArteryAorta: ROBO2 ----
# tissue <- "ArteryAorta"
# gene_name <- "ROBO2"
# gene <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
# sapply(traits, function(trait) gene %in% degs[[tissue]][[trait]])
# 
# # Gene TPM ----
# tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
# 
# # Gene residuals ----
# residuals <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residuals.rds"))
# 
# # get gene data ----
# data <- get.gene.data(gene_name, tpm, residuals)
cols <- c(traits_cols["Age"], traits_cols["Age"])
names(cols) <- c("[20-45)", "[45-70]")

p5 <- ggplot(data = data,
             aes(x = Age,
                 y = log2(1+TPM),
                 fill = Age),
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
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
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

#hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], ]
d <- data.frame("variable" = "Age",
                "value"  = 100*hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], "Age_abs"])
variable_col <- traits_cols["Age"]
names(variable_col) <- "Age"

p6 <- ggplot(d, aes(x=variable, y=value)) + 
  geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
  scale_fill_manual(values = variable_col) +
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
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
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank()) 

#pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/demographic_traits_and_cis_driven/figures/Figure_2D.Age_ArteryAorta_ROBO2.pdf",
#    width = 3, height = 4)
plot3 <- ggarrange(p5, p6, nrow = 2, heights =  c(9,2))
#dev.off()


# # AdiposeSubcutaneous: ROBO2 ----
# tissue <- "AdiposeSubcutaneous"
# gene_name <- "SLC27A2"
# gene <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
# sapply(traits, function(trait) gene %in% degs[[tissue]][[trait]])

# Gene TPM ----
# tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
# 
# # Gene residuals ----
# residuals <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residuals.rds"))
# 
# # get gene data ----
# data <- get.gene.data(gene_name, tpm, residuals)
cols <- c(traits_cols["BMI"], traits_cols["BMI"], traits_cols["BMI"])
names(cols) <- c("Normal", "Overweight", "Obese")

p7 <- ggplot(data = data,
             aes(x = BMI,
                 y = log2(1+TPM),
                 fill = BMI),
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
  ylab("TPM (log<sub>2</sub>") +
  scale_fill_manual(values = cols) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title=gene_name) +
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

#hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], ]
d <- data.frame("variable" = "BMI",
                "value"  = 100*hier.part.exprs[[tissue]][gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"], "BMI_abs"])
variable_col <- traits_cols["BMI"]
names(variable_col) <- "BMI"

p8 <- ggplot(d, aes(x=variable, y=value)) + 
  geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
  scale_fill_manual(values = variable_col) +
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
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
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank()) 

#pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/demographic_traits_and_cis_driven/figures/Figure_2D.BMI_AdiposeSubcutaenous_SLC27A2.pdf",
#    width = 3, height = 4)
plot4 <- ggarrange(p7, p8, nrow = 2, heights =  c(9,2))
#dev.off()

#pdf("~/GTEx_v8/Raquel/Draft/Analysis/expression/additive_effects_are_widespread_and_tissue_specific/figures/Figure_S3.genes_with_4_additive_effecsts.pdf",
      #width = 8, height = 10)
pdf('~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/Figure_S3D.pdf',
    width = 8, height = 10)
ggarrange(plot1, plot2, plot3, plot4, nrow=2, ncol = 2)
dev.off()
