rm(list=ls())

first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/"
# first_dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"

dir <- paste0(first_dir, "Jose/03_Models/Tissues/")

gene_annotation <- read.delim(paste0(first_dir, "Jose/00_Data/gencode.v26.GRCh38.genes.bed"), header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")



#Label the genes associated to diabetes in TWAS:
library("readxl")
type_2_files <- c()
for(file in list.files(paste0(first_dir, "Jose/03_Models/TWAS/type_2/"), full.names = T)){
  type_2_files <- rbind(type_2_files, read_excel(file))
}

type_1_files <- c()
for(file in list.files(paste0(first_dir, "Jose/03_Models/TWAS/type_1/"), full.names = T)){
  type_1_files <- rbind(type_1_files, read_excel(file))
}

#Diabetes in tibial nerve
type_2_genes <- unique(type_2_files$Gene)
type_2_files_nerve <- type_2_files[type_2_files$Tissue=="Nerve Tibial",]
type_2_files_nerve_genes <- type_2_files_nerve$Gene

metadata <- readRDS(paste0(dir, "NerveTibial/NerveTibial.SampleMetadata.MHT1D_MHT2D.rds"))

dea_res <- readRDS(paste0(dir, "NerveTibial/NerveTibial.MHT1D_MHT2D.voom_limma.results_PEER.rds"))


#Scatter plot

type_1 <- dea_res$MHT1D[dea_res$MHT1D$adj.P.Val<0.05,]
type_1 <- as.data.frame(cbind(rownames(type_1), type_1$logFC))
colnames(type_1) <- c("gene", "logFC")
type_2 <- dea_res$MHT2D[dea_res$MHT2D$adj.P.Val<0.05,]
type_2 <- as.data.frame(cbind(rownames(type_2), type_2$logFC))
colnames(type_2) <- c("gene", "logFC")
age <- dea_res$Age[dea_res$Age$adj.P.Val<0.05,]
age <- as.data.frame(cbind(rownames(age), age$logFC))
colnames(age) <- c("gene", "logFC")

#Type 2 and age
data <- merge(age, type_2, by="gene")
data$logFC.x <- as.numeric(data$logFC.x)
data$logFC.y <- as.numeric(data$logFC.y)

data$color <- "grey60"
data$color[data$logFC.x>0 & data$logFC.y>0] <- "#B2182B" #rev(brewer.pal(11, "RdBu")[c(2, 10)]) 
data$color[data$logFC.x<0 & data$logFC.y<0] <- "#2166AC" #library(RColorBrewer)

colors <- c("grey60", "#B2182B", "#2166AC") #I know there is a better way to do it, but this will do the job by now
names(colors) <- c("grey60", "#B2182B", "#2166AC")

data <- merge(data, gene_annotation, by="gene")

library(ggplot2)
library(ggrepel)

data_out <- data[data$symbol %in% type_2_files_nerve_genes,]

g <- ggplot(data) +
  geom_point(aes(logFC.x, logFC.y, col=color), alpha=0.6) + theme_bw() + 
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlab("Age (logFC)") +
  ylab("Type 2 diabetes (logFC)") +
  scale_color_manual(values=colors) +
  theme(legend.position = "none") + 
  geom_text_repel(data=data_out, aes(logFC.x, logFC.y, label=symbol), max.overlaps = 200) 
g

png(paste0(first_dir, "Jose/03_Models/Disease_Sharing/", "Age_type_2_TWAS", ".png"),
    units = "in", width = 5, height = 5, res = 200)
g
dev.off()
