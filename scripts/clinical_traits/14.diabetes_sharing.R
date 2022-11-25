#!/usr/bin/env Rscript
rm(list=ls())

library(WebGestaltR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE)
library(vroom)
library(tidyr)
library(msigdbr)
library(xCell)
library(GSEABase)


first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/"
# first_dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"

dir <- paste0(first_dir, "Jose/03_Models/Tissues/")


gene_annotation <- read.delim(paste0(first_dir, "Jose/00_Data/gencode.v26.GRCh38.genes.bed"), header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")



#Diabetes in tibial nerve

# load("~/Documents/mn4/Raquel/R_functions/Funcional_enrichments.RData")
metadata <- readRDS(paste0(dir, "NerveTibial/NerveTibial.SampleMetadata.MHT1D_MHT2D.rds"))

dea_res <- readRDS(paste0(dir, "NerveTibial/NerveTibial.MHT1D_MHT2D.voom_limma.results_PEER.rds"))


#Common DEG:      Upset
DEG <- list()
DEG[["Type 1 Diabetes"]] <- rownames(dea_res$MHT1D[dea_res$MHT1D$adj.P.Val<0.05,])
DEG[["Type 2 Diabetes"]] <- rownames(dea_res$MHT2D[dea_res$MHT2D$adj.P.Val<0.05,])

library(VennDiagram)

venn <- venn.diagram(
  x = DEG,
  fill=c(gray.colors(4)[4], gray.colors(4)[3]),
  category.names = names(DEG),
  filename = NULL,
  output=TRUE,
  lwd = 1.5,
  cex = 2.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.pos=c(182,182),
  cat.cex=2, 
  disable.logging	= T
)

pdf(file=paste0(first_dir, "Jose/03_Models/Disease_Sharing/Diabetes_Nerve_Genes_Venn.pdf"))
grid.draw(venn)
dev.off()


#Scatter plot
type_1 <- dea_res$MHT1D[dea_res$MHT1D$adj.P.Val<0.05,]
type_1 <- as.data.frame(cbind(rownames(type_1), type_1$logFC))
colnames(type_1) <- c("gene", "logFC")
type_2 <- dea_res$MHT2D[dea_res$MHT2D$adj.P.Val<0.05,]
type_2 <- as.data.frame(cbind(rownames(type_2), type_2$logFC))
colnames(type_2) <- c("gene", "logFC")

data <- merge(type_1, type_2, by="gene")
data$logFC.x <- as.numeric(data$logFC.x)
data$logFC.y <- as.numeric(data$logFC.y)

data$color <- "#B2182B"
colors <- c("#B2182B") 
names(colors) <- c("#B2182B")
#Use ggrepel

data <- merge(data, gene_annotation, by="gene")
og <- data
load(paste0(first_dir, "Raquel/R_functions/Funcional_enrichments.RData")) #Functions to do ORA enrichments

#Now, gene_annotation has the entrez.ids
data <- merge(data, gene_annotation, by.x="gene", by.y="ensembl.id")

#Do enrichments
common_genes <- DEG[["Type 1 Diabetes"]][DEG[["Type 1 Diabetes"]] %in% DEG[["Type 2 Diabetes"]]]
dea_res <- readRDS(paste0(dir, "NerveTibial/NerveTibial.MHT1D_MHT2D.voom_limma.results_PEER.rds")) #Reading again after loading code for enrichments
type_1 <- dea_res$MHT1D[dea_res$MHT1D$adj.P.Val<0.05,]
type_2 <- dea_res$MHT2D[dea_res$MHT2D$adj.P.Val<0.05,]

up_up <- common_genes[common_genes %in% rownames(type_1[type_1$logFC>0,]) & common_genes %in% rownames(type_2[type_2$logFC>0,])]
down_down <- common_genes[common_genes %in% rownames(type_1[type_1$logFC<0,]) & common_genes %in% rownames(type_2[type_2$logFC<0,])]

#Enrichment of the up_up on one side, and down_down on the other
gl <- sapply(up_up, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
gl <- gl[!is.na(gl)] #159

#Using as background expressed genes in the tissue:
bg <- rownames(dea_res$Age)
bg <- sapply(bg, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
bg <- bg[!is.na(bg)] #14748

dbs <- c("GO:BP", "GO:MF","GO:CC", "KEGG", "ReactomePA","human_phenotype","H", "DO","DisGeNET","OMIM", "gwas","cell_marker", "xCell")
ora_results_up_up <- lapply(dbs, function(db) {
  ora.fun( gene.list = gl,
           bg.list = bg,
           db = db,
           pvalueCutoff = 0.05)
}
)
names(ora_results_up_up) <- dbs
up_up_mf <- ora_results_up_up$`GO:MF`
up_term <- up_up_mf[3,2]
up_genes <- up_up_mf$geneID
up_genes <- unique(unlist(strsplit(up_genes, "/")))

#Down_down
gl <- sapply(down_down, function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "entrez.id"])
gl <- gl[!is.na(gl)] #159
ora_results <- lapply(dbs, function(db) {
  ora.fun( gene.list = gl,
           bg.list = bg,
           db = db,
           pvalueCutoff = 0.05)
}
)
names(ora_results) <- dbs
down_down_mf <- ora_results$`GO:MF`
down_term_1 <- "ion channel activity"
down_genes <- down_down_mf$geneID
down_genes <- unique(unlist(strsplit(down_genes, "/")))
#ankyrin anchor ion channels, spectrin is also involved in the citoskeleton placement in neurons
#they both belong to cytoskeletal protein binding
# down_term_2 <- "cytoskeletal protein binding"



to_label <- data[data$entrez.id %in% c(down_genes,up_genes),]


library(ggrepel)
g <- ggplot(data) +
  geom_point(aes(logFC.x, logFC.y, col=color)) + theme_bw() +  #, alpha=0.7
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlab("Type 1 diabetes (logFC)") +
  ylab("Type 2 diabetes (logFC)") +
  xlim(-2,2) + ylim(-2,2) +
  scale_color_manual(values=colors) +
  theme(legend.position = "none") + 
  geom_text_repel(data=to_label, aes(logFC.x, logFC.y, label=gene.name), max.overlaps = 50) #+
  # geom_label(x=1.2, y=1.9, label=up_term) +
  # geom_label(x=-1.1, y=-1.9, label=down_term_1)

g

# pdf(paste0(first_dir, "Jose/03_Models/Disease_Sharing/logFC.pdf"),
#     width = 5, height = 5)
png(paste0(first_dir, "Jose/03_Models/Disease_Sharing/logFC.png"),
    units = "in", width = 5, height = 5, res = 200)
g
dev.off()


#Save table with the enrichments
up_up_to_save <- head(ora_results_up_up$`GO:MF`)
down_down_to_save <- head(ora_results$`GO:MF`)
#Change geneID to gene name
up_up_to_save$enrichment_ratio <- sapply(up_up_to_save$GeneRatio, function(i) as.numeric(unlist(strsplit(i, split = "/")))[[1]]/as.numeric(unlist(strsplit(i, split = "/")))[[2]] )/sapply(up_up_to_save$BgRatio, function(i) as.numeric(unlist(strsplit(i, split = "/")))[[1]]/as.numeric(unlist(strsplit(i, split = "/")))[[2]] )
down_down_to_save$enrichment_ratio <- sapply(down_down_to_save$GeneRatio, function(i) as.numeric(unlist(strsplit(i, split = "/")))[[1]]/as.numeric(unlist(strsplit(i, split = "/")))[[2]] )/sapply(down_down_to_save$BgRatio, function(i) as.numeric(unlist(strsplit(i, split = "/")))[[1]]/as.numeric(unlist(strsplit(i, split = "/")))[[2]] )
up_up_to_save <- up_up_to_save[order(up_up_to_save$enrichment_ratio, decreasing = T),]
down_down_to_save <- down_down_to_save[order(down_down_to_save$enrichment_ratio, decreasing = T),]

#Getting gene name
gene_annotation_subset <- na.omit(gene_annotation)
for(i in 1:nrow(up_up_to_save)){
  row <- up_up_to_save$geneID[i]
  genes <- strsplit(row, split="/")[[1]]
  print(genes)
  genes <- sapply(genes, function(gene) gene_annotation_subset[gene_annotation_subset$entrez.id==gene, "gene.name"])
  gene_names <- paste0(genes, sep="/", collapse = "")
  gene_names <- substr(gene_names, 1, nchar(gene_names)-1)
  up_up_to_save$gene_name[i] <- gene_names
}
for(i in 1:nrow(down_down_to_save)){
  row <- down_down_to_save$geneID[i]
  genes <- strsplit(row, split="/")[[1]]
  print(genes)
  genes <- sapply(genes, function(gene) gene_annotation_subset[gene_annotation_subset$entrez.id==gene, "gene.name"])
  gene_names <- paste0(genes, sep="/", collapse = "")
  gene_names <- substr(gene_names, 1, nchar(gene_names)-1)
  down_down_to_save$gene_name[i] <- gene_names
}

final_data <- rbind(up_up_to_save, down_down_to_save)
final_data <- final_data[,-8]

#Mention that it is GO:MF
library("xlsx")
write.xlsx(as.data.frame(final_data), "~/Documents/mn4/Jose/Tables/Table_S6I.xlsx", 
           col.names = TRUE, row.names = F, append = FALSE) #To data frame so that the column names are mantained
