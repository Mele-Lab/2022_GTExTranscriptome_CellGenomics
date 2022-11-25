#!/usr/bin/env Rscript
rm(list=ls())

library(ggplot2)

first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/"
# first_dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"

folder <- c("03_Models/Tissues/")
acronyms <- c("MHT1D", "MHT2D")
diseases <- c("Type 1 diabetes", "Type 2 diabetes")
tissues <- list.dirs(paste0(first_dir, "Jose/", folder), recursive = F, full.names = F)

tissue_info <- readRDS(paste0(first_dir, "Jose/00_Data/Tissue_info.rds")) #To translate to colours later on
tissue_cols <- sapply(tissues, function(tissue) tissue_info$colcodes[tissue_info$tissue_ID==tissue])
abbreviations <- sapply(tissues, function(tissue) tissue_info$tissue_abbrv[tissue_info$tissue_ID==tissue])
names(tissue_cols) <- abbreviations

gene_annotation <- read.delim(paste0(first_dir,"Jose/00_Data/gencode.v26.GRCh38.genes.bed"), header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")


#Tissue Sharing
final_table <- data.frame(gene="1", tissue="1", logFC=1, disease="1")

#DEG
n_tissues <- c()
for(i in 1:length(acronyms)){ 
  print(acronyms[i])
  disease <- diseases[i]
  genes <- data.frame(gene="1", tissue="1", logFC=1, disease="1")

  c <- 0
  for(tissue in tissues){
    files <- list.files(paste0(first_dir, "/Jose/", folder, tissue, "/"))
    files <- files[grep("voom_limma", files)]
    files <- files[grep("interactions", files, invert = T)]
    file <- files[grep(acronyms[i], files)]
    if(identical(file, character(0))){
      next
    }
    file <- paste0(first_dir, "/Jose/", folder, tissue, "/", file)

    if(file.exists(file)){
      c <- c + 1
      dea_res <- readRDS(file)
      if(is.null(nrow(dea_res[[paste0(acronyms[i], 1)]]))){
        next
      }
      for(gene in rownames(dea_res[[paste0(acronyms[i], 1)]][dea_res[[paste0(acronyms[i], 1)]]$adj.P.Val<0.05,])){
        genes <- rbind(genes, c(gene=gene, tissue=tissue, logFC=dea_res[[paste0(acronyms[i], 1)]][gene,1], disease=disease))
      }
    }
  }
  n_tissues <- c(n_tissues, c)
  genes <- genes[-1,]
  final_table <- rbind(final_table, genes)
  
}

final_table <- final_table[-1,]

library(plyr) #Counting
counts <- ddply(final_table, .(gene, disease), nrow)
names(counts) <- c("gene", "disease", "number")

#Save the table for the supplements:
count_1 <- counts[counts$number>1,]
count_1 <- merge(count_1, gene_annotation, by="gene")
genes_ordered <- unique(count_1$symbol[order(count_1$number, decreasing = T)]) #Ordering gene levels so that the order is the number of tissues
data <- merge(final_table, count_1, by=c("gene", "disease"))
data$logFC <- as.numeric(data$logFC)
data$symbol <- factor(data$symbol, levels=genes_ordered)
data <- merge(data, tissue_info, by.x = "tissue", by.y = "tissue_ID") #To get the tissue abbreviation

gene_tissue <- data.frame("gene"="", "tissues"="", disease="")
for(i in 1:length(unique(count_1$gene))){
  gene <- unique(count_1$gene)[i] #i <- 388
  #iterate now over diseases
  for(d in unique(data[data$gene==gene,"disease"])){
    disease_subset <- data[data$disease==d,]
    my_tissues <- disease_subset[disease_subset$gene==gene,"tissue"]
    my_tissues <- paste0(my_tissues, collapse="", sep = ";")
    my_tissues <- substr(my_tissues,1,nchar(my_tissues)-1)
    gene_tissue <- rbind(gene_tissue, c(gene, my_tissues, d))
  }
}
colnames(gene_tissue) <- c("gene", "tissues", "disease")
to_save <- merge(count_1, gene_tissue, by=c("gene", "disease"))
to_save <- to_save[,c("gene", "symbol", "disease", "number", "tissues")]
colnames(to_save) <- c("gene", "gene_name", "disease", "n", "tissues")
to_save <- to_save[order(to_save$n, decreasing = T),] #Sort by n

#Add citations manually
to_save$reference <- ""
to_save$reference[1] <- "https://pubmed.ncbi.nlm.nih.gov/19296078/"
to_save$reference[3] <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4599009/"
to_save$reference[9] <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5414755/"
to_save$reference[11] <- "https://pubmed.ncbi.nlm.nih.gov/28065675/ https://pubmed.ncbi.nlm.nih.gov/26188370/"
to_save$reference[12] <- "https://pubmed.ncbi.nlm.nih.gov/34484123/"
to_save$reference[13] <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5605719/"
to_save$reference[14] <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8074411/"

#Order by n
to_save <- to_save[order(to_save$n, to_save$disease, decreasing = T),] #Sort by n

#Add mean TPM per gene and tissue:
#First, read tpm per tissue
for(tissue in unique(data$tissue)){
  print(tissue)
  #create a variable with the tissue name that has the tpm means of the tissue
  tpm <- readRDS(paste0(first_dir, "Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))
  colnames(tpm) <- sapply(colnames(tpm), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  
  files <- list.files(paste0(first_dir, "Jose/03_Models/Tissues/", tissue, "/"))
  metadata_file <- grep("SampleMetadata", files, value=T)
  metadata <- readRDS(paste0(first_dir, "Jose/03_Models/Tissues/", tissue, "/", metadata_file))
  
  means <- cbind(round(rowMeans(tpm[,colnames(tpm) %in% metadata[metadata$MHT1D==0, "Sample"]]), 2), 
                round(rowMeans(tpm[,colnames(tpm) %in% metadata[metadata$MHT1D==1, "Sample"]]), 2),
                round(rowMeans(tpm[,colnames(tpm) %in% metadata[metadata$MHT2D==0, "Sample"]]), 2),
                round(rowMeans(tpm[,colnames(tpm) %in% metadata[metadata$MHT2D==1, "Sample"]]), 2))
  colnames(means) <- c("Healthy_1", "T1D", "Healthy_2", "T2D")
  
  assign(tissue, means)
}

to_save$healthy <- NA
to_save$diabetes <- NA
for(i in 1:nrow(to_save)){
  gene <- to_save[i,"gene"]
  if(to_save[i,"disease"]=="Type 1 diabetes"){column <- 1}else{column <- 3}
  my_tissues <- to_save$tissues[i]
  my_tissues <- strsplit(my_tissues, ";")[[1]]
  tpms_healthy <- c()
  tpms_diabetes <- c()
  for(tissue in my_tissues){
    tpms_healthy <- c(tpms_healthy, eval(as.symbol(tissue))[gene,column])
    tpms_diabetes <- c(tpms_diabetes, eval(as.symbol(tissue))[gene,column+1])
  }
  to_save$healthy[i] <- paste(tpms_healthy, collapse = ";")
  to_save$diabetes[i] <- paste(tpms_diabetes, collapse = ";")
}

names(to_save) <- c("ensembl id", "gene name", "clinical trait", "number of tissues", "tissues", "reference", "mean healthy TPM", "mean diabetic TPM" )
library("xlsx")
write.xlsx(as.data.frame(to_save), paste0(first_dir, "Jose/Tables/Table_S6F.xlsx"),
           col.names = TRUE, row.names = F, append = FALSE)


#Plot
counts <- counts[counts$number>2,]
counts <- merge(counts, gene_annotation, by="gene")
genes_ordered <- counts$symbol[order(counts$number, decreasing = T)] #Ordering gene levels so that the order is the number of tissues
data <- merge(final_table, counts, by=c("gene", "disease"))
data$logFC <- as.numeric(data$logFC)
data$symbol <- factor(data$symbol, levels=genes_ordered)
data <- merge(data, tissue_info, by.x = "tissue", by.y = "tissue_ID") #To get the tissue abbreviation

colors <- tissue_cols[names(tissue_cols) %in% unique(data$tissue_abbrv)]
g1 <- ggplot(data=data, aes(symbol, logFC)) +
  geom_point(aes(col=tissue_abbrv)) + geom_hline(yintercept=0, lty=2) + theme_bw() +
  facet_grid(.~ disease, drop = T, scales = "free") +
  scale_color_manual(values=colors) +
  theme(#panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.5),
        axis.ticks = element_line(size = 0.2)) + xlab("")

g1
library(grid)
gt = ggplot_gtable(ggplot_build(g1))
gt$widths[7] = 1.75*gt$widths[7] #To change the facet width
grid.draw(gt)

#Shared in at least 3 tissues
pdf(paste0(first_dir, "Jose/03_Models/TSA/Top_Shared.pdf"),
    width = 6, height = 3.5)
png(paste0(first_dir, "Jose/03_Models/TSA/scatterplot_3.png"),
    units = "in", width = 6, height = 3.5, res=200)
grid.draw(gt)
dev.off()
