#!/usr/bin/env Rscript

path_to_data <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Jose/"
setwd(path_to_data)
plot_path <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Jose/03_Models/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(reshape2)

# metadata ----
dir <- paste0("03_Models/Tissues/")

tissue_info <- readRDS("00_Data/Tissue_info.rds")
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
sex_tissues <- c("Uterus","Vagina","Ovary","Testis","Prostate")

# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38. PAR genes excluded
gene_annotation <- read.delim("00_Data/gencode.v26.GRCh38.genes.bed", header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

#Diseases we analyse
main_data <- read.csv("00_Data/main_final.csv") #read from here the histology data, and then all signal tissues for both diabetes
acronyms <- main_data$Acronym[-c(1:2)]
tissues <- main_data$Tissue[-c(1:2)]

tissues_for_diabetes <- read.table("00_Data/tissues_to_keep.txt")
tissues <- c(tissues_for_diabetes[,1], tissues_for_diabetes[,1], tissues)
acronyms <- c(rep( "MHT1D", length(tissues_for_diabetes[,1])), rep( "MHT2D", length(tissues_for_diabetes[,1])), acronyms)

tissue_names <- sapply(tissues, function(tissue) tissue_info$tissue_abbrv[tissue_info$tissue_ID==tissue]) #acronyms
names(tissue_names) <- tissue_names
my_names <- paste0(acronyms, " - ", tissue_names)
#Changing the names to prettiers ones
my_names <- sapply(my_names, function(name) 
  if(grepl("MHT1D", name)){gsub("MHT1D", "Type 1 diabetes", name)
  } else if(grepl("MHT2D", name)){gsub("MHT2D", "Type 2 diabetes", name)
  } else if(grepl("Atherosclerotic_Arteries - ARTTBL", name)){"Atherosclerosis - ARTTBL"
  } else if(grepl("Pneumonia_Hist - LUNG", name)){"Pneumonia - LUNG"
  } else{ name
  })
names(my_names) <- NULL

# differential expression analysis ----
dea <- readRDS("Tables/00.differential_expression_analysis.clinical_traits.rds")

# list of DEGs de per tissue and trait
get_deg <- function(tissue,trait){
  if(length(dea[[tissue]][[trait]])==1){
    return(NA)
  }else{
    de_genes <- rownames(dea[[tissue]][[trait]])[ dea[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(de_genes)
  }
}

# genes DE with each trait in each tissue
traits <- c("Ancestry", "Sex", "Age", "BMI", unique(acronyms))
DEGs <- lapply(unique(tissues), function(tissue)
  lapply(traits, function(trait) 
    get_deg(tissue,trait)
  )
)
names(DEGs) <- unique(tissues)
for(tissue in unique(tissues)){
  names(DEGs[[tissue]]) <- traits
  DEGs[[tissue]] <- DEGs[[tissue]][!sapply(DEGs[[tissue]],is.null)]}

# Is there a bias towards a particular direction of change ----
Xsq_fun <- function(tissue, trait1, trait2){
  #print(paste0(tissue, ": ", trait1, "-", trait2))
  if(tissue %in% sex_tissues & trait1 == "Sex"){
    return(list("twenty"= F,
                "P-value" = NA, 
                "O/E" = rep(NA, 4),
                "counts" = rep(NA, 4)))
  }else{
    # Do we observe a higher than expected overlap of DEGs in a particular direction of change?
    # a numeric vector representing the observed proportions
    # a vector of probabilities (of the same length of the observed proportions) representing the expected proportions
    trait1.up <- rownames(dea[[tissue]][[trait1]][dea[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                        dea[[tissue]][[trait1]]$logFC > 0,])
    trait1.down <- rownames(dea[[tissue]][[trait1]][dea[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                          dea[[tissue]][[trait1]]$logFC < 0,])
    
    trait2.up <- rownames(dea[[tissue]][[trait2]][dea[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                        dea[[tissue]][[trait2]]$logFC > 0,])
    trait2.down <- rownames(dea[[tissue]][[trait2]][dea[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                          dea[[tissue]][[trait2]]$logFC < 0,])
    
    # Observed counts
    counts <- c(sum(trait1.up %in% trait2.up), # upup
                sum(trait1.down %in% trait2.up), # downup
                sum(trait1.up %in% trait2.down), # updown
                sum(trait1.down %in% trait2.down) # downdown
    )
    # Expected proportions
    trait1.up.p <- length(trait1.up)/length(c(trait1.up, trait1.down))
    trait1.down.p <- length(trait1.down)/length(c(trait1.up, trait1.down))
    trait2.up.p <- length(trait2.up)/length(c(trait2.up, trait2.down))
    trait2.down.p <- length(trait2.down)/length(c(trait2.up, trait2.down))
    expected_prob <- c(trait1.up.p * trait2.up.p, trait1.down.p * trait2.up.p, trait1.up.p * trait2.down.p, trait1.down.p * trait2.down.p)
    expected_counts <- round(sum(counts) * expected_prob)
    
    # Return results
    if(sum(counts) < 20){
      print(paste0(tissue, ": Fewer than 20 genes DE with ", trait1, " and ", trait2))
      return(list("twenty" = F,
                  "P-value" = NA,
                  "O/E" = rep(NA, 4),
                  "counts" = counts))
      break
    }else{
      if(min(expected_counts) < 5){
        print(paste0(tissue, ": Number of observations is not enough for Chi-Square Test\nUsing Monte Carlo simulations"))
        Xsq <- chisq.test(counts, 
                          p = expected_prob,
                          simulate.p.value = T)
      }else{
        Xsq <- chisq.test(counts,
                          p = expected_prob)
      }
      oe <- Xsq$observed/round(Xsq$expected)
      return(list("twenty" = T,
                  "P-value" = Xsq$p.value,
                  "O/E" = oe,
                  "counts" = counts))
    }
  }
}

# Enrichment analysis  --
variables <- c("Age", "Ancestry", "Sex", "BMI")
Xsq_results <- lapply(variables, function(pw)
  lapply(1:length(tissues), function(i) 
    Xsq_fun(tissues[i], pw, acronyms[i])
  ))
names(Xsq_results) <- variables
for(pw in variables){names(Xsq_results[[pw]]) <- my_names}

# Enrichment statistics --
# P-values
p_values <- lapply(variables, function(i)
  sapply(my_names, function(tissue)
    Xsq_results[[i]][[tissue]][["P-value"]]
  )
)
names(p_values) <- variables
Xsq_pvalues <- p_values

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p_values), method = "BH") #NAs are not considered

adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(my_names), ncol = length(variables),
                         byrow = F)
fdr <- -log10(adjPVal_matrix)
colnames(fdr) <- variables
rownames(fdr) <- my_names
fdr[fdr <= -log10(0.05)] <- NA # if not tested NA (grey in plot); if FDR >= 0.05, 0 (white in plot)
sum(apply(fdr, 2, function(x) sum(!is.na(x))))
Xsq_tissues <- apply(fdr, 2, function(x) my_names[which(!is.na(x))])


i <- "Age"
data1 <- as.matrix(sapply(my_names, function(tissue) Xsq_results[[i]][[tissue]][["counts"]]))
rownames(data1) <- c("old - disease", "young - disease", "old - healthy", "young - healthy")
trues <- colSums(data1)>=20 & !is.na(fdr[,1]) #Excluding tissues with less than 20 genes with additive effects, and the ones that are not significant (hyperplasia in thyroid is not significant)
data1 <- as.data.frame(data1[,trues])
data1$type <- rownames(data1)
data2 <- as.data.frame(sapply(my_names, function(tissue) Xsq_results[[i]][[tissue]][["counts"]]/sum(Xsq_results[[i]][[tissue]][["counts"]])))
rownames(data2) <- c("old - disease", "young - disease", "old - healthy", "young - healthy")
data2 <- data2[,trues]
data2$type <- rownames(data2)

data <- melt(data1)
data$y <- 100*melt(data2)[,3]
#Order by FDR
sorting <- rownames(fdr)[order(fdr[,1])]
data$variable <- factor(data$variable, levels=sorting, order=T)
data$variable <- droplevels(data$variable)
data$type <- factor(data$type, levels = rev(c("old - disease", "young - disease", "old - healthy", "young - healthy")), order = T)
 
#To main (the three we replicate when downsampling)
data_subset <- data[data$variable %in% c("Type 2 diabetes - NERVET", "Type 1 diabetes - NERVET", "Hashimoto - THYROID"),]
age_subset <- ggplot(data_subset,
              aes(x = variable, 
                  y = y, 
                  fill = type,
                  label = value) ) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = rev(c("#00b159", "#00aedb", "#f37735", "#ffc425"))) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(vjust = 0.5)) +
  ylab("DEGs (%)") + xlab("")

png("Figures/Bias/Age_main.png", units="in", res=200,
    height = 2, width=6)
pdf("Figures/Bias/Age_main.pdf",
    height = 1.5, width=6)
age_subset
dev.off()