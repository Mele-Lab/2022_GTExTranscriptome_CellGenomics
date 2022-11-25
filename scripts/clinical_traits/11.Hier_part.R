#!/usr/bin/env Rscript
library(ggplot2)

# ---- Data ---- ####
first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Jose/"
dir <- paste0(first_dir, "03_Models/Tissues/")
# dir <- "~/Documents/mn4/Jose/"
# Individual Traits ----
traits_cols <- c("Age" = "#56B4E9", "Ancestry" = "#E69F00", "Sex" = "#009E73", "BMI" = "#CC79A7", "Disease" = "#000000")
my_traits <- names(traits_cols)

# Tissues ----
tissue_info <- readRDS(paste0(first_dir, "00_Data/Tissue_info.rds")) # tissues ordered by sample size

tissues <- list.dirs(dir, full.names = F)[-1]

# Gene Annotation ----
gene_annotation <- read.delim(paste0(first_dir, "00_Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed"))
autosomal_genes <- gene_annotation[!gene_annotation$chr %in% c("chrX","chrY","chrM"),"ensembl.id"] 

# ---- Code ---- ####
folders <- tissues
# 1. Differential expression results ----
folders_with_expression <- c()
dea_res <- list()
names <- c()
for(folder in folders){
  print(folder)
  files <- list.files(paste0(dir, folder, "/"))
  files <- files[grep("voom_limma.results_PEER", files)]
  files <- grep(files, pattern = "no_disease", invert=TRUE, value = T)
  files <- grep(files, pattern = "interactions", invert=TRUE, value = T)
  for(file in files){ #Outdated, previously we had many models per tissue, now we don't
    data <- readRDS(paste0(dir, folder, "/", file))
    for(trait in names(data)){
      if(!trait %in% c("Age", "Ancestry1", "Sex1", "BMI")){
        data[[trait]]
        if(sum(data[[trait]]$adj.P.Val<0.05)>0){
          folders_with_expression <- c(folders_with_expression, folder)
          break
        }
      }
    }
    name <- paste(strsplit(file, ".", fixed="T")[[1]][1], strsplit(file, ".", fixed="T")[[1]][2], sep=".")
    dea_res[[name]] <- data
    names <- c(names, name)
  }
}


# folders_with_expression
write.table(folders_with_expression, paste0(first_dir,"/00_Data/tissues_to_keep.txt")) #Tissues with DEGs

# 2. Hier.part:expression ----
# All genes expressed in tissue 
hier.part.exprs <- list()
for(folder in folders_with_expression){
  print(folder)
  files <- list.files(paste0(dir, folder, "/"))
  files <- files[grep("hier_part", files)]
  files <- files[grep("spliced", files, invert = T)]
  for(file in files){
    data <- readRDS(paste0(dir, folder, "/", file))
    name <- paste(strsplit(file, ".", fixed="T")[[1]][1], strsplit(file, ".", fixed="T")[[1]][2], sep=".")
    hier.part.exprs[[name]] <- data
  }
}

#dea_res to the subset of dea res with signal:
names_with_expression <- sapply(names, function(name) strsplit(name, ".", fixed = T)[[1]][1] %in% folders_with_expression)
names <- names[names_with_expression]

dea_res <- dea_res[names]

for(name in names){   #Removing the 1 from the categorial variables (Sex1 -> Sex)
  labels <- names(dea_res[[name]])
  labels[substr(labels,nchar(labels),nchar(labels))==1] <-   substr(labels[substr(labels,nchar(labels),nchar(labels))==1], 1, nchar(labels[substr(labels,nchar(labels),nchar(labels))==1])-1)
  names(dea_res[[name]]) <- labels
}

# 3 DEG per tissue:trait  ----
# Keep only autosomal genes -- and DE
function_for_deg <- function(i, trait){
  genes <- rownames(dea_res[[names[i]]][[trait]][dea_res[[names[i]]][[trait]]$adj.P.Val<0.05,])
  genes <- genes[genes %in% autosomal_genes]
  return(genes)
}

de.genes <- lapply(1:length(names), function(i)
  lapply(names(dea_res[[names[i]]]), function(trait)
    function_for_deg(i, trait)
  )
)
names(de.genes) <- names

for(tissue in names){
  names(de.genes[[tissue]]) <- names(dea_res[[tissue]])
}

#Part 2: Barplots
get.deg <- function(i){
  genes <- unique(unlist(de.genes[[names[i]]]))
  return(genes)
}

deg <- lapply(1:length(names), function(i) get.deg(i))
names(deg) <- names

# 4.2 Subset hier.part data ---- This has already been done in the previous script
for(i in 1:length(names)){
  hier.part.exprs[[names[i]]] <- hier.part.exprs[[names[i]]][rownames(hier.part.exprs[[names[i]]]) %in% deg[[names[i]]],]
}


# 4.3 Create gene:trait adj.P.Val matrix ----
fdr.matrix <- lapply(1:length(names), function(i)
  do.call(cbind.data.frame,
          lapply(names(dea_res[[names[i]]]), function(trait)
            sapply(rownames(hier.part.exprs[[names[i]]]), function(gene)
              dea_res[[names[i]]][[trait]][gene,"adj.P.Val"]
            )
          )
  )
)
og_fdr.matrix <- fdr.matrix
# names(fdr.matrix) <- tissues[!tissues %in% sex_tissues]
names(fdr.matrix) <- names
for(i in 1:length(names)){
  colnames(fdr.matrix[[names[i]]]) <- names(dea_res[[names[i]]])
}

for(i  in 1:length(names)){
  fdr.matrix[[names[i]]][fdr.matrix[[names[i]]] >= 0.05] <- NA
}

for(i in 1:length(names)){
  hier.part.exprs[[names[i]]] <- hier.part.exprs[[names[i]]][,grepl("_abs", names(hier.part.exprs[[names[i]]]))]
}

for(i in 1:length(names)){
  boolean <- is.na(fdr.matrix[[names[i]]])
  final_names <- gsub('_abs', '', names(hier.part.exprs[[names[i]]]))
  names(hier.part.exprs[[names[i]]]) <- final_names
  hier.part.exprs[[names[i]]] <- hier.part.exprs[[names[i]]][names(hier.part.exprs[[names[i]]]) %in% colnames(boolean)]
  hier.part.exprs[[names[i]]][boolean] <- NA
}

##########################


# 4.5 Explained variance (ev) ----
get_expln_var <- function(tissue, hier.part.data, i){
  if(is.null(hier.part.data[[tissue]])){
    x <- rep(0, length(names(dea_res[[tissue]])))
    names(x) <- names(dea_res[[tissue]])
  } else{
    labels <- names(dea_res[[tissue]])

    x <- colSums(hier.part.data[[tissue]][,labels], na.rm = T)

  }
  return(x)
}

get_expln_var_demo <- function(tissue, hier.part.data, i){
  if(is.null(hier.part.data[[tissue]])){
    x <- rep(0, length(names(dea_res[[tissue]])[names(dea_res[[tissue]]) %in% c("Age", "Ancestry", "Sex", "BMI")]))
    names(x) <- names(dea_res[[tissue]])[names(dea_res[[tissue]]) %in% c("Age", "Ancestry", "Sex", "BMI")]
  } else{
    labels <- names(dea_res[[tissue]])[names(dea_res[[tissue]]) %in% c("Age", "Ancestry", "Sex", "BMI")]
    x <- colSums(hier.part.data[[tissue]][,labels], na.rm = T)
    
  }
  return(x)
}

ev.exprs.total <- lapply(1:length(names), function(i)
  get_expln_var(names[i], hier.part.exprs, i)
)
names(ev.exprs.total) <- names

# 4.6 Explained variance (ev) ; autosomal genes ----
hier.part.exprs.autosomal_genes <- lapply(1:length(names), function(i)
  hier.part.exprs[[names[i]]][rownames(hier.part.exprs[[names[i]]]) %in% autosomal_genes,])
names(hier.part.exprs.autosomal_genes) <- names

ev.exprs.total.aut_genes <- lapply(1:length(names), function(i)
  get_expln_var(names[i], hier.part.exprs.autosomal_genes, i)
)
names(ev.exprs.total.aut_genes) <- names

ev.exprs.total.aut_genes_demo <- lapply(1:length(names), function(i)
  get_expln_var_demo(names[i], hier.part.exprs.autosomal_genes, i)
)
names(ev.exprs.total.aut_genes_demo) <- names

# 4.8 Bar plot ordered by sample size ----
traits <- c()
my_tissues <- c()
for(i in 1:length(names)){
  my_names <- names(ev.exprs.total.aut_genes[[names[i]]])
  traits <- c(traits, my_names)
  tissue <- strsplit(names[i],".", fixed = T)[[1]][1]
  my_tissues <- c(my_tissues, rep(tissue, length(my_names)))
}


data_plot <- cbind.data.frame(traits,
                              my_tissues,
                              unlist(ev.exprs.total.aut_genes))   #The fraction (frac.) is not correctly computed for the diseases in the lung
colnames(data_plot) <- c("Trait","Tissue", "Value")                 
data_plot$Abbrv <- sapply(data_plot$Tissue, function(i) tissue_info[tissue_info$tissue_ID==i,"tissue_abbrv"])
data_plot$Trait <- as.factor(data_plot$Trait)#, levels=traits, order = T)

data_plot$Variable <- paste0(data_plot$Trait, " - ", data_plot$Abbrv)

data_plot$Trait_type <- data_plot$Trait
data_plot$Trait_type <- as.character(data_plot$Trait_type)
data_plot$Trait_type <- as.factor(data_plot$Trait_type)

data_plot <- data_plot[,-grep("Variable", colnames(data_plot))]

final_data <- reshape(data_plot, idvar = "Tissue", timevar = "Trait_type", direction = "wide")
final_data <- final_data[,c(1, grep("Value", colnames(final_data)))]
rownames(final_data) <- final_data[,1]
final_data <- final_data[,-1]
colnames(final_data) <- unique(data_plot$Trait)
final_data <- as.matrix(final_data)
final_data[is.na(final_data)] <- 0
final_data_frac <- final_data/rowSums(final_data)
final_data_frac <- final_data_frac[nrow(final_data_frac):1, ]

rownames(final_data_frac) <- sapply(rownames(final_data_frac), function(tissue) tissue_info$tissue_abbrv[tissue_info$tissue_ID==tissue])

ordering <- order(match(rownames(final_data_frac), tissue_info$tissue_abbrv), decreasing=T) #Ordering tissues by sample size
final_data_frac <- final_data_frac[ordering,]

final_data_frac <- final_data_frac[,c("Ancestry", "Sex", "Age","BMI", "MHT1D", "MHT2D",  "Atelectasis", "Atherosclerotic_Arteries",
                    "Atrophy", "Congestion", "Emphysema", "Fibrosis", "Gastritis", "Gynecomastoid", 
                    "Hashimoto", "Hyperplasia", "Pneumonia_Hist", "Saponification", "Sclerotic", "Spermatogenesis", "Steatosis")]

numbers <- rowSums(final_data_frac[,-c(1:6)]!=0) #These amounts of grey per row
to_plot <- sapply(1:nrow(final_data_frac), function(tissue) colnames(final_data_frac)[final_data_frac[tissue,]!=0]) #The ones that will be colored per tissue
#Now, create a table with the colors if we want to have the given color, and an NA if no color
to_plot <- as.data.frame(t(sapply(1:nrow(final_data_frac), function(tissue) final_data_frac[tissue,]!=0))) #The ones that will be colored per tissue
rownames(to_plot) <- rownames(final_data_frac)

to_plot$Age[to_plot$Age] <- traits_cols["Age"]
to_plot$Ancestry[to_plot$Ancestry] <- traits_cols["Ancestry"]
to_plot$Sex[to_plot$Sex] <- traits_cols["Sex"]
to_plot$BMI[to_plot$BMI] <- traits_cols["BMI"]

grays <- gray.colors(max(numbers)) 
old <- to_plot
to_plot <- old
for(i in 1:length(numbers)){ 
  print(numbers[i])
  if(to_plot[i,"MHT1D"]==T){
    to_plot[i,"MHT1D"] <- "#4C2C69"
  } 
  if(to_plot[i,"MHT2D"]==T){
    to_plot[i,"MHT2D"] <- "#805D93"
  }
  n <- sum(to_plot[i,]==T)
  if(n==0){next}
  else if(n==1){to_plot[i,][to_plot[i,]==T] <- grays[1]} else
  {to_plot[i,][to_plot[i,]==T] <- grays[1:(1+n-1)]}

}
to_plot[to_plot==FALSE] <- NA
colors <- unique(na.exclude(unlist(to_plot)))
names(colors) <- colors
#Heatmap with the colours:
suppressPackageStartupMessages(library(ComplexHeatmap))
colnames(to_plot)[5] <- "Type 1 diabetes"
colnames(to_plot)[6] <- "Type 2 diabetes"
colnames(to_plot)[8] <- "Atherosclerosis"
colnames(to_plot)[17] <- "Pneumonia"

ordering <- order(match(rownames(to_plot), tissue_info$tissue_abbrv)) #Ordering tissues by sample size, not reversed now
to_plot <- to_plot[ordering,]

#Sometimes the contribution to expression variation is 0, but the traits are included in the models anyways
to_plot$BMI <- "#CC79A7" #BMI is always included in the models
to_plot$`Type 1 diabetes` <- "#4C2C69" #T1D is not always included, but in these tissues it is.
to_plot$`Type 2 diabetes` <- "#805D93" #T2D is always included in the model
#The histology traits are always different than 0

#I will edit something in inkscape
final_colors <- sapply(colnames(to_plot), function(trait) na.omit(to_plot[,trait])[1]) #First non NA per column

sorted <- numbers[order(numbers, decreasing = T)]
sorted <- c(names(sorted[sorted>0]), "ADPSBQ") #These will be first in this order, and then the rest according to sample size

to_plot_sorted <- to_plot[sorted,]

to_plot_rest <- to_plot[rownames(to_plot)[!rownames(to_plot) %in% sorted],]

to_plot_final <- rbind(to_plot_sorted, to_plot_rest)

disease_type_main <- c(rep("white", 4), rep("#68c8cc", 2), rep("#ff8783", 15))
names(disease_type_main) <- c(rep("white", 4), rep("#68c8cc", 2), rep("#ff8783", 15))
top_anno <- HeatmapAnnotation(
  "Disease type" = disease_type_main,
  col = list("Disease type" = disease_type_main),
  border=T,
  show_legend = F,
  show_annotation_name =F,
  height = unit(5, "cm"))

#6A
png(paste0(first_dir, "03_Models/Hier_Part/Heatmap_legend.png"),
    units = "in", width = 2.8, height = 3.2, res = 200)
# pdf(paste0(first_dir, "03_Models/Hier_Part/Heatmap_legend.pdf"),
#     width = 2.8, height = 3.2)
Heatmap(as.matrix(to_plot_final),
        na_col = "white",
        col = colors,
        row_names_side = "left",
        column_names_side = "top",
        show_heatmap_legend = F,
        column_names_gp = gpar(fontsize = 6.5),
        top_annotation = top_anno,
        row_names_gp = gpar(fontsize = 6.5))
dev.off()

# png(paste0(first_dir, "03_Models/Hier_Part/Barplots_grays.png"),
#     units = "in", width = 5, height = 9, res = 200)
pdf(paste0(first_dir, "03_Models/Hier_Part/Barplots_grays.pdf"),
    width = 5, height = 9)
par(mfrow=c(1,1), oma = c(0,2,0,0))
barplot(t(final_data_frac),
        col = final_colors,
        horiz = T,
        border = NA,
        xlab = "% expression variance explained",
        las = 2,
        cex.names = 0.9,
        xaxt = 'n')
axis(1, at = axTicks(1))
dev.off()




#6C

final_data_frac_main <- final_data_frac[rev(rownames(to_plot_sorted)),]
png(paste0("~/Documents/mn4/Jose/03_Models/Hier_Part/Barplots_grays_main.png"),
    units = "in", width = 3.5, height = 4, res = 200) #Try 3.5 and 4
pdf(paste0("~/Documents/mn4/Jose/03_Models/Hier_Part/Barplots_grays_main.pdf"),
    width = 3.5, height = 4)
par(mfrow=c(1,1), oma = c(0,2,0,0))
barplot(t(final_data_frac_main),
        col = final_colors,
        horiz = T,
        border = NA,
        xlab = "% expression variance explained",
        las = 2,
        cex.names = 0.9,
        xaxt = 'n')
axis(1, at = axTicks(1))
dev.off()

