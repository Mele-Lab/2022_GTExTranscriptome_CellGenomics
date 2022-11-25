#!/usr/bin/env Rscript
#Code for Figure 6A. The input format is different from the zenodo table, for this code the input data is divided into subfolders, one per tissue, rather than using the whole zenodo table as input
rm(list=ls())

#Loading libraries
suppressPackageStartupMessages(library(ComplexHeatmap))
library(RColorBrewer)
suppressPackageStartupMessages(library(circlize))
library(optparse)

#Functions
my_pretty_num_function <- function(n){
  if(n==""){
    return(n)
  } else{
    prettyNum(n, big.mark = ",")
  }
}

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-d", "--dir"), type="character",
                     dest="dir", default="~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/Jose/",
                     help="Directory")

options=parse_args(parser)

first_dir=options$dir
# first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/"
dir <- paste0(first_dir, "03_Models/Tissues/")

main_data <- read.csv(paste0(first_dir, "00_Data/main_final.csv"))
folders <- main_data$Tissue
acronyms <- main_data$Acronym

tissue_info <- readRDS(paste0(first_dir,"00_Data/Tissue_info.rds"))[,1:2]
abbreviations <- merge(main_data, tissue_info, by.x = "Tissue", by.y = "tissue_ID", sort = F)
abbreviations <- abbreviations$tissue_abbrv[match(main_data$Tissue, abbreviations$Tissue)]

table <- c()
up <- c()
down <- c()
disease <- data.frame("Healthy"=1, "Diseased"=1)
n_samples <- c()

#Preprocessing. We iterate over subfolders, if zenodo is the input, we should iterate over the elements in the .rds file
for(i in 1:nrow(main_data)){
  path <- paste0(first_dir, "03_Models/Tissues/", folders[i])
  files <- list.files(path)[grepl("voom_limma",list.files(path))] #Models
  files <- files[grep(acronyms[i], files)]
  file <- files[grep("interactions", files, invert = T)]
  model <- readRDS(paste0(path, "/", file))[[paste0(acronyms[i],1)]]
  table <- c(table, sum(model$adj.P.Val<0.05))
  up <- c(up, sum(model$adj.P.Val<0.05 & model$logFC>0))
  down <- c(down, sum(model$adj.P.Val<0.05 & model$logFC<0))
  
  files_m <- list.files(path)[grepl("SampleMetadata",list.files(path))]
  file <- files_m[grep(acronyms[i], files_m)]
  metadata <- readRDS(paste0(path, "/", file))
  n_samples <- c(n_samples, nrow(metadata))
  to_append <- c("Healthy"=sum(metadata[[acronyms[i]]]==0), "Diseased"=sum(metadata[[acronyms[i]]]==1))
  disease <- rbind(disease, to_append) #Computing positives and negatives
}

disease <- disease[-1,]
positives <- disease$Diseased
disease_names <- paste0(acronyms, " - ", abbreviations)
disease_names <- sapply(disease_names, function(name) 
  if(grepl("MHT1D", name)){gsub("MHT1D", "Type 1 diabetes", name)
  } else if(grepl("MHT2D", name)){gsub("MHT2D", "Type 2 diabetes", name)
  } else if(grepl("Atherosclerotic_Arteries - ARTTBL", name)){"Atherosclerosis - ARTTBL"
  } else if(grepl("Pneumonia_Hist - LUNG", name)){"Pneumonia - LUNG"
  } else{ name
  })


names(positives) <- disease_names
names(n_samples) <- disease_names

disease_fr <- t(apply(disease, 1, function(x) x/sum(x)))

#Tissue Info, including colors
tissue_info <- readRDS(paste0(first_dir, "00_Data/Tissue_info.rds"))
size <- 12
tissues_cols <- sapply(folders, function(tissue) tissue_info[tissue_info$tissue_ID==tissue,3])
new_tissues <- sapply(folders, function(tissue) tissue_info[tissue_info$tissue_ID==tissue,2]) #acronyms
names(tissues_cols) <- new_tissues


# disease_type_main <- c(rep("#68c8cc", 11), rep("#ff8783", 20))
# names(disease_type_main) <- c(rep("#68c8cc", 11), rep("#ff8783", 20))
disease_type_main <- c(rep("#0fa3b1", 2), rep("#ff8783", 20))
names(disease_type_main) <- c(rep("#0fa3b1", 2), rep("#ff8783", 20))

names(table) <- disease_names


row_ha_left_main <- HeatmapAnnotation(
  tissue_names = anno_empty(border = FALSE,
                            width= -unit(2, "mm")), 
  "Disease type" = disease_type_main,
  col = list("Disease type" = disease_type_main), 
  "Samples" = anno_barplot(n_samples,  #positives if we want to show the number of positives instead of n
                                    gp = gpar(fill = tissues_cols,
                                              col = tissues_cols),
                                    border=F, width = unit(1.5, "cm")),
  "Diseased (%)" = anno_barplot(disease_fr,
                           gp = gpar(fill = c("#AAAAAA","#282828"),
                                     col = c("#AAAAAA","#282828")),
                           border=F, width = unit(1, "cm")),

  gap = unit(0.3,"cm"),
  show_legend = F,
  show_annotation_name = T,
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = size),
  which = "row")


# pdf(paste0(first_dir, "03_Models/Heatmaps/DEA/mock_main.pdf"),
#     width = 6, height = 7)
png(paste0(first_dir, "03_Models/Heatmaps/DEA/mock_main.png"), units="in",
    width = 6, height = 7, res=200)
Heatmap(as.matrix(table),
        col= colorRamp2( c(0,1,1000,1500,2500,3000),
                         brewer.pal(8, "BuPu")[c(1,2,4,5,6,7)]),
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    grid_width = unit(1, "cm"),
                                    labels_gp=gpar(fontsize=size),
                                    title_gp=gpar(fontsize=size, fontface=2)),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Number of DEGs",
        row_names_side = "left",
        column_labels="",
        row_names_gp = gpar(fontsize = size),
        left_annotation = row_ha_left_main,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(my_pretty_num_function(table[i]), x, y, gp = gpar(fontsize = size))}
)

dev.off()


#New order:
order <- c(3,9,12,18, #LUNG
           5,11,20,22, #TESTIS
           6,16,17, #THYROID
           4,13, #ARTTBL
           10,21, #LIVER
           7, #MSCLSK
           8, #ESPMCS
           15, #BREAST
           14, #STMACH
           19, 1, #PNCREAS
           2#ADPSBQ
           )


new_table <- table[order]
disease_fr_new <- disease_fr[order,]
n_samples_new <- n_samples[order]
disease_type_new <- disease_type_main[order]
tissues_cols_new <- tissues_cols[order]

row_ha_left_main <- HeatmapAnnotation(
  tissue_names = anno_empty(border = FALSE,
                            width= -unit(2, "mm")), 
  "Disease type" = disease_type_new,
  col = list("Disease type" = disease_type_new), 
  "Samples" = anno_barplot(n_samples_new,  #positives if we want to show the number of positives instead of n
                           gp = gpar(fill = tissues_cols_new,
                                     col = tissues_cols_new),
                           border=F, width = unit(1.5, "cm")),
  "Diseased (%)" = anno_barplot(disease_fr_new,
                                gp = gpar(fill = c("#AAAAAA","#282828"),
                                          col = c("#AAAAAA","#282828")),
                                border=F, width = unit(1, "cm")),
  
  gap = unit(0.3,"cm"),
  show_legend = F,
  show_annotation_name = T,
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = size),
  which = "row")

names(new_table)[names(new_table)=="Gynecomastoid - BREAST"] <- "Gynecomastia - BREAST"
pdf(paste0(first_dir, "03_Models/Heatmaps/DEA/mock_main_ordered.pdf"),
    width = 6, height = 7)
Heatmap(as.matrix(new_table),
        col= colorRamp2( c(0,1,1000,1500,2500,3000),
                         brewer.pal(8, "BuPu")[c(1,2,4,5,6,7)]),
        heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                    grid_width = unit(1, "cm"),
                                    labels_gp=gpar(fontsize=size),
                                    title_gp=gpar(fontsize=size, fontface=2)),
        na_col = "white",
        cluster_rows = F,
        cluster_columns = F,
        name = "Number of DEGs",
        row_names_side = "left",
        column_labels="",
        row_names_gp = gpar(fontsize = size),
        left_annotation = row_ha_left_main,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(my_pretty_num_function(new_table[i]), x, y, gp = gpar(fontsize = size))}
)

dev.off()

