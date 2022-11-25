#Code for checking how a gene's tpm change with respect to another variable (e.g., Age, BMI, Sex)
#This code uses the TPMs available from the GTEx portal

#!/usr/bin/env Rscript

first_dir <- "~/Documents/mn4/projects/bsc83/Projects/GTEx_v8/"
# first_dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"

outdir <- paste0(first_dir, "Jose/03_Models/Tissues/")

gene_annotation <- read.delim(paste0(first_dir,"Jose/00_Data/gencode.v26.GRCh38.genes.bed"), header=F)[,c(6,7)]
colnames(gene_annotation) <- c("gene", "symbol")

traits_cols <- c("Age" = "#56B4E9", "Ancestry" = "#E69F00", 
                 "Sex" = "#009E73", "BMI" = "#CC79A7", "Disease" = "#696969")
my_traits <- names(traits_cols)

# Libraries ----
library(ggplot2)
library(gridExtra)
library(gtools)
library(ggpubr)

# Functiions ----
get.gene.data <- function(acronym, tissue, gene, combo, g = 2){
  #  Expression data ----
  #tpm should be changed for the GTEx portal file
  tpm <- readRDS(paste0(first_dir, "/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))
  colnames(tpm) <- sapply(colnames(tpm), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  #res are the expression residuals of linear models with technical covariates, but we will only use TPMs, so the residuals are not needed
  res <- readRDS(paste0(first_dir, "/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/", tissue, "/", tissue, ".exprs_residuals.rds"))
  colnames(res) <- sapply(colnames(res), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )

  # metadata <- readRDS(paste0(outdir, disease, "/Tissues/", tissue, "/", tissue, ".SampleMetadata.", disease,".rds"))
  files <- list.files(paste0(outdir, tissue, "/"))
  metadata_file <- files[grep("SampleMetadata", files)]
  metadata <- readRDS(paste0(outdir, tissue, "/", metadata_file))
  metadata[[acronym]] <- as.factor(metadata[[acronym]])
  metadata$Sample <- sapply(metadata$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  
  #Subset with the gene data
  gene.name <- gene_annotation[gene_annotation$gene==gene, "symbol"]
  df <- cbind.data.frame(tpm[gene,])
  colnames(df) <- "TPM"
  df$Residuals <- res[gene,]
  
  #Subset for what we have metadata
  df <- df[rownames(df) %in% metadata$Sample,]
  
  df$Age_int <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Age"])
  df$BMI_int <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"BMI"])
  metadata$Age_Class <- cut(metadata$Age, breaks = c(20, 30, 40, 50, 60, 70), include.lowest=T, right=F) #Include_lowest allows 20 to be included
  levels(metadata$Age_Class) <- c("20-29", "30-39", "40-49", "50-59", "60-70")
  df$Age <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Age_Class"])
  df$Age <- factor(df$Age, levels = c("20-29",
                                      "30-39",
                                      "40-49",
                                      "50-59",
                                      "60-70"),
                   order = T)

  df$Age2 <- ifelse(df$Age_int < 45, "[20-45)", "[45-70]")
  df$Age2 <- factor(df$Age2,
                    levels = c("[20-45)", "[45-70]"),
                    order = T)
  df$Ancestry <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  metadata$BMI_Class <- cut(metadata$BMI, breaks = c(15, 25, 30, 36), include.lowest=T, right=F) #Include_lowest allows 20 to be included
  
  df$BMI <- sapply(rownames(df), function(i) metadata[metadata$Sample==i,"BMI_Class"])
  levels(df$BMI) <- c("Normal", "Overweight", "Obese")
  df$BMI2 <- ifelse(df$BMI_in < 25, "<25", ">=25")
  df$BMI2 <- factor(df$BMI2,
                    levels = c("<25", ">=25"),
                    order = T)
  df$Tissue <- rep(tissue, nrow(df))
  df$Gene <- rep(gene.name, nrow(df))
  
  df$Disease <- sapply(rownames(df), function(i) metadata[metadata$Sample==i, acronym])
  if(combo == "Age:Sex"){
    if(g==1){
      df$x_dummy <- paste0(df$Age, "_",df$Sex)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("20-29",
                                                    "30-39",
                                                    "40-49",
                                                    "50-59",
                                                    "60-70"),  function(i) paste0(i, "_", c("Male", "Female")))),
                           order = T)    
    }else{
      df$x_dummy <- paste0(df$Age2, "_",df$Sex)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("[20-39]", "[40-70]"),
                                                  function(i) paste0(i, "_", c("Male", "Female")))),
                           order = T)  
    }
  }else if(combo == "Age:Ancestry"){
    if(g==1){
      df$x_dummy <- paste0(df$Age, "_",df$Ancestry)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("20-29",
                                                    "30-39",
                                                    "40-49",
                                                    "50-59",
                                                    "60-70"),  function(i) paste0(i, "_", c("EA", "AA")))),
                           order = T)    
    }else{
      df$x_dummy <- paste0(df$Age2, "_",df$Ancestry)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("[20-39]", "[40-70]"),
                                                  function(i) paste0(i, "_", c("EA", "AA")))),
                           order = T)  
    }
  }else if(combo == "Ancestry:BMI"){
    if(g==1){
      df$x_dummy <- paste0(df$BMI, "_",df$Ancestry)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("Normal", "Overweight", "Obese"),  function(i) paste0(i, "_", c("EA", "AA")))),
                           order = T)    
    }else{
      df$x_dummy <- paste0(df$BMI, "_",df$Ancestry)
      df$x_dummy <- factor(df$x_dummy,
                           levels = unlist(lapply(c("<25", ">=25"),  function(i) paste0(i, "_", c("EA", "AA")))),
                           order = T)  
    }
  } else if(combo == "Disease:Age"){
    df$x_dummy <- paste0(df$Age2, "_", df$Disease)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("[20-45)", "[45-70]"),
                                                function(i) paste0(i, "_", c("0", "1")))),
                         order = T)  
  }else if(combo == "Disease:Sex"){
    df$x_dummy <- paste0(df$Sex, "_", df$Disease)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("Male", "Female"),
                                                function(i) paste0(i, "_", c("0", "1")))),
                         order = T)  
  }else if(combo == "Disease:Ancestry"){
    df$x_dummy <- paste0(df$Ancestry, "_", df$Disease)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("EA", "AA"),
                                                function(i) paste0(i, "_", c("0", "1")))),
                         order = T)  
  }else if(combo == "Disease:BMI"){
    df$x_dummy <- paste0(df$BMI2, "_", df$Disease)
    df$x_dummy <- factor(df$x_dummy,
                         levels = unlist(lapply(c("<25", ">=25"),
                                                function(i) paste0(i, "_", c("0", "1")))),
                         order = T)  
  }
  return(df)
}

get_box_stats <- function(y, upper_limit = max(log2(data$TPM)) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "N =", length(y), "\n"
    )
  ))
}

my.box_plot <- function(data, gene, variable_col, x_variable){
  cols <- c("white",  "white")
  facet_color <- "#805D93"
  if(x_variable == "Age_int"){ x_facet <- "Age"
  }else if(x_variable == "Age2"){ 
    x_facet <- "Age2"
    facet_color <- traits_cols["Age"]
  }else if(x_variable == "BMI_int"){ 
    x_facet <- "BMI"
    facet_color <- traits_cols["BMI"]
  }else if(x_variable == "BMI"){ x_facet <- "BMI2"}
  else{x_facet <- x_variable}
  gene.name <- gene_annotation[gene_annotation$ensembl.id == gene, "gene.name"]
  p <- ggplot(data = data,
              aes(x = x_dummy,
                  y = log2(TPM+1))
  ) +
    geom_violin(aes(fill = eval(parse(text=variable_col))),
                col = "black") +
    geom_boxplot(col = "black",
                 outlier.shape = NA,
                 notch = T,
                 width = 0.25) +
    geom_jitter(col = "black", 
                alpha = 0.1,
                size = 0.8) +
    theme_minimal() +
    xlab("") +
    ylab("log2(TPM)") +
    scale_fill_manual(values = cols) +
    labs(title=paste0(gene.name,
                      " (", tissue, ")")) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5,
                                    size = 10)) +
    stat_summary(aes(x = x_dummy, y = log2(TPM)),
                 fun.data = get_box_stats, geom = "text",
                 hjust = 0.5, vjust = 0.9, size = 2) +
    facet_grid(.~ eval(parse(text = x_facet)), drop = T, scales = "free") +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 0.1),
          strip.background = element_rect(fill=facet_color)) + ##B3B3B3
    theme(legend.position = "none")
  return(p)
}

acronym <- "MHT2D"
tissue <- "NerveTibial"
variables <- "Disease:Age"
gene <- "ENSG00000175445.14" #LPL

p1 <- my.box_plot(data, gene, "Age2", "Disease")
p1

hier_file <- list.files(paste0(outdir, tissue), pattern = ".hier_part", full.names = T)
hier_file <- hier_file[grep(hier_file, pattern = "spliced", invert = T)]
hier.part.exprs <- readRDS(hier_file)

demographic_trati <- "Age"
clinical_trait <- "Type 2 diabetes"

d <- data.frame("variable" = c(clinical_trait, demographic_trati),
                "value"  = as.numeric(hier.part.exprs[gene, c(paste0(acronym, "_abs"), paste0(demographic_trati, "_abs"))]))
variables_col <- c(traits_cols["Disease"], traits_cols[demographic_trati])
variables_col[1] <- "#805D93"
names(variables_col) <- c(clinical_trait, demographic_trati)
d$variable <- factor(d$variable, levels = rev(c(clinical_trait, demographic_trati)), order = T)
p3 <- ggplot(d, aes(x=1, y=value, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = variables_col) + 
  coord_flip() +
  xlab("") + ylab("Expression variation explained (%)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.x =  element_line(size = 0.5),
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 
p3

png(paste0(first_dir, "Jose/03_Models/Overlaps/Final_", gene_annotation[gene_annotation$gene == gene, "symbol"], "_", acronym, "_", tissue, "_hier_part.png"),
    units = "in", width = 4, height = 3.5, res=150)
pdf(paste0(first_dir, "Jose/03_Models/Overlaps/Final_", gene_annotation[gene_annotation$gene == gene, "symbol"], "_", acronym, "_", tissue, "_hier_part.pdf"),
    width = 4, height = 3.5)
ggarrange(p1, p3, nrow=2, heights = c(4,1))#, labels = paste0("LPL (", tissue, ")"), hjust = -1)
dev.off()
