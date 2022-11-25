
#!/usr/bin/env Rscript

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()
# ---------------------- #

# Libraries ####
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(hier.part))
suppressMessages(library(parallel))
library(gtools)

# Functions ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Tissue ----
tissue <- args[1]
#tissue <- "ArteryAorta"

# Paths ----
# To save: 
# -- hier.part of PCG and lincRNA
# -- expression residuals of PCG and lincRNA
# -- Differential expression analysis (DEA) results
inpath <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/"
#inpath <- "~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/"
inpath <- paste0(inpath, tissue, "/")
outpath <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/02.Interactions/Tissues/"
outpath <- paste0(outpath, tissue, "/")
#outpath <- "~/GTEx_v8/Raquel/Draft/01.DiffExprs/10.Interactions/Tissues/"
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

# Number of CPU (cores) to parallelize mclapply
n_cores <- args[3]
#n_cores <- 4

# --------------- #
print("#############################")
print(paste0("Tissue: ", tissue))
print(paste0("Outpath: ", outpath))
print(paste0("n: ", n_cores))
print("#############################")
print("")
# --------------- #

#### Code ####

# 1. Select PCG and lincRNA expressed in tissue ####

# 1.1 Read in data ----
counts <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
#counts <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))
#tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))
metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
#metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")
degs <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/",tissue,"/",tissue,".voom_limma.covariates_and_traits.results.rds"))
#degs <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/",tissue, "/", tissue, ".voom_limma.covariates_and_traits.results.rds"))

# Age --
pseudo_categorize_age <- function(age){
  if(age < 45){
    return("Young")
  }else{
    return("Old")
  }
}
metadata$Age_Class <-  sapply(metadata$Age, function(age) pseudo_categorize_age(age))
metadata$Age_Class <- factor(metadata$Age_Class, levels = c("Young", "Old"), order = T)

# BMI --
pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}
metadata$BMI_Class <-  sapply(metadata$BMI, function(bmi) pseudo_categorize_bmi(bmi))
metadata$BMI_Class <- factor(metadata$BMI_Class, levels = c("Normal", "Overweight", "Obese"), order = T)

# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
gene_annotation  <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#gene_annotation  <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
# 26,196 genes
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# 1.2 Genes expressed per tissue ----
# 1.2.1 TPM>=1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  

# 1.2.2 Count >=10 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  

# 1.2.3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

# Exclude chrY genes in female-only tissues
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$Ensembl_ID
if(tissue %in% c("Uterus","Ovary","Vagina","BreastMammaryTissue_Female")){
  exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
}

# voom object ----
v <- readRDS(paste0(inpath,tissue,".voom_limma.covariates_and_traits.data.rds"))[["v"]]

# 1. Differential expression analysis with interaction ####
print("# ---- Running differential expression analysis ---- #")

# 1.1 Variables ----
if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  individual_traits <- c("Age","Ancestry", "Sex","BMI")
}else{
  individual_traits <- c("Age","Ancestry", "BMI")
}

individual_traits.interaction <- individual_traits[sapply(individual_traits, function(trait) sum(degs[[trait]]$adj.P.Val < 0.05))>0]
pw.traits <- combinations(length(individual_traits.interaction), r = 2, individual_traits.interaction)
interaction.terms <- sapply(1:nrow(pw.traits), function(row) paste(c(pw.traits[row,1], 
                                                                     pw.traits[row,2]),
                                                                   collapse = ":"))
if("Age:BMI" %in% interaction.terms){
  interaction.terms <- interaction.terms[-which(interaction.terms == "Age:BMI")]
}
if("BMI:Sex" %in% interaction.terms){
  interaction.terms[which(interaction.terms == "BMI:Sex")] <- "Sex:BMI"
}

keep.term.fun <- function(term){
  if(term == "Age:Ancestry"){
    ifelse(min(table(metadata$Age_Class, metadata$Ancestry)) < 20, 0, 1)
  # }else if(term == "Age:BMI"){
  #   ifelse(min(table(metadata$Age_Class, metadata$BMI_Class)) < 20, 0, 1)
  }else if(term == "Age:Sex"){
    ifelse(min(table(metadata$Age_Class, metadata$Sex)) < 20, 0, 1)
  }else if(term == "Ancestry:BMI"){
    ifelse(min(table(metadata$Ancestry, metadata$BMI_Class)) < 20, 0, 1)
  }else if(term == "Ancestry:Sex"){
    ifelse(min(table(metadata$Ancestry, metadata$Sex)) < 20, 0, 1)
  }else if(term == "Sex:BMI"){
    ifelse(min(table(metadata$BMI_Class, metadata$Sex)) < 20, 0, 1)
  }else{
    return(NA)
  }
}

interaction.terms <- names(sapply(interaction.terms, function(term) keep.term.fun(term))[sapply(interaction.terms, function(term) keep.term.fun(term))==1])
metadata <- metadata[, -which(colnames(metadata) %in% c("Age_Class", "BMI_Class"))]
covariates <- colnames(metadata)[! colnames(metadata) %in% c("Donor", "Sample", individual_traits)]

# Limma test using sum2zero for individual traits ----
if(!tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate")){
  metadata$Sex <- factor(metadata$Sex, levels=c(2,1)) #1 will be the reference (male),
  contrasts(metadata$Sex) <- contr.sum(2)
}
metadata$Ancestry <- factor(metadata$Ancestry, levels=c("AFR","EUR")) #Reference is European, sobreescribiendo el relevel
contrasts(metadata$Ancestry) <- contr.sum(2)
contrasts(metadata$HardyScale) <- contr.sum(length(levels(metadata$HardyScale)))
if("Cohort" %in% colnames(metadata)){
  contrasts(metadata$Cohort) <- contr.sum(length(levels(metadata$Cohort)))
}
if("NucAcIsoBatch" %in% colnames(metadata)){
  contrasts(metadata$NucAcIsoBatch) <- contr.sum(length(levels(metadata$NucAcIsoBatch)))
}


# 1.2 limma fit : expression ~ covariates + traits ----
fml_args_mod <- paste(c(covariates, individual_traits, interaction.terms), collapse = " + ")
mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# 1.3 Limma fit ----
fit <- lmFit(v, mod)
colnames(fit$design)

# 1.4 Limma test ----
dea_res <- lapply(interaction.terms, function(phenotype) limma_lm.interactions(fit,phenotype) )
names(dea_res) <- interaction.terms

# 1.5. Save table with results
saveRDS(dea_res,
        paste0(outpath,tissue,".voom_limma.interaction_terms.rds"))

# # 1.6 Limma test using sum 2 zero for individual traits ----
# if(!tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
#   metadata$Sex <- factor(metadata$Sex, levels=c(2,1)) #1 will be the reference (male),
#   contrasts(metadata$Sex) <- contr.sum(2)
# }
# metadata$Ancestry <- factor(metadata$Ancestry, levels=c("AFR","EUR")) #Reference is European, sobreescribiendo el relevel
# contrasts(metadata$Ancestry) <- contr.sum(2)
# contrasts(metadata$HardyScale) <- contr.sum(length(levels(metadata$HardyScale)))

if(tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate")){
  dea_res.sum2zero <- list()
  dea_res.sum2zero[[1]] <- NA
  dea_res.sum2zero <- c(dea_res.sum2zero, 
                        lapply(individual_traits, function(phenotype) limma_lm.sum2zero(fit,phenotype,metadata) ))
  names(dea_res.sum2zero) <- c("Sex",individual_traits)
  dea_res.sum2zero <- dea_res.sum2zero[c("Age","Ancestry", "Sex","BMI")]
}else{
  dea_res.sum2zero <- lapply(individual_traits, function(phenotype) limma_lm.sum2zero(fit,phenotype,metadata ) )
  names(dea_res.sum2zero) <- individual_traits
}

saveRDS(dea_res.sum2zero,
        paste0(outpath,tissue,".voom_limma.indivual_traits.rds"))

# # 2. Computing avrg TPM and exprs var ####
# print("# ---- Calculating avrg TPM and var ---- #")
# 
# # 5.1 Compute average TPM expression and variance for each event
# avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
# median_TPM <- apply(tpm, 1, function(x) median(log2(x+1)))
# var_TPM <- apply(tpm, 1, function(x) var(x))
# 
# # Add AvgExprs & ExprsVar and order data.frame
# for(trait in interaction.terms){
#   print(trait)
#   dea_res[[trait]]$AvgTPM <- sapply(rownames(dea_res[[trait]]), function(gene) avrg_TPM[gene])
#   dea_res[[trait]]$MedianTPM <- sapply(rownames(dea_res[[trait]]), function(gene) median_TPM[gene])
#   dea_res[[trait]]$VarTPM <- sapply(rownames(dea_res[[trait]]), function(gene) var_TPM[gene])
#   gene_names <- sapply(rownames(dea_res[[trait]]), function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name"])
#   names(gene_names) <- NULL
#   dea_res[[trait]][["gene.name"]] <- gene_names
# } 
# 


# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Session info ---- #")
sessionInfo()
print("# ---------------------- #\n")
print("")
# ---------------------- #

