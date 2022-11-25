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

# Functions ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Tissue ----
tissue <- args[1]
#tissue <- "BrainHypothalamus"

# Paths ----
# To save: 
# -- hier.part of PCG and lincRNA
# -- expression residuals of PCG and lincRNA
# -- Differential expression analysis (DEA) results
outpath <- args[2] #paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/03_DEA/Tissues/",tissue,"/")
outpath <- paste0(outpath, tissue, "/")
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

tissue_info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissue_info.rds")
#tissue_info <- readRDS("~/GTEx_v8/Raquel/00_Data/Tissue_Info.rds")
tissue_id <- tissue_info[tissue_info$tissue_ID==tissue,]$tissue_id

# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
gene_annotation  <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#gene_annotation  <- read.delim("~/bsc83_Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
# 26,196 genes
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# 1.2 Genes expressed per tissue ----
# 1.2.1 TPM>=0.1 in at least 20% of the tissue samples
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]  

# 1.2.2 Count >=6 in at least 20% of the tissue samples
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]  

# 1.2.3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts) 

# Exclude chrY genes in female-only tissues
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$Ensembl_ID
if(tissue %in% c("Uterus","Ovary","Vagina","BreastMammaryTissue_Female")){
    exprs_genes <- exprs_genes[!exprs_genes %in% Y_genes]
}

# 1.3 Save gene list ----
# saveRDS(exprs_genes,
#         paste0(outpath,tissue,".SelectedSamples.Expressed_genes.rds"))

# 2. Calculate gene residuals ####
print("# ---- Calculating residuals ---- #")

# 2.1 Variables ----
if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    individual_traits <- c("Age","Ancestry", "Sex","BMI")
}else{
    individual_traits <- c("Age","Ancestry", "BMI")
}
covariates <- colnames(metadata)[! colnames(metadata) %in% c("Donor", "Sample", individual_traits)]

# 2.2 Create DGEList object ----
dge <- DGEList(counts[exprs_genes,])

# 2.3 Calculate normalization factors (does not do the normalization, only computes the factors) ----
dge <- calcNormFactors(dge)

# 2.4 Voom ----
v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates

# 2.5 limma fit: Expression ~ covariates  ----
fml_args_mod <- paste(c(covariates), collapse = " + ")
mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
fit <- lmFit(v, mod)

# 2.6 Calculate expression residuals ----
exprs_residuals <- resid(fit, v)


# 3. hier.part ####
print("# ---- Running hier.part ---- #")

# -------------- #
print(Sys.time())
# -------------- #

# 3.1 Calculate explained variance ----
hier.part.results <- mclapply(rownames(exprs_residuals), function(gene)
    hier.part.mod(y=exprs_residuals[gene,], x=metadata[,individual_traits], fam = "gaussian", gof = "Rsqu"),
    mc.cores = n_cores)
names(hier.part.results) <- rownames(exprs_residuals)

# -------------- #
print(Sys.time())
# -------------- #

# 3.2 Exclude genes with negative estimates ----
print(paste0("Genes with unestimable contributions: ", sum(is.na(hier.part.results))))
if(sum(is.na(hier.part.results)) > 0){
    exprs_residuals <- exprs_residuals[-which(is.na(hier.part.results)),]
    hier.part.results <- hier.part.results[-which(is.na(hier.part.results))]
}

# 3.3 Parse results ----
rsq <- sapply(rownames(exprs_residuals), function(gene) 
    sum(hier.part.results[[gene]]$IJ[,1])) 
names(rsq) <- rownames(exprs_residuals)
rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part.results), function(gene) 
                        as.numeric(unlist(hier.part.results[[gene]]$I.perc))))
rownames(rel_perc) <- rownames(exprs_residuals)
colnames(rel_perc) <- individual_traits
abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part.results), function(gene) 
                         hier.part.results[[gene]]$IJ[,1])
)
rownames(abs_perc) <- rownames(exprs_residuals)
colnames(abs_perc) <- individual_traits
hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))


# 3.4. Save results ----
saveRDS(exprs_residuals,paste0(outpath, tissue,".exprs_residuals.rds"))
saveRDS(hier_data,paste0(outpath,tissue, ".hier_part.rds"))


# 4. Differential expression analysis ####
print("# ---- Running differential expression analysis ---- #")

# 4.1 limma fit : expression ~ covariates + traits ----
my_data <- list()
fml_args_mod <- paste(c(covariates, individual_traits), collapse = " + ")
mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# 4.2 Limma fit ----
fit <- lmFit(v, mod)

# Add objects to data
my_data[["dge"]] <- dge
my_data[["v"]] <- v
my_data[["fit"]] <- fit

# 4.3 Limma test ----
if(tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    dea_res <- list()
    dea_res[[1]] <- NA
    dea_res <- c(dea_res, 
                 lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata) ))
    names(dea_res) <- c("Sex",individual_traits)
    dea_res <- dea_res[c("Age","Ancestry", "Sex","BMI")]
}else{
    dea_res <- lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) )
    names(dea_res) <- individual_traits
}

# 4.3 Limma test using sum 2 zero ----
if(!tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    metadata$Sex <- factor(metadata$Sex, levels=c(2,1)) #1 will be the reference (male),
    contrasts(metadata$Sex) <- contr.sum(2)
}
metadata$Ancestry <- factor(metadata$Ancestry, levels=c("AFR","EUR")) #Reference is European, sobreescribiendo el relevel
contrasts(metadata$Ancestry) <- contr.sum(2)
contrasts(metadata$HardyScale) <- contr.sum(length(levels(metadata$HardyScale)))

if(tissue %in% c("Vagina","Uterus","Ovary","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    dea_res.sum2zero <- list()
    dea_res.sum2zero[[1]] <- NA
    dea_res.sum2zero <- c(dea_res.sum2zero, 
                 lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata) ))
    names(dea_res.sum2zero) <- c("Sex",individual_traits)
    dea_res.sum2zero <- dea_res.sum2zero[c("Age","Ancestry", "Sex","BMI")]
}else{
    dea_res.sum2zero <- lapply(individual_traits, function(phenotype) limma_lm(fit,phenotype,metadata ) )
    names(dea_res.sum2zero) <- individual_traits
}

#sapply(individual_traits, function(trait) sum(dea_res[[trait]]$adj.P.Val < 0.05))
#sapply(individual_traits, function(trait) sum(dea_res.sum2zero[[trait]]$adj.P.Val < 0.05))

# 5. Computing avrg TPM and exprs var ####
print("# ---- Calculating avrg TPM and var ---- #")

# 5.1 Compute average TPM expression and variance for each event
avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
median_TPM <- apply(tpm, 1, function(x) median(log2(x+1)))
var_TPM <- apply(tpm, 1, function(x) var(x))

# Add AvgExprs & ExprsVar and order data.frame
for(trait in individual_traits){
    print(trait)
    if(trait == "Sex"){
        if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
            dea_res[[trait]] <- NA
        }else{
            if(length(dea_res[[trait]])>1){ # KidneyCortex: Ancestry -< NA
                # Add average expression TPM
                dea_res[[trait]]$AvgTPM <- sapply(rownames(dea_res[[trait]]), function(gene) avrg_TPM[gene])
                dea_res[[trait]]$MedianTPM <- sapply(rownames(dea_res[[trait]]), function(gene) median_TPM[gene])
                dea_res[[trait]]$VarTPM <- sapply(rownames(dea_res[[trait]]), function(gene) var_TPM[gene])
                gene_names <- sapply(rownames(dea_res[[trait]]), function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name"])
                names(gene_names) <- NULL
                dea_res[[trait]][["gene.name"]] <- gene_names
            }
        }
    }else{
        if(length(dea_res[[trait]])>1){
            # Add average expression TPM and variance
            dea_res[[trait]]$AvgTPM <- sapply(rownames(dea_res[[trait]]), function(gene) avrg_TPM[gene])
            dea_res[[trait]]$MedianTPM <- sapply(rownames(dea_res[[trait]]), function(gene) median_TPM[gene])
            dea_res[[trait]]$VarTPM <- sapply(rownames(dea_res[[trait]]), function(gene) var_TPM[gene])
            gene_names <- sapply(rownames(dea_res[[trait]]), function(gene) gene_annotation[gene_annotation$ensembl.id==gene, "gene.name"])
            names(gene_names) <- NULL
            dea_res[[trait]][["gene.name"]] <- gene_names
        }
    }
} 


# 5.2 Save table with results
saveRDS(dea_res,
        paste0(outpath,tissue,".voom_limma.covariates_and_traits.results.rds"))
saveRDS(dea_res.sum2zero,
        paste0(outpath,tissue,".voom_limma.covariates_and_traits.results.sum2zero_contrast.rds"))
saveRDS(my_data,
        paste0(outpath,tissue,".voom_limma.covariates_and_traits.data.rds"))

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

