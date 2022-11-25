#!/usr/bin/env Rscript

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()
# ---------------------- #

# Libraries ####
library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(caret)
library(tidyr)
library(broom)
library(multcomp)
library(hier.part)

# Functions ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/Expression_R_functions.R")
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/Expression_R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Tissue
tissue <- args[1]
#tissue <- "AdiposeSubcutaneous"

# Paths 
# To save: 
# ---- ancestry DEG classified as cis-driven and not cis-driven
inpath <- args[2]
#inpath <- "~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/"
inpath <- paste0(inpath, "/", tissue, "/")

#outpath <- "~/GTEx_v8/Raquel/Draft/01.DiffExprs/04.Classify_AncestryDEG/Tissues/"
outpath <- args[3] 
outpath <- paste0(outpath, "/", tissue, "/")
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

print("#### ----- ####")
print(paste0("Tissue: ", tissue))
print(paste0("Inpath: ", inpath))
print(paste0("Outpath: ", outpath))
print("#### ----- ####")

# Number of CPU (cores) to parallelize mclapply
#n_cores <- args[3]

# Metadata ----
metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
#metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")

# Covariates ----
if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis", "BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  individual_traits <- c("Ancestry","Age", "Sex","BMI")
  contrast_names <- c('AncestryEUR','Age','Sex2','BMI')
}else{
  individual_traits <- c("Ancestry","Age", "BMI")
  contrast_names <- c('AncestryEUR','Age','BMI')
}
covariates <- colnames(metadata)[! colnames(metadata) %in% c("Donor", "Sample", individual_traits)]

#### Code ####
# Model gene as:
# exprs_residual ~ ieVariants_in_eGene + traits
# The genotypes of the sVariants are modelled as a factor, where: 0|0 -> 0, 0|1 -> 1; 1|1 ->2

### Prepare data ####
print("# ---- Reading in data ---- #")

# Tissue info ----
tissue_info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissue_info.rds")
#tissue_info <- readRDS("~/GTEx_v8/Raquel/Draft/Data/Tissue_info.rds")
tissue_id <- tissue_info[tissue_info$tissue_ID==tissue,]$tissue_id

# eGenes ----
inpath_eqtls <- "/gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/"
#inpath_eqtls <- "~/GTEx_v8_data/cisQTLs/"
egene_data <- as.data.frame(data.table::fread(paste0(inpath_eqtls,"eGenes/GTEx_Analysis_v8_eQTL/",tissue_id,".v8.egenes.txt.gz")))

# Select eGenes (q.val <= 0.05) --
eGenes <- egene_data[egene_data$qval <= 0.05,"gene_id"]

# eGenes:ieQTLs ----
eqtls_dir <- "/gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_eQTL_independent/"
#eqtls_dir <- "~/GTEx_v8_data/cisQTLs/GTEx_Analysis_v8_eQTL_independent/"
eqtls.tissue <-  as.data.frame(data.table::fread(paste0(eqtls_dir, tissue_id,".v8.independent_eqtls.txt.gz")))

# sGenes with at leats 1 independent sQTL ----
eGenes.ieqtl <- unique(eqtls.tissue$gene_id)

# Differential expression analysis: dataframe with results ----
dea_res <- readRDS(paste0(inpath, tissue,".voom_limma.covariates_and_traits.results.rds"))[["Ancestry"]]

# Save summary data ----
d <- list()
d[["eGenes"]] <- length(eGenes)
d[["eGenes:ieQTL"]] <- length(eGenes.ieqtl)
d[["Ancestry:DEG"]] <- sum(dea_res$adj.P.Val < 0.05)
d[["Ancestry:DEG:eGenes"]] <- sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% eGenes)
d[["Ancestry:DEG:eGenes-ieQTL"]] <- sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% eGenes.ieqtl)

# Print to std.out summary data ----
print(paste0("No. of eGenes: ", length(eGenes)))
print(paste0("No. of eGenes-ieQTLs:  ", length(eGenes.ieqtl)))
print(paste0("No. of ancestry DEG:  ", sum(dea_res$adj.P.Val < 0.05)))
print(paste0("No. of ancestry DEG:eGenes:  ", sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% eGenes) ))
print(paste0("No. of ancestry DEG:eGenes:  ", sum(rownames(dea_res[dea_res$adj.P.Val < 0.05,])  %in% eGenes.ieqtl) ))

# Goal: select differentially expressed genes that are eGenes with at least 1 genotyped independent eQTL (eVariant) with MAF > 0.01  ####
# Select differentially expressed genes ----
dea_res <- dea_res[dea_res$adj.P.Val < 0.05,] # significant genes
dea_res <- dea_res[rownames(dea_res) %in% eGenes.ieqtl,] # i-eGenes

# Subset eQTLs in eGenes differentially expressed with at least 1 independent sQTL ----
eqtls.tissue <- eqtls.tissue[eqtls.tissue$gene_id  %in% rownames(dea_res),]

# eVariants in eGenes differentially expressed with at leats 1 independent eQTL ----
eVariants <- unique(eqtls.tissue$variant_id) # ieQTLs
print(paste0("No. of eVariants (i-eQTL) in ancestry DEG eGenes: ", length(eVariants)))

# Genotyped eVariants ----
# Read the genotypes for the sVariants which are also independent cis-eQtls 
# Job run in cluster
#vcftools --gzvcf ${vcf_file} --snps <(zcat /gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/${tissue_id}.v8.independent_sqtls.txt.gz  | cut -f7 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.txt.gz
gt_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Fst/ieQTL_GT/"
#gt_path <- "~/GTEx_v8/Raquel/Draft/Data/Fst/ieQTL_GT/"
if(!tissue %in% c("BreastMammaryTissue_Female", "BreastMammaryTissue_Male")){
  variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, tissue, ".ieQTLs.Donor_genotypes.txt.gz")))  
}else{
  variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, "BreastMammaryTissue", ".ieQTLs.Donor_genotypes.txt.gz"))) 
}

# Set proper colnames --
colnames.variants <- c('CHROM','POS','variant_id','REF','ALT')
colnames.individuals <- str_split_fixed(as.character(variants_genotype[1,-c(1:5)]),'=',2)[,1]
colnames(variants_genotype) <- c(colnames.variants,colnames.individuals)

# Retain selected tissue samples --
variants_genotype <- variants_genotype[,c(colnames.variants, metadata$Donor)]
# Donors in same order as in metadata
if(!identical(metadata$Donor, colnames(variants_genotype)[6:ncol(variants_genotype)])){
  print("Donors not in the same order in metadata and genotype tables")
  quit()
}

if(length(eVariants) > 1){
  # Subet genotype file and only keep eVariants in eGenes differentially expressed --
  variants_genotype <- variants_genotype[variants_genotype$variant_id %in% eVariants,] 
}

# Customize the genotypes files --
vg.tissue_mod <- cbind(variants_genotype[1:5], apply(variants_genotype[6:ncol(variants_genotype)],2, FUN=function(x){str_split_fixed(x,'=',2)[,2]}))
vg.tissue_wide <- cbind(vg.tissue_mod[1:5], apply(vg.tissue_mod[6:ncol(vg.tissue_mod)],2, FUN=function(x){ifelse(x=='0|0','0',
                                                                                                                 ifelse(x=='0|1','1',
                                                                                                                        ifelse(x=='1|1','2',NA)))}))
vg.tissue_long <- reshape2::melt(vg.tissue_wide,
                                 id.vars = colnames.variants,
                                 variable.name = 'Individual_ID',
                                 value.name = 'genotype')
if(length(eVariants) == 1){
  # Subet genotype file and only keep eVariants in eGenes differentially expressed --
  variants_genotype <- variants_genotype[variants_genotype$variant_id %in% eVariants,]
  vg.tissue_mod <- vg.tissue_mod[vg.tissue_mod$variant_id %in% eVariants,]
  vg.tissue_wide <- vg.tissue_wide[vg.tissue_wide$variant_id %in% eVariants,]
  vg.tissue_long <- vg.tissue_long[vg.tissue_long$variant_id %in% eVariants,]
}

# Subset tissue samples ----
vg.tissue_long <- vg.tissue_long[vg.tissue_long$Individual_ID %in% metadata$Donor,]

# eVariants in VCF -> eVariants in eGenes with at leats 1 independent eQTL with MAF > 0.001 ----
eVariants_lost <- eVariants[!eVariants %in% unique(vg.tissue_long$variant_id)] # MAF < 0.001?
eVariants <- eVariants[eVariants %in% unique(vg.tissue_long$variant_id)]
print(paste0("No. of eVariants (i-eQTL) with MAF > 0.01 in ancestry DEG eGenes: ", length(eVariants)))
print(paste0("Minimum MAF value of genotyped eVariants: ",  min(eqtls.tissue[eqtls.tissue$variant_id %in% eVariants, "maf"])))

# Add info to summary data ----
d[["ieQTL:MAF_b_01"]] <- length(eVariants_lost)
d[["ieQTL:MAF_01"]] <- length(eVariants)
d[["minMAF"]] <-  min(eqtls.tissue[eqtls.tissue$variant_id %in% eVariants, "maf"])

# Subset eQTLs in eGenes differentially spliced with at least 1 genotyped independent eQTL with MAF > 0.001 ----
eqtls.tissue <- eqtls.tissue[eqtls.tissue$variant_id %in% eVariants,]
ieGenes.de.maf01 <- unique(eqtls.tissue$gene_id)

# Subset genes with ieQTL with MAF 001 (in vcf file) ----
dea_res <- dea_res[rownames(dea_res) %in% ieGenes.de.maf01,]
d[["Ancestry:DEG:eGenes-ieQTL:MAF01"]] <- length(ieGenes.de.maf01)


# Create 'dictionary' (named list)  ----
gene_variants.list <- sapply(ieGenes.de.maf01, function(x){
  variant <- unique(eqtls.tissue[eqtls.tissue$gene_id==x,]$variant_id)
  return(variant)
}, simplify=F)

############### Linear model for each event #########################

## To build the model:
# 1. For each gene, recover its eVariants (ieQTLS with MAF 001) ----
# 2. Compare 2 models:
# residuals ~ Age + Sex + BMI + cis-effects (ieQTL(s))
# residuals ~ Age + Sex + BMI + cis-effects (ieQTL(s)) + Ancestry 
# Is Ancestry adding someting that the sQTLs did not capture, a not cis-driven effect?


# Genes  ----
genes <- rownames(dea_res)

# Modify to match ----
colnames(metadata)[1] <- "Individual_ID"
colnames(metadata)[2] <- "Sample_ID"
metadata <- metadata[,c("Individual_ID", "Sample_ID",  individual_traits)]
it <- "Ancestry"
if(!tissue %in% c("Uterus","Ovary","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  its <- c("Age","Sex","BMI","Ancestry")    
}else{
  its <- c("Age","BMI","Ancestry")
}

# Recover the expression residuals associated with DEGs ----
exprs_residuals <- readRDS(paste0(inpath, tissue,".exprs_residuals.rds"))

# Susbet residuals of genes to be modelled ----
if(length(rownames(dea_res))==1){
  exprs_residuals <- as.data.frame(t(as.matrix(exprs_residuals[rownames(dea_res),])))
  rownames(exprs_residuals) <- rownames(dea_res)
}else{
  exprs_residuals <- exprs_residuals[rownames(dea_res),]    
}

#identical(metadata$Sample, colnames(exprs_residuals))
if(!identical(metadata$Sample_ID, colnames(exprs_residuals))){
  print("The samples in the metadata and the residual data do not coincide")
  q()
}
if(!identical(metadata$Individual_ID, as.character(unique(vg.tissue_long$Individual_ID)))){
  print("The donors in the metadata and the genotyped data do not coincide")
  q()
}

# Is the Ancestry effect on the eGene cis-driven or not cis-driven ----
# Is the event cis-driven or not cis-driven ----
print('*****************************************************')
print(paste0(tissue, ' has a total of ', length(genes), ' DEGs with cis-independent eQTLs with MAF > 0.01'))
print('Applying glm() function to each of those events')
print('*****************************************************')
cat('\n')
cat('\n')

# lm per event ----
g.lm_models <- sapply(genes, function(g) lm.cis_models(g), simplify = F)
names(g.lm_models) <-  genes
# genes that  cannot modelled ----
d[["Ancestry:DEG:NotModelled"]] <- sum(is.na(g.lm_models))
# genes modelled ----
g.lm_models <- g.lm_models[!is.na(g.lm_models)] # SNP with no variance or no SNP not correlated
d[["Ancestry:DEG:Modelled"]] <- length(g.lm_models)

# Compare modelA vs modelB ----
deg <- names(g.lm_models)
# Anova
anova.P.Value <- sapply(deg, function(g) anova(g.lm_models[[g]][["lm.modelA"]],g.lm_models[[g]][["lm.modelB"]], test="F")[2,6])
# Multiple testing correction
anova.adj.P.Val <- p.adjust(anova.P.Value, method = "BH")

# Parse results data ----
results.df <- cbind.data.frame(deg, anova.adj.P.Val)
results.df$Class <- ifelse(results.df$anova.adj.P.Val < 0.05, "Not_cis-driven","Cis-driven")

# Parse hier.part results ----
if(length(unique(sapply(deg, function(g) length(g.lm_models[[g]][["hier.part"]])))) > 1){
  print("Some genes do not have trait variance")
  q()
}

hier.part.results <- do.call(rbind.data.frame,
        lapply(deg, function(g) g.lm_models[[g]][["hier.part"]]))
colnames(hier.part.results) <- c(its, "eVariants",  "n.eVariants", "EA","AA")
rownames(hier.part.results) <- deg

if(identical(rownames(hier.part.results), rownames(results.df))){
  results.df <- cbind.data.frame(results.df, hier.part.results)
}

# Save ancestry-DEG classified ----
saveRDS(results.df,
        paste0(outpath,tissue,'.Ancestry_DEG.Classified.rds'))
saveRDS(g.lm_models,
        paste0(outpath,tissue,'.Ancestry_DEG.g.lm_models.rds'))


# Report summary ----
report_df <- do.call(rbind.data.frame,
                     lapply(deg, function(gene)
                       g.lm_models[[gene]][["report_summary"]]))
saveRDS(report_df,
        paste0(outpath,tissue,'.ieQTLs.Report_summary.rds'))

# d ----
saveRDS(d,
        paste0(outpath,tissue,'.Ancestry_DEG.Classification_summary.rds'))

# eVariant not in VCF
saveRDS(eVariants_lost,
        paste0(outpath, tissue, ".eVariants_not_in_vcf.rds"))

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
