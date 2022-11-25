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
library(gtools)
library(hier.part)

# Functions ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/Splicing_R_functions.R")
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/Splicing_R_functions.R")
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
#inpath <- "~/GTEx_v8/Raquel/Draft/02.DiffSplic/01.DSA/Tissues/"
inpath <- paste0(inpath, tissue, "/")
outpath <- args[3] 
outpath <- paste0(outpath, tissue, "/")
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

# Number of CPU (cores) to parallelize mclapply
#n_cores <- args[3]

print("#### ----- ####")
print(paste0("Tissue: ", tissue))
print(paste0("Inpath: ", inpath))
print(paste0("Outpath: ", outpath))
print("#### ----- ####")

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

# ---- Code ----
# Model each event as:
# psi_residual ~ sVariants_in_sGene + trait
# The genotypes of the sVariants are modelled as a factor, where: 0/0 -> 0, 0/1 -> 1; 1/1 ->2

### Prepare data ####

# Tissue info ----
tissue_info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissue_info.rds")
#tissue_info <- readRDS("~/GTEx_v8/Raquel/Draft/Data/Tissue_info.rds")
tissue_id <- tissue_info[tissue_info$tissue_ID==tissue,]$tissue_id

# sGenes ----
inpath_sqtls <- "/gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/"
#inpath_sqtls <- "~/GTEx_v8_data/cisQTLs/"
sgene_data <- as.data.frame(data.table::fread(paste0(inpath_sqtls,"sGenes/GTEx_Analysis_v8_sQTL/",tissue_id,".v8.sgenes.txt.gz")))

# Select sGenes (q.val < 0.05) --
sGenes <- sgene_data[sgene_data$qval <= 0.05,"gene_id"]

# sGenes:isQTL ----
# independent sQTLs -> sVariants (associated which sGenes, qval<=0.05)) 
sqtls_dir <- "/gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/"
#sqtls_dir <- "~/GTEx_v8_data/cisQTLs/GTEx_Analysis_v8_sQTL_independent/"
if(tissue %in% c("BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    sqtls.tissue <-  as.data.frame(data.table::fread(paste0(sqtls_dir, "Breast_Mammary_Tissue", ".v8.independent_sqtls.txt.gz")))
}else{
    sqtls.tissue <-  as.data.frame(data.table::fread(paste0(sqtls_dir, tissue_id,".v8.independent_sqtls.txt.gz")))    
}

# sGenes with at leats 1 independent sQTL ----
sGenes.isqtl <- unique(sqtls.tissue$group_id)

# Differential splicing analysis: dataframe with results ----
dsa_res <- readRDS(paste0(inpath, tissue,".fractional_regression.covariates_and_traits.results.rds"))[["Ancestry"]]

# Save summary data ----
d <- list()
d[["sGenes"]] <- length(sGenes)
d[["sGenes:isQTL"]] <- length(sGenes.isqtl)
d[["Ancestry:DSE"]] <- sum(dsa_res$adj.P.Val < 0.05)
d[["Ancestry:DSG"]] <- length(unique(dsa_res[dsa_res$adj.P.Val < 0.05,'Ensembl_id']))
d[["Ancestry:DSE:sGenes"]] <- sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"]  %in% sGenes)
d[["Ancestry:DSE:sGenes-isQTL"]] <- sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"]  %in% sGenes.isqtl)

# Print to std.out summary data ----
print(paste0("No. of sGenes: ", length(sGenes)))
print(paste0("No. of sGenes-isQTLs:  ", length(sGenes.isqtl)))
print(paste0("No. of ancestry DSE:  ", sum(dsa_res$adj.P.Val < 0.05)))
print(paste0("No. of ancestry DSG:  ", length(unique(dsa_res[dsa_res$adj.P.Val < 0.05,'Ensembl_id']))))
print(paste0("No. of ancestry DSE (SE, MX, A5, A3, AF, AL):  ", nrow(dsa_res[dsa_res$adj.P.Val < 0.05 & 
                                                                                 dsa_res$Type %in% c("SE","MX","A5","A3","AF","AL"),]))) 
print(paste0("No. of ancestry DSE (RI):  ", nrow(dsa_res[dsa_res$adj.P.Val < 0.05 & 
                                                             dsa_res$Type %in% c("RI"),]))) 
print(paste0("No. of ancestry DSE in sGenes: ", sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"]
                                                    %in% 
                                                        sGenes    
),
" (", round(  sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"] %in% sGenes)/sum(dsa_res$adj.P.Val < 0.05),2) ,")"
))
print(paste0("No. of ancestry DSE in sGenes-isQTL: ", sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"]
                                                          %in% 
                                                              sGenes    
),
" (", round(  sum(dsa_res[dsa_res$adj.P.Val < 0.05,"Ensembl_id"] %in% sGenes)/sum(dsa_res$adj.P.Val < 0.05),2) ,")"
))

# Goal: select differentially spliced events (SE|MX|A5|A3|AF|AL) that affect sGenes with at least 1 genotyped independent sQTL (sVariant) with MAF > 0.001 ----
dsa_res <- dsa_res[dsa_res$adj.P.Val < 0.05 &
                       dsa_res$Ensembl_id %in% sGenes.isqtl,]
print(paste0( nrow(dsa_res), 
              " (", round(nrow(dsa_res)/d$`Ancestry:DSE`, 2)*100, "%) of Ancestry DSE are Events in sGenes-isQTL"))

# Subset sQTLs in sGenes differentially spliced with at least 1 independent sQTL ----
sqtls.tissue <- sqtls.tissue[sqtls.tissue$group_id %in% unique(dsa_res$Ensembl_id),]

# sVariants in sGenes differentially spliced with at leats 1 independent sQTL ----
sVariants <- unique(sqtls.tissue$variant_id) # isQTLs
print(paste0("No. of sVariants (i-sQTL) in ancestry DSG sGenes: ", length(sVariants)))

# Genotyped sVariants ----
# Read the genotypes for the sVariants which are also independent cis-sQtls 
# Job run in cluster
#vcftools --gzvcf ${vcf_file} --snps <(zcat /gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/${tissue_id}.v8.independent_sqtls.txt.gz  | cut -f7 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.txt.gz

# Read in table data --
gt_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Fst/isQTL_GT/"
#gt_path <- "~/GTEx_v8/Raquel/Draft/Data/Fst/isQTL_GT/"
if(!tissue %in% c("BreastMammaryTissue_Female", "BreastMammaryTissue_Male")){
    variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, tissue, ".isQTLs.Donor_genotypes.txt.gz")))  
}else{
    variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, "BreastMammaryTissue", ".isQTLs.Donor_genotypes.txt.gz"))) 
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

if(length(sVariants) > 1){
    # Subet genotype file and only keep sVariants in eGenes differentially expressed --
    variants_genotype <- variants_genotype[variants_genotype$variant_id %in% sVariants,] 
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
if(length(sVariants) == 1){
    # Subet genotype file and only keep sVariants in eGenes differentially expressed --
    variants_genotype <- variants_genotype[variants_genotype$variant_id %in% sVariants,]
    vg.tissue_mod <- vg.tissue_mod[vg.tissue_mod$variant_id %in% sVariants,]
    vg.tissue_wide <- vg.tissue_wide[vg.tissue_wide$variant_id %in% sVariants,]
    vg.tissue_long <- vg.tissue_long[vg.tissue_long$variant_id %in% sVariants,]
}

# Subset tissue samples ----
vg.tissue_long <- vg.tissue_long[vg.tissue_long$Individual_ID %in% metadata$Donor,]

# sVariants in VCF -> sVariants in eGenes with at leats 1 independent eQTL with MAF > 0.001 ----
sVariants_lost <- sVariants[!sVariants %in% unique(vg.tissue_long$variant_id)] # MAF < 0.001?
sVariants <- sVariants[sVariants %in% unique(vg.tissue_long$variant_id)]
print(paste0("No. of sVariants (i-eQTL) with MAF > 0.01 in ancestry DSE sGenes: ", length(sVariants)))
print(paste0("Minimum MAF value of genotyped sVariants: ",  min(sqtls.tissue[sqtls.tissue$variant_id %in% sVariants, "maf"])))

# Add info to summary data ----
d[["isQTL:MAF_b_01"]] <- length(sVariants_lost)
d[["isQTL:MAF_01"]] <- length(sVariants)
d[["minMAF"]] <-  min(sqtls.tissue[sqtls.tissue$variant_id %in% sVariants, "maf"])

# Subset sQTLs in sGenes differentially spliced with at least 1 genotyped independent sQTL with MAF > 0.001 ----
sqtls.tissue <- sqtls.tissue[sqtls.tissue$variant_id %in% sVariants,]
isGenes.ds.maf01 <- unique(sqtls.tissue$group_id)

# Add info to summary data ----
d[["isQTL:MAF_b_001"]] <- length(sVariants_lost)
d[["isQTL:MAF_001"]] <- length(sVariants)
d[["minMAF"]] <-  min(sqtls.tissue[sqtls.tissue$variant_id %in% sVariants, "maf"])

# Subset sQTLs in sGenes differentially spliced with at least 1 genotyped independent sQTL with MAF > 0.001 ----
sqtls.tissue <- sqtls.tissue[sqtls.tissue$variant_id %in% sVariants,]
isGenes.ds.maf01 <- unique(sqtls.tissue$group_id)

# Subset event with isQTL with MAF 001 (in vcf file) ----
dsa_res <- dsa_res[dsa_res$Ensembl_id %in% isGenes.ds.maf01,]
d[["Ancestry:DSE:sGenes-isQTL:MAF001"]] <- nrow(dsa_res)
d[["Ancestry:DSG:sGenes-isQTL:MAF001"]] <- length(isGenes.ds.maf01)

# Create 'dictionary' (named list):  ----
gene_variants.list <- sapply(isGenes.ds.maf01, function(x){
    variant <- unique(sqtls.tissue[sqtls.tissue$group_id==x,]$variant_id)
    return(variant)
}, simplify=F)

############### Generalize linear model for each event #########################

## To build the model:
# 1. For each event, recover the gene it affects  ----
# 2. For each gene, recover its sVariants (isQTLS with MAF 001) ----
# 3. Compare 2 models:
# residuals ~ Age + Sex + BMI + cis-effects (isQTL(s))
# residuals ~ Age + Sex + BMI + cis-effects (isQTL(s)) + Ancestry 
# Is Ancestry adding someting that the sQTLs did not capture, a not cis-driven effect?

# Events ----
events <- rownames(dsa_res)

# Events' gene -----
events.genes <- lapply(rownames(dsa_res), function(event) dsa_res[event,"Ensembl_id"])
names(events.genes) <- rownames(dsa_res)

# Modify to match  --
colnames(metadata)[1] <- "Individual_ID"
colnames(metadata)[2] <- "Sample_ID"
metadata <- metadata[,c("Individual_ID", "Sample_ID",  individual_traits)]
it <- "Ancestry"
if(!tissue %in% c("Uterus","Ovary","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
    its <- c("Age","Sex","BMI","Ancestry")    
}else{
    its <- c("Age","BMI","Ancestry")
}

# Recover the PSI residuals associated to those events ----
psi_residuals <- readRDS(paste0(inpath, tissue,".Alternatively_spliced_events.psi_residuals.rds"))

# Susbet residuals of events to be modelled ----
psi_residuals <- psi_residuals[rownames(dsa_res),]    

if(!identical(metadata$Sample_ID, colnames(psi_residuals))){
    print("The samples in the metadata and the residual data do not coincide")
    q()
}
if(!identical(metadata$Individual_ID, as.character(unique(vg.tissue_long$Individual_ID)))){
    print("The donors in the metadata and the genotyped data do not coincide")
    q()
}

# Is the event cis-driven or not cis-driven ----
print('*****************************************************')
print(paste0(tissue, ' has a total of ', length(events), ' DS events with cis-independent sQTLs with MAF > 0.001'))
print('Applying glm() function to each of those events')
print('*****************************************************')
cat('\n')
cat('\n')

# glm per event ----
e.glm_models <- sapply(events, function(e) glm.cis_models(e), simplify = F)
names(e.glm_models) <-  events
# events that  cannot modelled ----
d[["Ancestry:DSE:NotModelled"]] <- sum(is.na(e.glm_models))
# events modelled ----
e.glm_models <- e.glm_models[!is.na(e.glm_models)] # SNP with no variance or no SNP not correlated
d[["Ancestry:DSE:Modelled"]] <- length(e.glm_models)

# Compare modelA vs modelB ----
dse <- names(e.glm_models)
# Anova
anova.P.Value <- sapply(dse, function(e) anova(e.glm_models[[e]][["glm.modelA"]],e.glm_models[[e]][["glm.modelB"]], test="F")[2,6])
# Multiple testing correction
anova.adj.P.Val <- p.adjust(anova.P.Value, method = "BH")

# Parse results data ----
results.df <- cbind.data.frame(dse, anova.adj.P.Val)
results.df$Class <- ifelse(results.df$anova.adj.P.Val < 0.05, "Not_cis-driven","Cis-driven")
results.df$ensembl.id <- sapply(rownames(results.df), function(i) 
    unlist(strsplit(i, split = ";"))[[1]]
)
results.df$Type <- sapply(rownames(results.df), function(i) 
    unlist(strsplit(unlist(strsplit(i, split = ";"))[[2]], split = ":"))[[1]]
) 

# Parse hier.part results ----
if(length(unique(sapply(dse, function(e) length(e.glm_models[[e]][["hier.part"]])))) > 1){
    print("Some genes do not have trait variance")
    q()
}

hier.part.results <- do.call(rbind.data.frame,
                             lapply(dse, function(e) e.glm_models[[e]][["hier.part"]]))
colnames(hier.part.results) <- c(its, "sVariants",  "n.sVariants", "EA","AA")
rownames(hier.part.results) <- dse

if(identical(rownames(hier.part.results), rownames(results.df))){
    results.df <- cbind.data.frame(results.df, hier.part.results)
}


# library(vcd)
# mosaic(ct, 
#        split_vertical = TRUE,
#        main = tissue)
# vcd::assoc(table(results.df$Type, results.df$Class),
#            shade=T, main="Types of DSE in ribosomal proteins")

# Save ancestry-DSE classified ----
saveRDS(results.df,
        paste0(outpath,tissue,'.Ancestry_DSE.Classified.rds'))
# saveRDS(g.df,
#         paste0(outpath,tissue,'.Ancestry_DSG.Classified.rds'))
saveRDS(e.glm_models,
        paste0(outpath,tissue,'.Ancestry_DSE.e.glm_models.rds'))

# Report summary ----
report_df <- do.call(rbind.data.frame,
                     lapply(dse, function(event)
                         e.glm_models[[event]][["report_summary"]]))
saveRDS(report_df,
        paste0(outpath,tissue,'.isQTLs.report_summary.rds'))

# d ----
saveRDS(d,
        paste0(outpath,tissue,'.Ancestry_DSE.Classification_summary.rds'))
#saveRDS(gene.counts,
#        paste0(outpath,tissue,'.Ancestry_DSG.Classification_summary.rds'))

# eVariant not in VCF
saveRDS(sVariants_lost,
        paste0(outpath, tissue, ".sVariants_not_in_vcf.rds"))


# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #
