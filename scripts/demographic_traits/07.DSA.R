#!/usr/bin/env Rscript

Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()
# ---------------------- #

# Libraries ####
library(caret)
library(parallel)
#print(detectCores())
library(hier.part)
library(sandwich)
library(lmtest)

# Functions ####
source("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")
#source("~/GTEx_v8/Raquel/R_functions/DEA_and_DSA.R_functions.R")

# Command line arguments ####
args <- commandArgs(trailingOnly=TRUE)

# Tissue
tissue <- args[1]
#tissue <- "BrainHypothalamus"

# Paths 
# To save: 
# -- PSI & TPM of alternatively spliced events (ASE) 
# -- hier.part of ASE
# -- PSI residuals of ASE
# -- Differential splicing analysis (DSA) results

outpath <- args[2] #paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/04_DSA/Tissues/",tissue,"/")
outpath <- paste0(outpath, tissue, "/")
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

# Number of CPU (cores) to parallelize mclapply
n_cores <- as.numeric(args[3])

# --------------- #
print("#############################")
print(paste0("Tissue: ", tissue))
print(paste0("Outpath: ", outpath))
print(paste0("n: ", n_cores))
print("#############################")

# Data ####
# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
gene_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#gene_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#length(unique(gene_annotation$ensembl.id)) # 26,196 genes
sex.biased.genes <- gene_annotation[gene_annotation$chr=="chrY" |
                                    gene_annotation$gene.name=="XIST", "ensembl.id"]

# Transcript annotation ----
transcript_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.transcripts.bed", header = F)
#transcript_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "g.biotype", "transcript.name","t.biotype")
transcript_annotation <- transcript_annotation[transcript_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
#length(unique(transcript_annotation$ensembl.id))
#length(unique(transcript_annotation$transcript.id))

# Exon annotation ----
exon_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.exons.bed", header = F)
#exon_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.exons.bed", header = F)
colnames(exon_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "transcript.name","exon.id","exon.numer","g.biotype", "t.biotype")
exon_annotation <- exon_annotation[exon_annotation$ensembl.id %in% gene_annotation$ensembl.id,]
#length(unique(exon_annotation$ensembl.id))

# Event annotation ----
events.info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")
#events.info <- readRDS("~/GTEx_v8/Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")

# Metadata ----
metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
#metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")

# Transcript TPM ----
if(tissue %in% c("BreastMammaryTissue_Female", "BreastMammaryTissue_Male")){
  transcript.tpm <- read.delim(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/SUPPA/TranscriptExpressionFiles/", "BreastMammaryTissue", ".transcript_TPM.txt"))
}else{
  transcript.tpm <- read.delim(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue, ".transcript_TPM.txt"))
}  
#transcript.tpm <- read.delim(paste0("~/GTEx_v8/Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue,".transcript_TPM.txt"))
colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))
# Subset to tissue samples --
transcript.tpm <- transcript.tpm[, metadata$Sample]
identical(colnames(transcript.tpm), metadata$Sample)

# Reading in PSI and TPM values for alternative splicing events annotated in PC and lincRNA genes ----
# SUPPA report values like 1.0000000000000002  & 0.9999999999999998 that when read into R appear as 1 but are not recognized internally as 1
# These types of events would not count as 1 for instance when you try to count the number of samples with PSI == 1
# Included round(psi,2) in previous step to prevent this issue.
psi <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
#psi <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))
#tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))

# Expressed genes ----
exprs.path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/"
#exprs.path <- "~/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/"
exprs_genes <- rownames(readRDS(paste0(exprs.path,tissue,"/",tissue,".voom_limma.covariates_and_traits.results.rds"))[["Age"]])
#exprs_genes <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.Expressed_genes.rds"))


# ----Code ---- #####
is.event.exprs <- function(event.id, tpm.threshold = 0.5){
  if(which(rownames(psi) == event.id) %% 1000 == 0){ system(paste("echo 'Processed: ", which(rownames(psi) == event.id)," out of ", nrow(psi), "events'"))}
  #event.id <- rownames(psi)[1]
  # The two most abundant isforms in numeratos and denominator[!numerator] median TPM >= 1
  # All anotated isoforms
  #transcripts.id <- transcript_annotation[transcript_annotation$ensembl.id == events.info[event.id, "ensembl.id"], "transcript.id"]
  
  # Isoforms that include the event
  isoforms.in <- unlist(strsplit(events.info[event.id, "isoforms.spliced_in"],split = ","))
  # Isoforms that excluded the event
  isoforms.out <- unlist(strsplit(events.info[event.id, "isoforms.spliced_out"],split = ","))
  # Isoforms exprs value for each tissue sample    
  isoforms.tpm <- lapply(metadata$Sample, function(sample)
    sapply(c(isoforms.in, isoforms.out), function(isoform)
      transcript.tpm[isoform,sample]
    ))
  names(isoforms.tpm) <- metadata$Sample
  
  # ---- #
  # Option 1: Most abundant isoform median TPM >= 1 ----
  # Most abundant isoform expression (median per isoform across samples)
  # isoforms.median.tpm <- sapply(c(isoforms.in, isoforms.out), function(isoform)
  #     median(sapply(metadata$Sample, function(sample)
  #         isoforms.tpm[[sample]][isoform]
  #     ))
  # )
  # isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.in]
  # isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.out]
  # 
  # condition <- max(isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.in]) >  tpm.threshold & 
  #     max(isoforms.median.tpm[names(isoforms.median.tpm) %in% isoforms.out]) > tpm.threshold
  # return(condition)
  # ---- #
  
  # ---- #
  # Option 2: At least 20% of the samples express the most abundant isofoorm above 1 TPM ----
  # Most abundant isoform expression (median per isoform across samples)
  iso.in.most_abundant <- names(which.max(sapply(isoforms.in, function(isoform)
    median(sapply(metadata$Sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))
  iso.out.most_abundant <- names(which.max(sapply(isoforms.out, function(isoform)
    median(sapply(metadata$Sample, function(sample)
      isoforms.tpm[[sample]][isoform]
    ))
  )))
  # 20% of the samples express the most abundant above 1 TPM
  condition1 <- sum(sapply(metadata$Sample, function(sample)         
    isoforms.tpm[[sample]][iso.in.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))
  condition2 <- sum(sapply(metadata$Sample, function(sample)         
    isoforms.tpm[[sample]][iso.out.most_abundant] >= tpm.threshold
  )) >= round(0.2*nrow(metadata))
  
  return(condition1 & condition2)
  
}


# 1. Filter AS events ----
print("# ---- Selecting alternatively spliced events ---- #")
print(paste0("Number of ASE in PC and lincRNA genes: ", nrow(psi)))

# * Track number of alternative splicing events ----
no.ase <- vector() 
no.ase <- nrow(psi)

# 1.1 Subset events in PCG and lincRNA expressed in tissue ----
print("# ---- Subset events in PCG and lincRNA expressed in tissue ---- #")

# Keep events in expressed PC and lincRNA --
psi$ensembl_id <- sapply(rownames(psi), function(gene) unlist(strsplit(gene,split = ";"))[[1]])
tpm$ensembl_id <- sapply(rownames(tpm), function(gene) unlist(strsplit(gene,split = ";"))[[1]])
#length(unique(psi$ensembl_id))
psi <- psi[psi$ensembl_id %in% exprs_genes,]
tpm <- tpm[tpm$ensembl_id %in% exprs_genes,]
psi <- psi[!psi$ensembl_id %in% sex.biased.genes,]
tpm <- tpm[!tpm$ensembl_id %in% sex.biased.genes,]
#length(unique(psi$ensembl_id))
psi <- psi[,-ncol(psi)]
tpm <- tpm[,-ncol(tpm)]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.2 Retain events with no samples having an NA ----
print("# ---- Retain events with no samples having an NA ---- #")

# Retain events with no samples having an NA
number_na <- rowSums(is.na(psi),na.rm=F) # vector with number of samples with NA per event
noNA_events.psi <- names(number_na)[number_na==0] # events with 0 NAs
psi <- psi[noNA_events.psi, ]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.3 Calculate event residuals ####
print("# ---- Calculating residuals ---- #")

# 1.3.1 Variables ----
if(! tissue %in% c("Vagina","Uterus","Ovary","Prostate","Testis","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
  individual_traits <- c("Age","Ancestry", "Sex","BMI")
  contrast_names <- c('Age','AncestryAFR','Sex2','BMI')
}else{
  individual_traits <- c("Age","Ancestry", "BMI")
  contrast_names <- c('Age','AncestryAFR','BMI')
}
covariates <- colnames(metadata)[! colnames(metadata) %in% c("Donor", "Sample", individual_traits)]

# 1.3.2 Compute residuals ----
# -------------- #
print(Sys.time())
# -------------- #
# Residuals
fr <- mclapply(rownames(psi), function(event_id) get_residuals(event_id, metadata), mc.cores = n_cores )
names(fr) <- rownames(psi)
# -------------- #
print(Sys.time())
# -------------- #

# 1.3.3 Create dataframe ----
psi_residuals <- do.call(rbind.data.frame,
                         fr)
colnames(psi_residuals) <- colnames(psi)
rownames(psi_residuals) <- rownames(psi)
psi_residuals <- round(psi_residuals, 2)

# 1.4 Exclude events with low complexity ----
print("# ---- Exclude event with low complexity ---- #")
# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi.complexity <- apply(psi, 1, function(x) length(unique(x)))
# exclude events with fewer than max(10, 0.1n) unique values, where n is the sample size
psi_residuals.complexity <- apply(psi_residuals, 1, function(x) length(unique(x)))

psi <- psi[intersect(names(psi.complexity[psi.complexity >= 10]),
                              names(psi_residuals.complexity[psi_residuals.complexity >= 10])
),]
psi_residuals <- psi_residuals[intersect(names(psi.complexity[psi.complexity >= 10]),
                                         names(psi_residuals.complexity[psi_residuals.complexity >= 10])
                                         ),]


# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.5 Exclude events with with insufficient variability ----
print("# ---- Exclude events with with insufficient variability ---- #")
event_freq <- apply(psi, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
event_freq_residuals <- apply(psi_residuals, 1, function(x) sort(table(x),decreasing = T)[1]/sort(table(x),decreasing = T)[2] < 80/20)
#excluded.psi <- psi[which(event_freq==F),]
psi <- psi[event_freq & event_freq_residuals,]
psi_residuals <- psi_residuals[event_freq & event_freq_residuals,]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 1.6 Exclude events not sufficiently expressed ----
print("# ---- Exclude events not sufficiently expressed ---- #")
# -------------- #
print(Sys.time())
# -------------- #
events_exprs <- unlist(mclapply(rownames(psi), function(i) is.event.exprs(i, 0.5),  mc.cores = n_cores  ))
# -------------- #
print(Sys.time())
# -------------- #

# Subset ASE sufficiently exprs ----
psi <- psi[events_exprs,]
psi_residuals <- psi_residuals[rownames(psi),]
tpm <- tpm[rownames(psi),]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 2. hier.part ----
print("# ---- Running hier.part ---- #")

# 2.1 Calculate explained variance ----
# -------------- #
print(Sys.time())
# -------------- #
hier.part.results <- mclapply(rownames(psi_residuals), function(event)
  hier.part.mod(y=as.numeric(psi_residuals[event,]), x=metadata[,individual_traits], 
                fam = "quasibinomial", link = "logit", gof = "Rsqu",control = list(maxit = 100)), mc.cores = n_cores)
names(hier.part.results) <- rownames(psi_residuals)
# -------------- #
print(Sys.time())
# -------------- #

# 2.2 Exclude events with negative estimates ----
print(paste0("ASE events with unestimable contributions: ", sum(is.na(hier.part.results))))
psi <- psi[-which(is.na(hier.part.results)),]
psi_residuals <- psi_residuals[-which(is.na(hier.part.results)),]
hier.part.results <- hier.part.results[-which(is.na(hier.part.results))]

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi))

# 2.3 Parse results ----
rsq <- sapply(rownames(psi_residuals), function(event) 
  sum(hier.part.results[[event]]$IJ[,1])) 
names(rsq) <- rownames(psi_residuals)
rel_perc <- do.call(rbind.data.frame,
                    lapply(names(hier.part.results), function(event) 
                      as.numeric(unlist(hier.part.results[[event]]$I.perc))))
rownames(rel_perc) <- rownames(psi_residuals)
colnames(rel_perc) <- individual_traits
abs_perc <-  do.call(rbind.data.frame,
                     lapply(names(hier.part.results), function(event) 
                       hier.part.results[[event]]$IJ[,1])
)
rownames(abs_perc) <- rownames(psi_residuals)
colnames(abs_perc) <- individual_traits
hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))


# 3. Differential splicing analysis ----
print("# ---- Running differential splicing analysis ---- #")

# 3.1 Run PSI models ----
# -------------- #
print(Sys.time())
# -------------- #
# One model per event
fr <- mclapply(rownames(psi), function(event_id) model_psi(event_id, metadata), mc.cores = n_cores )
names(fr) <- rownames(psi)
# -------------- #
print(Sys.time())
# -------------- #

# 3.2 Parse results table ----
results <- do.call(rbind.data.frame,
                   lapply(rownames(psi), function(event_id) fr[[event_id]][['res']]))
dsa_res <- lapply(individual_traits, function(trait)
  results[results$Trait==trait,])
names(dsa_res) <- individual_traits

for(trait in individual_traits){
  # Set event ID as as rownames
  rownames(dsa_res[[trait]]) <- dsa_res[[trait]]$Event_id
  dsa_res[[trait]] <- dsa_res[[trait]][,-1]
}

# 3.3 Exclude events with warnings ----
# -- Beware of warnings -- #
# The coeftest function raises warning at particular events when the variance-covariance matrix cannot be computed
# In extreme instances it returns NA, for example, instances of events completely stratified that we excluded before-hand
# Yet, some events might remain that cannot be modelled
# If coef.test raises a warning we consider those events are not properly modelled and assign a P.Value of NA
# Multiple testing pvalue correction ->
# if in the p-values vector (not corrected) there is an NA, the p.adjust() does not consider it in the n (number of observations)
# Remove events with warnings in coeftest function
# These warnings appear in tissues with low sample size for events with low variance
# We consider these events cannot be modelled
print(paste0("ASE that raised glm warning: ", sum(dsa_res$Age$glm_Warning == "1") ))
print(paste0("ASE that raised coef.test warning: ", sum(dsa_res$Age$coeftest_Warning == "1") ))

# * Track number of alternative splicing events ----
no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1"))
no.ase <- c(no.ase, nrow(psi) - sum(dsa_res$Age$glm_Warning == "1") - sum(dsa_res$Age$coeftest_Warning == "1"))
no.ase <- c(no.ase, nrow(dsa_res[["Age"]]))

if(sum(dsa_res$Age$glm_Warning == "1") > 0 ){
  for(trait in individual_traits){
    dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$glm_Warning=="1")] <- NA
  }
}
if(sum(dsa_res$Age$coeftest_Warning == "1") > 0 ){
  for(trait in individual_traits){
    dsa_res[[trait]]$P.Value[which(dsa_res[[trait]]$coeftest_Warning=="1")] <- NA
  }
}
for(trait in individual_traits){
  dsa_res[[trait]] <- dsa_res[[trait]][!is.na(dsa_res[[trait]]$P.Value),]
}

# 3.4 FDR correction ---
for(trait in individual_traits){
  dsa_res[[trait]]$adj.P.Val <- p.adjust(dsa_res[[trait]]$P.Value, method = "BH")
}
print(paste0("ASE modelled: ", nrow(dsa_res$Age)))

# 3.5 Save results ----
trait <- "Age"
psi <- psi[rownames(dsa_res[[trait]]),]
psi_residuals <- psi_residuals[rownames(dsa_res[[trait]]),]
tpm <- tpm[rownames(dsa_res[[trait]]),]
hier_data <- hier_data[rownames(dsa_res[[trait]]),]

# PSI values
saveRDS(psi, 
        paste0(outpath,tissue,".Alternatively_spliced_events.PSI_values.rds"))
# TPM values
saveRDS(tpm, 
        paste0(outpath,tissue,".Alternatively_spliced_events.TPM.rds"))
# Save residuals
saveRDS(psi_residuals,paste0(outpath, tissue,".Alternatively_spliced_events.psi_residuals.rds"))
# Save hier.part
saveRDS(hier_data,paste0(outpath,tissue, ".Alternatively_spliced_events.hier_part.rds"))
# Save filtering
saveRDS(no.ase,
        paste0(outpath,tissue,".Alternatively_spliced_events.Filtering.rds"))
saveRDS(dsa_res,
        paste0(outpath,
               tissue,".fractional_regression.covariates_and_traits.results.tmp.rds")	
)

# 4. DeltaPSI and average TPM ----
print("# ---- Calculating deltaPSI ---- #")

# 5.2 Compute average TPM expression and variance for each event
avrg_TPM <- apply(tpm, 1, function(x) mean(log2(x+1)))
var_TPM <- apply(tpm, 1, function(x) var(x))

AFR_samples <- metadata[metadata$Ancestry=="AFR", "Sample"]
EUR_samples <- metadata[metadata$Ancestry=="EUR", "Sample"]
Female_samples <- metadata[metadata$Sex=="2", "Sample"]
Male_samples <- metadata[metadata$Sex=="1", "Sample"]
obese_samples <- metadata[metadata$BMI>=30, "Sample"]
normal_samples <- metadata[metadata$BMI<25, "Sample"]

# 5.2 Compute deltaPSI ----
get.deltaPSI <- function(trait){
  print(trait)
  # -------------- #
  print(Sys.time())
  # -------------- #)
  
  if(trait == "Ancestry"){
    # delta.psi <- unlist(mclapply(rownames(dsa_res[[trait]]), function(event)
    #     mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T),
    #     mc.cores = n_cores
    # ))
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,EUR_samples]),na.rm=T) - mean(as.numeric(psi[event,AFR_samples]),na.rm=T))
  }else if(trait == "Sex"){
    if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
      delta.psi <- NA
    }else{
      # delta.psi <-  unlist(mclapply(rownames(dsa_res[[trait]]), function(event)
      #     mean(as.numeric(psi[event,Female_samples]),na.rm=T) - mean(as.numeric(psi[event,Male_samples]),na.rm=T),
      #     mc.cores = n_cores
      # ))
      delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
        mean(as.numeric(psi[event,Female_samples]),na.rm=T) - mean(as.numeric(psi[event,Male_samples]),na.rm=T) )
    }
  }else if(trait == "Age"){
    younger_samples <- metadata[metadata$Age<40,"Sample"] # [20-39]
    older_samples <- metadata[metadata$Age>=50,"Sample"] # [50-70]
    # delta.psi <- unlist(mclapply(rownames(dsa_res[[trait]]), function(event)
    #     mean(as.numeric(psi[event,older_samples]),na.rm=T) - mean(as.numeric(psi[event,younger_samples]),na.rm=T),
    #     mc.cores = n_cores
    # ))
    delta.psi <- sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,older_samples]),na.rm=T) - mean(as.numeric(psi[event,younger_samples]),na.rm=T))
  }else if(trait == "BMI"){
    # delta.psi <-  unlist(mclapply(rownames(dsa_res[[trait]]), function(event)
    #     mean(as.numeric(psi[event,obese_samples]),na.rm=T) - mean(as.numeric(psi[event,normal_samples]),na.rm=T),
    #     mc.cores = n_cores
    # ))
    delta.psi <-  sapply(rownames(dsa_res[[trait]]), function(event)
      mean(as.numeric(psi[event,obese_samples]),na.rm=T) - mean(as.numeric(psi[event,normal_samples]),na.rm=T))
  }else{
    delta.psi <- NA
  }
  # -------------- #
  print(Sys.time())
  # -------------- #
  return(delta.psi)
}

# Add deltaPSI, AvgExprs, ExprsVar and order data.frame
for(trait in individual_traits){
  if(trait == "Sex"){
    if(tissue %in% c("Ovary","Uterus","Vagina","Testis","Prostate","BreastMammaryTissue_Female","BreastMammaryTissue_Male")){
      dsa_res[[trait]] <- NA
    }else{
      # Add deltaPSI
      dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
      # Add average expression TPM
      dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
      dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
      # Order by FDR
      dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
    }
  }else{
    # Add deltaPSI
    dsa_res[[trait]]$deltaPSI <- get.deltaPSI(trait)
    # Add average expression TPM and variance
    dsa_res[[trait]]$AvgExprs <- sapply(rownames(dsa_res[[trait]]), function(event) avrg_TPM[event])
    dsa_res[[trait]]$ExprsVar <- sapply(rownames(dsa_res[[trait]]), function(event) var_TPM[event])
    # Order by FDR
    dsa_res[[trait]] <- dsa_res[[trait]][order(dsa_res[[trait]]$adj.P.Val),]
  }
}

# 5.3 Save data ----
saveRDS(dsa_res,
        paste0(outpath,
               tissue,".fractional_regression.covariates_and_traits.results.rds")	
)

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

