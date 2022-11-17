#!/usr/bin/env Rscript

# Load libraries ####
suppressMessages(library(edgeR)) #Already includes limma
library(optparse)
options(warn=-1)

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-t", "--tissue"), type="character",
                     dest="tissue",
                     help="Tissue")
parser <- add_option(parser, opt_str=c("-d", "--diseases"), type="character",
                     dest="diseases",
                     help="Acronyms separated by _")
parser <- add_option(parser, opt_str=c("-o", "--outdir"), type="character", default = "~/Documents/mn4/Jose/03_Models/Tissues/",
                     dest="out", help="Output directory where all the diseases folders are")

options=parse_args(parser)
tissue=options$tissue
tissue <- gsub("\302","",tissue)
diseases=options$diseases
out=options$out

# tissue <- "NerveTibial"
# diseases <- "MHT1D_MHT2D"
# out <- "~/Documents/mn4/Jose/03_Models/Tissues/"

# tissue <- "MuscleSkeletal"
# diseases <- "Atrophy_MHT1D_MHT2D"
# tissue <- "AdiposeSubcutaneous"
# diseases <- "MHT1D_MHT2D"
# tissue <- "Lung"
# diseases <- "Atelectasis_Emphysema_Fibrosis_Pneumonia_Hist_MHT1D_MHT2D"
# out <- "~/Documents/mn4/Jose/03_Models/Testing_Diabetes/final/"

name <- diseases
diseases <- strsplit(diseases, "_")[[1]]
if("Atherosclerotic" %in% diseases){
  where <- grep("Atherosclerotic", diseases)
  diseases[where] <- "Atherosclerotic_Arteries" 
  diseases <- diseases[!diseases=="Arteries"]
}
if("Ischemic" %in% diseases){
  where <- grep("Ischemic", diseases)
  diseases[where] <- "Ischemic_changes" 
  diseases <- diseases[!diseases=="changes"]
}
if("Pneumonia" %in% diseases){
  where <- grep("Pneumonia", diseases)
  diseases[where] <- "Pneumonia_Hist" 
  diseases <- diseases[!diseases=="Hist"]
}
if("Post" %in% diseases){
  where <- grep("Post", diseases)
  diseases[where] <- "Post_menopausal" 
  diseases <- diseases[!diseases=="menopausal"]
}
folder_name <- tissue

print(name)
print(tissue)

# Reading metadata
#dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"

# dir <- sub("/[^/]+$", "", out) #Not working anymore
dir <- strsplit(out, "/")[[1]]
dir <- dir[1:(length(dir)-2)]
dir <- paste0(dir, collapse="/")


###### CHANGED
outpath <- paste0(out, folder_name,"/") #Only moment where variable out is used, 

##Now, I change folder_name for diabetes testing
# folder_name <- "balanced"

metadata <- readRDS(paste0(outpath,tissue,".SampleMetadata.", name,".rds")) 
#Moment to exclude variables if we want to check one diseases with the n samples as if many diseases
#metadata <- metadata[,-23]
#diseases <- diseases[-5]

if(length(unique(metadata$HardyScale))<5){metadata$HardyScale <- droplevels(metadata$HardyScale)}
if(length(unique(metadata$Ancestry))<2){metadata$Ancestry <- droplevels(metadata$Ancestry)}
if("Sex" %in% colnames(metadata)){
  if(length(unique(metadata$Sex))<2){metadata$Sex <- droplevels(metadata$Sex)}
}

main_data <- read.csv(paste0(dir, "/00_Data/main_final.csv"))
tissues <- main_data$Tissue

# Functions ####

limma_lm2 <- function(fit, covariate, covariate_data){
#  print(covariate)
  if(covariate %in% nonEstimable(fit$design)){
    print("Covariate non estimable")
    print(covariate)
    return(0)
  }
  v.contrast <- rep(0,ncol(fit$design))
  v.contrast[ which( colnames(fit$design) == covariate) ] <- 1
  contrast.matrix <- cbind( "C1" = v.contrast)
  fitConstrasts <- contrasts.fit(fit,contrast.matrix)
  eb = eBayes(fitConstrasts)
  tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
  return(tt.smart.sv)
}


# Counts

# dir <- sub("/[^/]+$", "", dir)

dir <- strsplit(dir, "/")[[1]]
dir <- dir[-(length(dir))]
dir <- paste0(dir, collapse="/")


counts <- readRDS(paste0(dir, "/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
tpm <- readRDS(paste0(dir, "/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.PC_lincRNA.rds"))


colnames(counts) <- sapply(colnames(counts), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
colnames(tpm) <- sapply(colnames(tpm), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )

if(nchar(metadata$Sample[1])>18){ #The donor disease have the whole sample name. 
  metadata$Sample <- sapply(metadata$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
}

# Subset samples with disease annotation data
counts <- counts[,colnames(counts) %in% metadata$Sample]
tpm <- tpm[,colnames(tpm) %in% metadata$Sample]


# 2. Remove genes that are lowly expressed ###
# 1. TPM>=0.1 in at least 20% of the tissue samples. Changed to >=1
exprs_genes.tpm <- rownames(tpm)[apply(tpm, 1, function(x) sum(x>=1) ) >= 0.2*ncol(tpm)  ]

# 2. Count >=6 in at least 20% of the tissue samples. Changed to 10
exprs_genes.counts <- rownames(counts)[ apply(counts, 1, function(x) sum(x>=10) ) >= 0.2*ncol(counts)  ]

# 3. Intersect gene lists
exprs_genes <- intersect(exprs_genes.tpm,
                         exprs_genes.counts)

# Exclude chrY genes in female-only tissues
# Gene Annotation
gene_annotation  <- read.delim(paste0(dir, "/Jose/00_Data/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")) #It only contains protein_coding and lincRNA, and has an extra column with ensembl id without ".". No PAR genes


Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id
#The only disease with tissue of origin in a sex tissue is spermatogenesis and it is only in testis, so I don't need to exclude Y_genes from any disease now

# Save gene list
saveRDS(exprs_genes, paste0(outpath,tissue,".SelectedSamples.Expressed_genes.rds"))


#  Normalising gene expression distributions
# library(edgeR)

#We need the counts to have the same sample order as metadata!
#change order of counts samples
# numbers <- match(metadata$Sample, colnames(counts))
# counts_new <- counts[,numbers]
# counts <- counts_new
# Create DGEList object
dge <- DGEList(counts[exprs_genes,])

# Calculate normalization factors (does not do the normalization yet, only computes the factors)
dge <- calcNormFactors(dge)

# 4. voom + limma. 

# Voom
v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates

#New model

sex_tissues <- c("Vagina","Uterus","Ovary","Prostate","Testis", "BreastMammaryTissue_Female", "BreastMammaryTissue_Male")

if(tissue %in% sex_tissues){
  individual_traits <- c("Age","Ancestry", "BMI", diseases)
} else {
  individual_traits <- c("Age","Ancestry", "Sex","BMI", diseases)
}

covariates <- c(colnames(metadata)[!colnames(metadata) %in% c("Donor", "Sample", individual_traits)])

if("Atherosclerotic_Arteries" %in% diseases | "Atherosclerotic_Arteries" %in% covariates){
  covariates <- covariates[!covariates %in% c("Atherosclerosis", "Atherosis", "Calcification", "Sclerotic", "Steatosis")]
}

fml_args_mod <- paste(c(covariates, individual_traits), collapse = " + ")
print(fml_args_mod)

mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)

# Limma fit
fit <- lmFit(v, mod)

# Limma test
individual_traits <- sapply(individual_traits, function(trait) 
  if(is.factor(metadata[[trait]])){
    return(paste0(trait, "1"))
  } else{return(trait)}
, USE.NAMES = F)

print(individual_traits)

dea_res <- lapply(individual_traits, function(phenotype) limma_lm2(fit,phenotype,metadata ) )
names(dea_res) <- individual_traits

# Save table with results

print(paste0(outpath,tissue,".", name,".voom_limma.results_PEER.rds"))


saveRDS(dea_res, paste0(outpath,tissue,".", name,".voom_limma.results_PEER.rds"))
