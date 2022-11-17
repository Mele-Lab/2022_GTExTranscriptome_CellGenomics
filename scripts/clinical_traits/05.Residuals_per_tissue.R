#!/usr/bin/env Rscript
rm(list=ls())

# Load libraries ####
suppressMessages(library(edgeR))
library(hier.part)
library(ggplot2)
library(ggpubr)


first_dir <- "~/Documents/mn4/"
# first_dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"
dir <- paste0(first_dir, "Jose/03_Models/Tissues/")

tissues <- list.dirs(dir, full.names = F)[-1]

#Functions
current.model.mod <- function (y, current.comb, xcan, SStot=0,family = c("gaussian","quasibinomial"), 
                               link = c("logit"), gof = c("Rsqu","RMSPE"), ...)  {
  comb.data <- data.frame(xcan[, current.comb])
  colnames(comb.data) <- colnames(xcan)[current.comb]
  data <- data.frame(y, comb.data)
  depv <- names(data)[1]
  n.comb <- dim(comb.data)[2]
  xs <- vector("character", n.comb)
  for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
                                           "+", sep = "")
  xs[n.comb] <- names(comb.data)[n.comb]
  xss <- paste(xs, collapse = " ", sep = "")
  formu <- stats::formula(paste(depv, "~", xss, sep = ""))
  if (gof == "RMSPE") gf <- sqrt(sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu") {
    if (family == "quasibinomial") 
      gf <- (SStot-sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))/SStot
    if (family == "gaussian") 
      gf <- summary(stats::lm(formu, data = data))$r.squared
  }
  gf
}
all.regs.mod <- function (y, xcan, family = c("gaussian", "quasibinomial"), link = c("logit"), gof = c("Rsqu","RMSPE"),...) { 
  if (sum(is.na(xcan)) > 0) {
    missing <- is.na(apply(xcan, 1, FUN = sum))
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
            call. = FALSE)
  }
  if (sum(is.na(y)) > 0) {
    missing <- is.na(y)
    xcan <- xcan[!missing, ]
    y <- y[!missing]
    warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
            call. = FALSE)
  }
  pcan <- dim(xcan)[2]
  n <- (2^pcan) - 1
  combs <- combos1(pcan)$ragged
  SStot <- sum((y-mean(y))^2)
  
  if (gof == "RMSPE")  gfs <- sqrt(sum((stats::glm(y ~ 1, family = family,...)$fitted.values - y)^2))
  if (gof == "Rsqu")   gfs <- 0
  
  for (i in 1:n) {
    if (i%%500 == 0) 
      cat(i, "regressions calculated:", n - i, "to go...\n")
    current.comb <- as.vector(combs[i, ][combs[i, ] > 0]) 
    combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse = "")
    if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
    if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
    gfs <- c(gfs, new.line)
  }
  gfs
}
hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2] #5
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
}




#Do a for loop 

for(i in 1:length(unique(tissues))){
  tissue <- unique(tissues)[i]
  folder <- tissue
  print(paste0("Computing hier part for ", folder))
  
  outdir <- paste0(dir, folder)
  files <- list.files(outdir)
  metadata_files <- files[grep("SampleMetadata", files)]
  voom_files <- files[grep("voom", files)]
  for(file in metadata_files){ #In the previous version we had more than one metadata per tissue

    metadata <- readRDS(paste0(outdir, "/", file))
    name <- strsplit(file, ".", fixed = T)[[1]][3] #These next lines are just to retrieve expressed genes
    dea_file <- voom_files[grep(name, voom_files)]
    dea_file <- dea_file[grep("interactions", dea_file, invert=T)]
    dea_res <- readRDS(paste0(outdir, "/", dea_file))
    exprs_genes <- rownames(dea_res$Age)
    
    # Counts
    counts <- readRDS(paste0("~/Documents/mn4/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
    # counts <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.counts.PC_lincRNA.rds"))
    
    # Subset samples with disease annotation data
    colnames(counts) <- sapply(colnames(counts), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
    if(nchar(metadata$Sample[1])>18){ #The donor disease have the whole sample name.
      metadata$Sample <- sapply(metadata$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
    }
    
    counts <- counts[,colnames(counts) %in% metadata$Sample]
    
    
    # Create DGEList object
    dge <- DGEList(counts[exprs_genes,])
    
    # Calculate normalization factors (does not do the normalization, only computes the factors)
    dge <- calcNormFactors(dge)
    
    # Voom
    v <- voom(dge, design = NULL, normalize="quantile", save.plot=F, plot = F) # samples are treated as replicates
    
    covariates <- c("HardyScale","IschemicTime", "RIN", "Cohort", "NucAcIsoBatch", "ExonicRate", "PEER1", "PEER2")
    covariates <- covariates[covariates %in% names(metadata)]
    
    
    #Individual traits
    individual_traits <- names(metadata)[!names(metadata) %in% covariates]
    individual_traits <- individual_traits[!individual_traits %in% c("Donor", "Sample")]
    if("Atherosclerotic_Arteries" %in% individual_traits){ #Metadata with too many variables
      individual_traits <- individual_traits[!individual_traits %in% c("Atherosclerosis", "Atherosis", "Calcification", "Sclerotic")]
    }
    
    # Expression ~ covariates 
    fml_args_mod <- paste(c(covariates), collapse = " + ")
    mod <- model.matrix( as.formula(paste(" ~  ", paste(fml_args_mod,collapse = " "))), data =  metadata)
    
    # Limma fit
    fit <- lmFit(v, mod)
    
    # Calculate expression residuals
    exprs_residuals <- resid(fit, v)
    
    saveRDS(exprs_residuals, paste0(outdir, "/", tissue, ".", name, ".expression_residuals.rds"))

    # Differentially expressed genes
    traits <- names(dea_res) #In the model for T1D we only want the DEGs with T1D, so it is surprisingly convinient
    
    de_genes <- unique(unlist(lapply(traits, function(trait) rownames(dea_res[[trait]][dea_res[[trait]][,"adj.P.Val"] < 0.05,]) ) ))
    saveRDS(de_genes, paste0(outdir, "/", tissue, ".", name, ".genes_DE.rds"))
    
    # DE residuals
    DE_residuals <- exprs_residuals[de_genes,]
    
    # hier.part
    print("Hier Part")
    hier.part <- lapply(rownames(DE_residuals), function(gene)
      hier.part.mod(y=DE_residuals[gene,], x=metadata[,individual_traits], fam = "gaussian", gof = "Rsqu"))
    names(hier.part) <- rownames(DE_residuals)
    
    # Parse information
    rsq <- sapply(rownames(DE_residuals), function(gene) 
      sum(hier.part[[gene]]$IJ[,1])) 
    names(rsq) <- rownames(DE_residuals)
    
    rel_perc <- do.call(rbind.data.frame,
                        lapply(names(hier.part), function(gene) 
                          as.numeric(unlist(hier.part[[gene]]$I.perc))))
    rownames(rel_perc) <- rownames(DE_residuals)
    colnames(rel_perc) <- individual_traits
    
    abs_perc <-  do.call(rbind.data.frame,
                         lapply(names(hier.part), function(gene) 
                           hier.part[[gene]]$IJ[,1]*100)
    )
    rownames(abs_perc) <- rownames(DE_residuals)
    colnames(abs_perc) <- individual_traits
    
    hier_data <- cbind.data.frame(rsq,rel_perc,abs_perc)
    colnames(hier_data) <- c("R2",paste0(individual_traits,"_rel"),paste0(individual_traits,"_abs"))
    
    saveRDS(hier_data, paste0(outdir, "/", tissue, ".", name, ".hier_part.PEER.rds"))
    
    print("saved")
  }
    
}


# # -------------- #
# print(Sys.time())
# -------------- #