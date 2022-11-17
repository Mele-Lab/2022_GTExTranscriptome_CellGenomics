# Options ####
options(stringsAsFactors = F)

# Libraries ####
library(optparse)

# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-d", "--dir"), type="character",
                     dest="dir",
                     help="Folder where all data is stored")

options=parse_args(parser)
dir=options$dir

suffix <- ".voom_limma.results_PEER.rds"


#dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"
# dir <- "~/Documents/mn4/"

diseases_info <- read.delim(paste0(dir, "/Jose/00_Data/DonorDisease.Minimum_6_Donors.txt"))
diseases <- read.delim(paste0(dir, "Jose/00_Data/diseases.txt"))

# Tissues
tissue_info <- readRDS(paste0(dir, "/Jose/00_Data/Tissue_info.rds"))
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissues <- tissue_info$tissue_ID

# Number of affected samples ####
affected_samples <- list()
total_samples <- list()
for(dd in diseases[,1]){ ##dd is the disease name (e.g., Chronic_Respiratory_Disease), and disease will be the acronym MHCOPD
  print(dd)
  if(dd %in% diseases_info[,2]){ #If donor disease
    disease <- diseases_info[diseases_info$Description==dd, "Acronym"]
  } else{
    evaluation <- any(sapply(tissues, function(tissue) grepl(tissue, dd)))
    if(evaluation){ #If sample disease of the type Edinema_Thyroid
      disease <- sub("_[^_]+$", "", dd)
      if(grepl("BreastMammaryTissue", disease)){
        disease <- sub("_[^_]+$", "", disease)
      }
    } else{disease <- dd} #Smoking and Atherosclerotic_Arteries
  }
 
  v <- list()
  t <- list()
  outdir <- paste0(dir, "/Jose/01_Overview/Final_Diseases/")
  for(tissue in tissues){
    if(file.exists(paste0(outdir, dd, "/Tissues/", tissue, "/", 
                          tissue, ".SampleMetadata.", dd, ".rds"))){
      mdata <- readRDS(paste0(outdir, dd, "/Tissues/", tissue, "/", 
                              tissue, ".SampleMetadata.", dd, ".rds"))
    
      # 0s and 1s
      d <- sapply(c("0","1"), function(i) sum(mdata[,disease]==i))
      v[[tissue]] <- d["1"]
      t[[tissue]] <- d["1"] + d["0"]
        
    }else{
      v[[tissue]] <- NA
      t[[tissue]] <- NA
    }
  }
  v <- unlist(v)
  t <- unlist(t)
  names(v) <- tissues
  names(t) <- tissues
  affected_samples[[dd]] <- v
  total_samples[[dd]] <- t
}

number_disease_samples <- do.call(cbind.data.frame, affected_samples)
number_total_samples <- do.call(cbind.data.frame, total_samples)

trues <- number_disease_samples<6 #New threshold  
number_disease_samples[trues] <- NA
number_total_samples[trues] <- NA

# Number of DEG ####
deg <- list()
# DEG ####
for(dd in diseases[,1]){
  print(dd)
  og <- dd
  #Getting acronyms
  if(dd %in% diseases_info[,2]){ #If donor disease
    disease <- diseases_info[diseases_info$Description==dd, "Acronym"]
  } else{
    evaluation <- any(sapply(tissues, function(tissue) grepl(tissue, dd)))
    if(evaluation){ #If sample disease of the type Adenoma_Thyroid
      disease <- sub("_[^_]+$", "", dd)
      if(grepl("BreastMammaryTissue", disease)){
        disease <- sub("_[^_]+$", "", disease)
      }
    } else{disease <- dd} #Atherosclerotic_Arteries
  }
  
  v <- list()
  for(tissue in tissues){
    if(file.exists(paste0(dir, "Jose/01_Overview/Final_Diseases/", dd, "/Tissues/",tissue,"/", 
                          tissue,".", dd, suffix))){
      # DEA res
      dea_res <- readRDS(paste0(dir, "Jose/01_Overview/Final_Diseases/", dd,"/Tissues/",tissue,"/", 
                                tissue,".", dd, suffix))
      
      if(is.null(nrow(dea_res[[paste0(disease, 1)]]))){  #Because of congestion_spleen
        v[[tissue]] <- NA
        next
      }
      n <- sum(dea_res[[paste0(disease, 1)]]$adj.P.Val < 0.05, na.rm=TRUE)
      v[[tissue]] <- n
    }else{
      v[[tissue]] <- NA
    }
  }
  v <- unlist(v)
  names(v) <- tissues
  deg[[dd]] <- v
}

print("done part 2")

number_deg <- do.call(cbind.data.frame, deg)
print("These NAs before removing info:")
print(sum(is.na(number_deg)))
number_deg[is.na(number_disease_samples)] <- NA
print("These NAs after removing info:")
print(sum(is.na(number_deg)))

save.image(file = paste0(dir, "/Jose/01_Overview/Final_Diseases/heatmap_data.RData"))
