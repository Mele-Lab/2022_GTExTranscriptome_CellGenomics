# Options
options(stringsAsFactors = F)
library(purrr) #for function is_empty

# Command line arguments
args = commandArgs(trailingOnly=TRUE)

out <- "Tissues"    #folder_name

# dir <- "/gpfs/projects/bsc83/Projects/GTEx_v8/"
dir <- "~/Documents/mn4/"

#Retrieving acronyms:
diseases_info <- read.delim(paste0(dir, "Jose/00_Data/DonorDisease.Minimum_6_Donors.txt"), header = T)


#Diseases with DEG in the given tissue:
load(paste0(dir, "Jose/01_Overview/Final_Diseases/heatmap_data.RData")) 
dir <- "~/Documents/mn4/"

data <- number_deg[,!colSums(number_deg, na.rm = T)==0]
#Excluding the histology phenotypes that are not diseases or that we are not interested in
data <- data[,!colnames(data) %in% c("Atherosis_ArteryTibial", "Calcification_ArteryAorta", "Calcification_ArteryTibial", "Macrophages_Lung", "Monckeberg_ArteryTibial", "Pigment_Lung", "Sclerotic_ArteryTibial", "Atherosclerotic_Arteries_Final")]
#Excluding the record diseases with no signal in the tissue of origin or diseases not present at time of death
data <- data[,!colnames(data) %in% c("Alzheimers_OR_Dementia", "Arthritis", "Asthma", "Major_depression", "Liver_Disease", "Renal_Failure", "Schizophrenia", "Chronic_Lower_Respiratory_Disease", "Multiple_Sclerosis", "Ischemic_Heart_Disease", "Rheumatoid_Arthritis", "Bacterial_Infections", "Heart_Disease", "Cerebrovascular_Disease", "Pneumonia", "Heart_attack")]
#Excluding COPD and the record diseases with less than 50
data <- data[,!colnames(data) %in% c("Smoking","Systemic_Lupus", "Chronic_Respiratory_Disease", "Hypertension", "Dementia_With_Unknown_Cause")]

tissue <- args[1] 
# tissue <- "NerveTibial"
print(tissue)

diseases_all <- data[tissue,]
diseases <- diseases_all[,!is.na(diseases_all)]   #What happens if there is only 1 disease?
if(length(diseases)==1){
  names(diseases) <- colnames(diseases_all)[!is.na(diseases_all)]
}
record_disease <- diseases[names(diseases) %in% c( "Diabetes_mellitus_type_1", "Diabetes_mellitus_type_2")]
record_disease <- names(record_disease)


histology_disease <- diseases[!colnames(diseases) %in% c("Dementia_With_Unknown_Cause","Systemic_Lupus", "Diabetes_mellitus_type_1", "Diabetes_mellitus_type_2", "Smoking")]
histology_disease <- colnames(histology_disease)[histology_disease>5] #Checking if the tissue of origin has more than 5
diseases <- c(histology_disease, record_disease)
if(length(diseases)==0){
  stop(paste0("There are no diseases affecting the ", tissue))
}

acronyms <- c()
name <- ""
for(disease in diseases){
  #This is needed just to go from disease name to acronym in the case of diseases like Hashimoto_Thyroid
  tissues <- list.dirs(paste0(dir, "Jose/01_Overview/Final_Diseases/", disease,"/Tissues/"), full.names = F)[-1] 
  evaluation <- any(sapply(tissues, function(tissue) grepl(tissue, disease)))
  if(name==""){name <- disease
  } else{name <- paste0(name, "_", disease)}
  if(evaluation){
    disease <- sub("_[^_]+$", "", disease)
    if(grepl("BreastMammaryTissue", disease)){
      disease <- sub("_[^_]+$", "", disease)
    }
    acronyms <- c(acronyms, disease)
  }else if(disease %in% diseases_info$Description){
    acronyms <- c(acronyms, diseases_info[diseases_info$Description==disease,"Acronym"])
  }else{acronyms <- c(acronyms, disease)}
}


#Reading data
metadata <- list()
for(i in 1:length(diseases)){
  dd <- diseases[i]
  metadata[[dd]] <- readRDS(paste0(dir, "Jose/01_Overview/Final_Diseases/", dd, "/Tissues/", tissue, "/", tissue, ".SampleMetadata.",
diseases[i], ".rds"))
  if(nchar(metadata[[dd]]$Sample[1])>18){ #The donor disease have the whole sample name, I shorten it
    metadata[[dd]]$Sample <- sapply(metadata[[dd]]$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
  }
  if("SAMPID" %in% colnames(metadata[[dd]])){
    metadata[[dd]] <- metadata[[dd]][,-1]}
}

merge_two <- function(x, y) {
    names <- intersect(colnames(x), colnames(y))
    merge(x, y, by=names)}

output <- Reduce(merge_two, metadata)

# Save data

outpath <- paste0(dir, "Jose/03_Models/", out, "/", tissue,"/")
dir.create(outpath, recursive = T, showWarnings = F)
final_name <- paste0(acronyms[!two_models], collapse="_")
saveRDS(output, paste0(outpath,tissue,".SampleMetadata.", final_name, ".rds"))

