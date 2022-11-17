# Set options and load libraries ####
options(stringsAsFactors = F)
library(optparse)

#Reading arguments
# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-a", "--annotation_path"), type="character",
                     dest="annotation_path",
                     help="Folder containing the donor annotations")
parser <- add_option(parser, opt_str=c("-t", "--tissue_file"), type="character",
                     dest="tissue_file",
                     help="File including the information about tissues")
parser <- add_option(parser, opt_str=c("-d", "--metadata_dir"), type="character",
                     dest="dir",
                     help="Folder containing the metadata per tissue")
options=parse_args(parser)

annotation_path=options$annotation_path 
tissue_file=options$tissue 
dir=options$dir 

# Reading data ####
DonorInfo <- read.delim(paste0(annotation_path, "GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"),header =T, sep=",") #Data
DonorPhenotypes <- read.csv(paste0(annotation_path, "GTEx.SubjectPhenotypesDD.v8.csv"), header = T) #Metadata

tissue_info <- readRDS(tissue_file)
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissues <- tissue_info$tissue_ID

metadata <- lapply(tissues, function(tissue) #This reads the metadata generated in previous steps of the pipeline
  readRDS(paste0(dir, "/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds")))
names(metadata) <- tissues

donors <- unique(unlist(sapply(tissues, function(tissue) metadata[[tissue]]$Donor))) #Donors with gene expression data
DonorInfo <- DonorInfo[DonorInfo$SUBJID %in% donors,] #From 980 to 781 donors 

# Related to disease Autoimmune, Degenerative, Neurological and General Medical History
AutoimmuneDegenerativeNeurological <- unique(DonorPhenotypes[DonorPhenotypes$DOCFILE=="Autoimmune, Degenerative, Neurological","VARNAME"])
names(AutoimmuneDegenerativeNeurological) <- unique(DonorPhenotypes[DonorPhenotypes$DOCFILE=="Autoimmune, Degenerative, Neurological","VARDESC"])

GeneralMedicalHistory <- unique(DonorPhenotypes[DonorPhenotypes$DOCFILE=="General Medical History","VARNAME"])
names(GeneralMedicalHistory) <- unique(DonorPhenotypes[DonorPhenotypes$DOCFILE=="General Medical History","VARDESC"])
DonorDiseasePhenotypes <- c(AutoimmuneDegenerativeNeurological,GeneralMedicalHistory)

counts <- sapply(DonorDiseasePhenotypes, function(i) table(DonorInfo[,i]))
names(counts) <- DonorDiseasePhenotypes

# Selecting phenotypes with more than 5 donors with the disease
threshold <- 5
min0 <- sapply(1:length(counts), function(i) ifelse(counts[[i]]["0"]>threshold,i,NA)) # All have at least 10 donors annotated as 0
min1 <- sapply(1:length(counts), function(i) ifelse(counts[[i]]["1"]>threshold,i,NA)) # Not all have at least 10 donors annotated as 1

SelectedDonorDiseasePhenotypes <- DonorDiseasePhenotypes[which(!is.na(min1))]
DonorDiseasePhenotypesData <- cbind.data.frame(SelectedDonorDiseasePhenotypes,
                                               gsub(",","",gsub("'","" ,gsub("_$","",gsub(" ","_",sapply(names(SelectedDonorDiseasePhenotypes), function(i) unlist(strsplit(i,split = "\\("))[[1]]  )),perl = T))),
                                               names(SelectedDonorDiseasePhenotypes))
rownames(DonorDiseasePhenotypesData) <- NULL
colnames(DonorDiseasePhenotypesData) <- c("Acronym","Description","Extended_Description")
DonorDiseasePhenotypesData$Description[DonorDiseasePhenotypesData$Description=="Diabetes_mellitus_type_II"] <- "Diabetes_mellitus_type_2"
DonorDiseasePhenotypesData$Description[DonorDiseasePhenotypesData$Description=="Heart_attack_acute_myocardial_infarction_acute_coronary_syndrome"] <- "Heart_attack"

#This file will contain the diseases for downstream analysis
write.table(DonorDiseasePhenotypesData,
            paste0(dir, "/Jose/00_Data/DonorDisease.Minimum_6_Donors.txt"),
            col.names = T, row.names = F,
            quote = F, sep = "\t")


#Part 2: Preparing Metadata
whole_list_of_diseases <- c()

for(disease in DonorDiseasePhenotypesData[,1]){
  disease_name <- DonorDiseasePhenotypesData[DonorDiseasePhenotypesData$Acronym==disease, "Description"]
  print(disease_name)
  whole_list_of_diseases <- c(whole_list_of_diseases, disease_name)
  #Saving path
  outpath <- paste0(dir, "/Jose/01_Overview/Final_Diseases/", disease_name,"/")
  dir.create(outpath, recursive = T, showWarnings = FALSE)
  
  # Phenotype of interest
  DonorInfo_disease <- DonorInfo[,c("SUBJID",disease)]
  Donors_to_exclude <- c()
  
  #Removing missanotations from Diabetes (Individuals annotated as type 1 expression normal insulin levels)
  if(disease=="MHT1D"){
    DonorInfo_tmp <- DonorInfo_disease
    #Pancreas tpms to exclude individuals expressing insulin
    tpms <- readRDS(paste0(dir, "Raquel/Draft/Data/Tissues/Pancreas/Pancreas.SelectedSamples.TPM.PC_lincRNA.rds"))
    
    #Checking expression of insulin
    tpms <- data.frame(Sample=colnames(tpms), TPM=tpms["ENSG00000254647.6",]) # ENSG00000254647.6 is Insulin
    tpms$SUBJID <-  sapply(tpms$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:2],collapse="-" ) ) 

    tpms <- tpms[,c(3,2)]
    DonorInfo_tmp <- merge(DonorInfo_tmp, tpms, by="SUBJID")
    
    #New_annotation 
    Donors_to_exclude <- c(Donors_to_exclude, DonorInfo_tmp[(DonorInfo_tmp$MHT1D==1 & log10(DonorInfo_tmp$TPM+1)>2),"SUBJID"]) #Excluding diabetics with too much INS expression (6 donors)
    Donors_to_exclude <- c(Donors_to_exclude, DonorInfo_tmp[(DonorInfo_tmp$MHT1D==0 & log10(DonorInfo_tmp$TPM+1)<2),"SUBJID"]) #Excluding non-diabetics with poor INS expression (1 donor)
    Donors_to_exclude <- unique(Donors_to_exclude)
  }
  DonorInfo_disease <- DonorInfo_disease[!DonorInfo_disease$SUBJID %in% Donors_to_exclude,]
  
  # Add information regarding phenotype of interest to tissue metadata
  counts <- list()
  for(tissue in tissues){
    # Original number of samples
    n0 <- nrow(metadata[[tissue]])
    # Samples with metadata
    d <- merge(metadata[[tissue]], DonorInfo_disease, by.x = "Donor", by.y = "SUBJID")
    n1 <- nrow(d)
    # Samples that are either 0 o 1
    d <- d[d[,disease]!=99,]
    d[,disease] <- as.factor(d[,disease])
    n2 <- nrow(d)
    number_99s <- n1-n2
    # Count 0s and 1
    number_0s <- sum(d[,disease]=="0")
    number_1s <- sum(d[,disease]=="1")
    
    
    c <- c(n0, n1, number_99s, number_0s, number_1s, n2)
    counts[[tissue]] <- c
    
    # Save new metadata  
    if(number_1s>threshold){
      print(tissue)
      out <- paste0(outpath, "/Tissues/",tissue, '/')
      dir.create(out, recursive = T, showWarnings = FALSE)
      saveRDS(d, paste0(out,tissue,".SampleMetadata.",disease_name,".rds"))
    }
    
  }
  
  count_data <- do.call(rbind.data.frame,
                        counts)
  colnames(count_data) <- c("N","N_annotated","99","0","1","real_N")
  rownames(count_data) <- sapply(tissues, function(tissue) tissue_info[tissue_info$tissue_ID==tissue,]$tissue_abbrv)
  count_data <- count_data[order(rownames(count_data)),]
  # Save data

  # Save summary of samples with disease
  write.table(count_data,
              paste0(outpath,"Samples_with_", disease_name,"_metadata.txt"),
              col.names = T, row.names = T,
              quote = F, sep = "\t")

}


whole_list_of_diseases <- c(whole_list_of_diseases)
write.table(whole_list_of_diseases, paste0(dir, "Jose/00_Data/donor_diseases.txt"),
            row.names = F, col.names=F, sep = "\t", quote = F)

# write.table(tissue_info$tissue_ID, "~/Documents/mn4/Jose/00_Data/Tissues.txt", quote=F, row.names=F, col.names=F)
