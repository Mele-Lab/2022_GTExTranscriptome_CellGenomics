# Set options and load libraries ####
options(stringsAsFactors = F)
library(optparse)
suppressMessages(library(Hmisc)) #for function capitalize

#Reading arguments
# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-f", "--file"), type="character",
                     dest="file",
                     help="File downloaded from the GTEx Portal: histological data")
parser <- add_option(parser, opt_str=c("-t", "--tissue_file"), type="character",
                     dest="tissue_file",
                     help="File including the information about tissues")
parser <- add_option(parser, opt_str=c("-d", "--metadata_dir"), type="character",
                     dest="dir", default="/gpfs/projects/bsc83/Projects/GTEx_v8/",
                     help="Folder containing the metadata per tissue")
options=parse_args(parser)

file=options$file
tissue_file=options$tissue
dir=options$dir

# Data ####
# Metadata 
tissue_info <- readRDS(tissue_file)
tissue_info <- tissue_info[!grepl("BreastMammaryTissue_", tissue_info$tissue_ID),]
tissues <- tissue_info$tissue_ID

# Reading metadata per tissue
metadata <- lapply(tissues, function(tissue) #These files were output of previous steps of the pipeline
  readRDS(paste0(dir, "Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
)
names(metadata) <- tissues

# List all samples in our dataset with a shorter name (as in histological data)
samples <- unique(unlist(sapply(tissues, function(tissue) sapply(metadata[[tissue]]$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ))))  #unique because the shorter name is the same for the samples in some brain tissues (Only in brain), but that does not affect us

# Read histological data
final_histology_annotations <- read.csv(file)
  
# Keep the samples from the histological data that also have our metadata annotation
sample_phenotypes <- final_histology_annotations[final_histology_annotations$Tissue.Sample.ID %in% 
                                                   samples, ]

# Tissue nameing correspondance
tissue_names <- sort(unique(sample_phenotypes$Tissue))
actual_names <- sapply(tissue_names, function(name) tissue_info[tissue_info$tissue_name==name,"tissue_ID"], USE.NAMES = F)

#Printing number of samples per tissue in total, and with our annotation metadata (samples with expression data)
summary_samples_with_disease_phenotypes <- cbind.data.frame(tissue_names, actual_names,
                                                            sapply(1:length(tissue_names), function(i) 
                                                              nrow(final_histology_annotations[final_histology_annotations$Tissue==tissue_names[i],])
                                                            ),
                                                            sapply(1:length(tissue_names), function(i) 
                                                              nrow(sample_phenotypes[sample_phenotypes$Tissue==tissue_names[i],])
                                                            ))
colnames(summary_samples_with_disease_phenotypes) <- c("Tissue_name","Tissue_ID","N","n")

sample_phenotypes$Tissue <- sapply(sample_phenotypes$Tissue, function(i) actual_names[which(tissue_names==i)], USE.NAMES = F)

# Keeping disease phenotypes with more than 5 positives in total (not per tissue yet)
threshold <- 5
my_diseases <- colnames(sample_phenotypes)[10:ncol(sample_phenotypes)][sapply(colnames(sample_phenotypes)[10:ncol(sample_phenotypes)], 
function(i)
  sum(table(sample_phenotypes[sample_phenotypes[,i]==1,"Tissue"])>threshold)
  )>0]

sample_phenotypes <- sample_phenotypes[,c("Tissue.Sample.ID","Subject.ID","Tissue",my_diseases)]

#Capitalizing disease names
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}
my_disease_phenotypes <- sapply(my_diseases, CapStr)
colnames(sample_phenotypes)[4:ncol(sample_phenotypes)] <- my_disease_phenotypes


my_data <- t(sapply(sort(unique(sample_phenotypes$Tissue)), function(tissue) 
  sapply(my_disease_phenotypes, function(i)
    sum(sample_phenotypes[sample_phenotypes$Tissue==tissue,i]==1, na.rm=TRUE)
    )
  ))
colnames(my_data) <- sapply(colnames(my_data), CapStr, USE.NAMES = F)
my_data <- my_data[,-which(colnames(my_data)=="No_abnormalities" | colnames(my_data)=="Clean_specimens")]
write.table(my_data,
            paste0(dir, "Jose/00_Data/Sample_disease_phenotypes.txt"),
            col.names = T, row.names = T,
            quote = F, sep = "\t")


#Tissues with more than 5 positives per disease, the number just corresponds to the disease number
tissue_disease_combs <- apply(my_data, 2, function(x) which(x>threshold))
my_list <- vector()
for(i in 1:length(tissue_disease_combs)){
  x <- names(tissue_disease_combs[[i]])
  names(x) <- rep(names(tissue_disease_combs)[i],length(tissue_disease_combs[[i]]))
  my_list <- c(my_list,x)  
}

my_combs <- cbind.data.frame(my_list, names(my_list))
#Changing penumonia name to differentiate donor disease vs sample disease
my_combs[my_combs[,2]=="Pneumonia", 2] <- "Pneumonia_Hist"   
colnames(sample_phenotypes)[colnames(sample_phenotypes)=="Pneumonia"] <- "Pneumonia_Hist"


write.table(my_combs,
            paste0(dir, "Jose/00_Data/Sample_disease_phenotypes.Tissue_and_disease.txt"),
            col.names = F, row.names = F,
            quote = F, sep = "\t")



#Saving metadata for the disease subsets
whole_list_of_diseases <- c()

colnames(sample_phenotypes)[1]<- "Sample" #Changing variable name

# Add Sample ID to metadata and subset samples to those with histological annotation
for(disease in unique(my_combs[,2])){
  counts <- list()
  for(tissue in my_combs[my_combs[,2]==disease,1]){
    new_disease <- paste0(disease, "_", tissue)   #Fibrosis on testis is one phenotype and fibrosis in liver is another one
    print(new_disease)
    n0 <- nrow(metadata[[tissue]])
    if(is.null(metadata[[tissue]])){next}
    d <- metadata[[tissue]]

    d$Sample <- sapply(d$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #Changing sample name
    d <- merge(d, sample_phenotypes[,c("Sample",disease)], by="Sample") #Adding disease info to metadata
    d <- na.omit(d)
    if(disease=="Gynecomastoid"){
      d <- d[!(d$Gynecomastoid==1 & d$Sex==2),] #Removing females with gynecomastoid (2 individuals)
    }
    colnames(d) <- capitalize(colnames(d))
    whole_list_of_diseases <- c(whole_list_of_diseases, new_disease)

    #Counting
    n1 <- nrow(d)
    d <- d[d[,disease]!=99,]
    n2 <- nrow(d)
    number_99s <- n1-n2
    number_0s <- sum(d[,disease]=="0")
    number_1s <- sum(d[,disease]=="1")
    c <- c(tissue, n0, n1, number_99s, number_0s, number_1s, n2)
    counts[[tissue]] <- c
    if(number_1s>threshold){
      out <- paste0(dir, "Jose/01_Overview/Final_Diseases/", new_disease, "/Tissues/",tissue, '/')
      dir.create(out, recursive = T, showWarnings=FALSE)
      saveRDS(d, paste0(out,tissue,".SampleMetadata.", new_disease,".rds"))
    }

    #Saving count data
    count_data <- do.call(rbind.data.frame,
                          counts)
    rownames(count_data) <- count_data[,1]
    count_data <- count_data[,-1]
    colnames(count_data) <- c("N","N_annotated","99","0","1","real_N")
    count_data <- count_data[order(rownames(count_data)),]

    write.table(count_data,
                paste0(dir, "Jose/01_Overview/Final_Diseases/", new_disease,"/Samples_with_", new_disease,"_metadata.txt"),
                col.names = T, row.names = T,
                quote = F, sep = "\t")
  }
}




#Creating Atherosclerotic_Arteries, as a new variable combining four other similar phenotypes
disease <- "Atherosclerotic_Arteries"
artery_tissues <- grep("Artery",tissue_info$tissue_ID, value=T) #Three arteries
cvd_traits <-  c("Atherosclerosis","Atherosis", "Calcification","Sclerotic") #sclerotic includes Monckeberg
print(disease)

counts <- list()
DonorInfo <- c()
for(tissue in artery_tissues){ #For Arteries
  if(is.null(metadata[[tissue]])){next}
  d <- metadata[[tissue]]
  n0 <- nrow(d)
  d$Sample <- sapply(d$Sample, function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) ) #Changing Sample Names
  d <- merge(d, sample_phenotypes, by="Sample")
  
  d$Atherosclerotic_Arteries <- ifelse(apply(d[,cvd_traits], 1, function(x) sum(x))==0,0,1) #Creating new phenotype
  d <- d[,c(1:14, 20, 21, 23, 59, 64)] #I could also remove 19, 20, 22 and 56
  
  d <- na.omit(d)
  #Counting
  n1 <- nrow(d)
  d <- d[d[,disease]!=99,]
  n2 <- nrow(d)
  number_99s <- n1-n2
  number_0s <- sum(d[,disease]=="0")
  number_1s <- sum(d[,disease]=="1")
  c <- c(tissue, n0, n1, number_99s, number_0s, number_1s, n2)
  counts[[tissue]] <- c
  
  #Saving
  if(number_1s>threshold){
    out <- paste0(dir, "Jose/01_Overview/Final_Diseases/", disease, "/Tissues/",tissue, '/')
    dir.create(out, recursive = T, showWarnings=FALSE)
    saveRDS(d, paste0(out,tissue,".SampleMetadata.",disease,".rds"))
  }
  DonorInfo <- rbind(DonorInfo, d[,c("Donor", disease)]) #To keep track of the donors
}

#Saving count data
count_data <- do.call(rbind.data.frame, counts)
rownames(count_data) <- count_data[,1]
count_data <- count_data[,-1]
colnames(count_data) <- c("N","N_annotated","99","0","1","real_N")
count_data <- count_data[order(rownames(count_data)),]

write.table(count_data,
                paste0(dir, "Jose/01_Overview/Final_Diseases/", disease,"/Samples_with_", disease,"_metadata.txt"),
                col.names = T, row.names = T,
                quote = F, sep = "\t")

whole_list_of_diseases <- c(whole_list_of_diseases, disease)
write.table(whole_list_of_diseases, paste0(dir, "Jose/00_Data/sample_diseases.txt"),
            row.names = F, col.names=F, sep = "\t", quote = F)
 

