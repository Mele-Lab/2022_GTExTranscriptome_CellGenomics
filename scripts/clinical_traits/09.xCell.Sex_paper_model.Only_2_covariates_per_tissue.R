#!/usr/bin/env Rscript
#Sin filtro en general con los que están validados

# Libraries ----
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(hier.part)
library(yarrr)

first_dir <- "~/Documents/mn4/"

source(paste0(first_dir, "/Raquel/R_functions/DEA_and_DSA.R_functions.R"))

main_data <- read.csv(paste0(first_dir, "Jose/00_Data/main_final.csv"))
acronyms <- main_data[,1]
tissues <- main_data[,2]
folders <- tissues
names <- paste0(acronyms, "_", tissues)

# Traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" =  "#009E73",
                 "BMI" = "#CC79A7",
                 "Disease" = "#000000")
traits <- names(traits_cols)

# Metadata ----
metadata_list <- list()
for(i in 1:length(folders)){
  folder <- folders[i]
  tissue <- tissues[i]
  acronym <- acronyms[i]

  outdir <- paste0(first_dir, "Jose/03_Models/Tissues/", folder, "/")
  file <- list.files(outdir, full.names = T)[grepl("SampleMetadata",list.files(outdir))]
  metadata <- readRDS(file)
  if(tissue=="Testis"){
    metadata$Sex <- "1"
  }
  metadata$Tissue <- tissue
  name <- paste0(acronym, "_", tissue)
  metadata$name <- name  #Should I keep this?
  # metadata <- metadata[, !colnames(metadata) %in% c("HardyScale","Cohort", "NucAcIsoBatch","ExonicRate",
  #                           "PEER1", "PEER2")] #Why these? Ask Raquel, and she was also excluding RIN and Ischemic Time even though she told me not to.
  metadata <- metadata[, !colnames(metadata) %in% c("PEER1")] #Why these? Ask Raquel, and she was also excluding RIN and Ischemic Time even though she told me not to.
  metadata_list[[name]] <- metadata
}

full.metadata <- metadata_list

# xCell scores ----
xCell <- readRDS(paste0(first_dir, "Jose/00_Data/xCell.scored_adjusted.rds"))
colnames(xCell) <- sapply(colnames(xCell), function(i) paste(unlist(strsplit(i,split = "-"))[1:3],collapse="-" ) )
xCell.tissues <- lapply(names, function(name) xCell[, full.metadata[[name]]$Sample])
names(xCell.tissues) <- names



# Benchmarked cell types ----
cell_types <- c("Neurons", 
                "Epithelial cells",
                "Keratinocytes",
                "Adipocytes",
                "Myocytes",
                "Hepatocytes",
                "Neutrophils")


macrophages <- F
if(macrophages==T){
  cell_types <- c("Macrophages")
  names <- c("Atelectasis_Lung", "Emphysema_Lung", "Fibrosis_Lung", "Pneumonia_Hist_Lung")
}

# Select only tissues with benchmarked cell types ----
#Not using the threshold now: Using a smaller onw now
# tissue.cellTypes <- lapply(names, function(name) cell_types[apply(xCell.tissues[[name]][cell_types,], 1, function(x) median(x)) > 0.1])
tissue.cellTypes <- lapply(names, function(name) cell_types[apply(xCell.tissues[[name]][cell_types,], 1, function(x) median(x)) > 0.01])
names(tissue.cellTypes) <- names

if(macrophages==T){
  tissue.cellTypes <- lapply(names, function(name) cell_types[TRUE]) #The expression is very low, but let's see the plot
  names(tissue.cellTypes) <- names
}
#NerveTibial in both Diabetes have 0.01 as median in neurons, while epithelial cells have 0.049
tissue.cellTypes <- tissue.cellTypes[sapply(names, function(name) length(tissue.cellTypes[[name]]) > 0)]
tissues <- names(tissue.cellTypes) #Testis and NerveTibial are gone   

# Inverse normal transformed xCell scores ----
i.xCell <-  t(apply(xCell,1, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))) )
i.xCell.tissues <- lapply(names, function(name) i.xCell[, full.metadata[[name]]$Sample])
names(i.xCell.tissues) <- names

parse.data <- function(tissue){
  d <- as.data.frame(i.xCell.tissues[[tissue]])
  d$CellType <- rownames(d)
  m <- melt(d)
  colnames(m) <- c("CellType","Sample","ixCell")
  m <- m[m$CellType %in% tissue.cellTypes[[tissue]],]
  m$Sample <- as.character(m$Sample)
  df <- as.data.frame(t(sapply(full.metadata[[tissue]]$Sample, function(sample) 
    sapply(tissue.cellTypes[[tissue]], function(i)
      m[m$Sample==sample & m$CellType==i,"ixCell"]
    )
  )
  ))
  
  if(length(tissue.cellTypes[[tissue]])==1){
    df <- cbind.data.frame(unlist(df), gsub(paste0(".",tissue.cellTypes[[tissue]]) , "", names(df)))
    colnames(df) <- c(tissue.cellTypes[[tissue]], "Sample")
    colnames(df) <- gsub(" ","_", colnames(df))
  }else{
    df$Sample <- rownames(df)
    colnames(df) <- gsub(" ","_", colnames(df))
  }
  
  return(df)
}

ixCell <- lapply(names, function(name) parse.data(name))
names(ixCell) <- names

# Merge with tissue metadata ----
full.data <- lapply(names, function(name) 
  merge(full.metadata[[name]],
        ixCell[[name]],
        by = "Sample")
)
names(full.data) <- names

# Analysis ----
# Analysis function --

hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
  pcan <- dim(xcan)[2]
  gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
  hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
  
  params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                 family = family, link = link, gof = gof)
  if(TRUE %in% (hp$IJ$I<0)){
    hp$IJ$I[hp$IJ$I<0] <- NA
    list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
  }else{
    list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
  }
}

compositional.changes.fun <- function(tissue){
  print(tissue)
  
  data <- full.data[[tissue]] #I already included RIN and IschemicTime since the beggining

  # mdata <- full.metadata[[tissue]]
  # covariates <- colnames(mdata)[!colnames(mdata) %in% colnames(data)]
  # covariates <- covariates[covariates %in% c("IschemicTime",  "RIN")]
  # data <- merge(data, mdata[,c("Sample", covariates)], by="Sample")
  
  # Cell types in tissue
  # cell_types <- colnames(data)[!colnames(data) %in% c("Donor", "Sample", "Age", "Ancestry", "Sex", "BMI", "Tissue", covariates)]
  covariates_og <- colnames(data)[colnames(data) %in% c("IschemicTime", "RIN", "HardyScale", "Cohort", "NucAcIsoBatch", "ExonicRate", "PEER2")]
  covariates <- covariates_og[covariates_og %in% c("IschemicTime",  "RIN")]
  cell_types <- colnames(data)[!colnames(data) %in% c(covariates_og, "Donor", "Sample", "Age", "Ancestry", "Sex", "BMI", "Tissue", "name", "Emphysema", acronyms)] #Emphysema is the only disease not in acronyms
  
  # Traits
  traits <- names(data)[names(data) %in% c("Age", "Ancestry", "Sex", "BMI", acronyms)]
  if(grepl("Testis", tissue)){
    traits <- traits[-grep("Sex", traits)]
    data <- data[,-grep("Sex", colnames(data))]
  }
  traits.names <- c()
  for(trait in traits){
    if(trait=="Age" | trait == "BMI"){
      traits.names <- c(traits.names, trait)
    } else if(trait=="Ancestry"){
      traits.names <- c(traits.names, "AncestryAFR")
    } else if(trait=="Sex"){
      traits.names <- c(traits.names, "Sex2")
    } else{ #Assuming the rest are diseases
      data[[trait]] <- as.factor(data[[trait]])
      traits.names <- c(traits.names, paste0(trait, 1))
    }
  }
  
  # cell_type ~ Age + Ancestry + Sex + BMI
  lm.obj <- lapply(cell_types, function(cell_type) 
    lm(as.formula(paste(cell_type," ~  ", paste(paste(c(covariates,traits), collapse = " + "), collapse = " "))),
       data)
  )
  names(lm.obj) <- cell_types
  
  # P.value matrix: rows=traits; columns=cell_types
  m <- sapply(cell_types, function(cell_type) summary(lm.obj[[cell_type]])$coefficients[,4])
  m <- m[-1,]
  if(length(cell_types) > 1){
    m <- m[traits.names,]
    M <- t(m)#t(apply(m, 1, function(x) p.adjust(x, method = "BH")))  
  }else{
    M <- m[traits.names]
  }
  
  # coefficients
  coeff <- sapply(cell_types, function(cell_type) summary(lm.obj[[cell_type]])$coefficients[,1])
  coeff <- coeff[-1,]
  if(length(cell_types) > 1){
    coeff <- coeff[traits.names,]
  }else{
    coeff <- coeff[traits.names]
  }
  # R2
  R2.values <- lapply(cell_types, function(cell_type) hier.part.mod(y=data[,cell_type], 
                                                                    x = data[,traits], 
                                                                    fam = "gaussian", 
                                                                    gof = "Rsqu")$IJ[,1]*100 )
  names(R2.values) <- cell_types
  for(cell_type in cell_types){
    names(R2.values[[cell_type]]) <- traits
  }
  R2 <- as.matrix(do.call(cbind.data.frame, R2.values))
  #r2 <- R2
  #R2[is.na(M)] <- NA
  
  return(list("p.value" = M,
              "Coefficients" = coeff,
              "R2" = R2))
}

# Results --
final_names <- unique(names(tissue.cellTypes))
results <- lapply(final_names, function(name) compositional.changes.fun(name))
names(results) <- final_names

# Parse data 4 dotplot ----
# Tissue CellType Trait fdr R2(Size) Coeff | FDR
get.tissue.table <- function(tissue){
  print(tissue)
  if(length(tissue.cellTypes[[tissue]])==1){
    df <- cbind.data.frame(rep(tissue, length(results[[tissue]]$p.value)),
                           gsub("1", "", gsub("2", "", gsub("AFR","", names(results[[tissue]]$p.value)))),
                           rep(gsub(" ", "_", tissue.cellTypes[[tissue]]),length(results[[tissue]]$p.value)),
                           results[[tissue]]$p.value,
                           results[[tissue]]$Coefficients,
                           results[[tissue]]$R2)
    colnames(df) <- c("Tissue","Trait","CellType","p.value","Coeff","R2")
  }else{
    df <- cbind.data.frame(melt(t(results[[tissue]]$p.value)),
                           melt(results[[tissue]]$Coefficients),
                           melt(results[[tissue]]$R2)
    )
    colnames(df) <- c("Var1","CellType","p.value","Var2","CellType2","Coeff","Trait","CellType2","R2")
    df$Tissue <- rep(tissue, nrow(df))
    df <- df[,c("Tissue","Trait","CellType","p.value","Coeff","R2")]
  }
  #df$Significant <- as.factor(ifelse(df$fdr < 0.05, "Yes","No"))
  df$P.value <- ifelse(df$Coeff > 0, -log10(df$p.value), log10(df$p.value))
  return(df)
}

data <- do.call(rbind.data.frame, 
                lapply(final_names, function(tissue) get.tissue.table(tissue)))

data_final <- data[1,]
for(row in 1:nrow(data)){
  if(grepl(data[row,2],data[row,1])){
    data_final <- rbind(data_final, data[row,])
  }
}
data_final <- data_final[-1,]

# data_final <- data[c(5,10,15,20,25,30,35,40,45,50,55,60,61,67:71,76:77,82,83,88,93,98,103,108,113,118,123,128),]
# data_final <- data[c(5,15,24,31,46,55,63,71,81,89,96,103,111,121,131,138,146,154,162,170,176,182,190,196,202),]

if(macrophages==T){
  data_final <- data[c(5,16,27,38),]
}

data_final$adj.P.Val <- p.adjust(data_final$p.value, method = "BH")
data_final$FDR <- ifelse(data_final$Coeff > 0, -log10(data_final$adj.P.Val), log10(data_final$adj.P.Val))


# Data
# plot.data <- do.call(rbind.data.frame, data_final)
plot.data <- data_final
plot.data$Significant <- as.factor(ifelse(plot.data$adj.P.Val < 0.05, "Yes","No"))


plot.data$TissueCell <- paste0(plot.data$Tissue, ":", plot.data$CellType)
plot.data$FakeFDR <- ifelse(plot.data$FDR > 3, 3,
                            ifelse(plot.data$FDR < -3, -3,
                                   plot.data$FDR))
# plot.data$TissueAbbrv <- sapply(plot.data$Tissue, function(i) tissue_info[tissue_info$tissue_ID==i,"tissue_abbrv"])
# tissue.order <- c(grep("BRN", tissue_info$tissue_abbrv, value = T),
#                   "PTTARY",
#                   "THYROID",
#                   "LUNG", 
#                   "PNCREAS",
#                   "STMACH",
#                   "CLNTRN" ,
#                   "SNTTRM",
#                   "PRSTTE" ,
#                   "SLVRYG" ,
#                   "ESPMCS",
#                   "VAGINA",
#                   "SKINNS",
#                   "SKINS",
#                   "LIVER",
#                   "BREAST",
#                   #"BREAST_FEMALE",
#                   #"BREAST_MALE",
#                   "ADPVSC",
#                   "ADPSBQ",
#                   "HRTAA",
#                   "HRTLV",
#                   "MSCLSK",
#                   "SPLEEN",
#                   "WHLBLD"
# )
# plot.data$TissueAbbrv <- factor(plot.data$TissueAbbrv, levels = rev(tissue.order), order = T)
plot.data$CellType <- factor(plot.data$CellType, levels = c("Neurons", 
                                                            "Epithelial_cells",
                                                            "Keratinocytes",
                                                            "Hepatocytes",
                                                            "Adipocytes",
                                                            "Myocytes",
                                                            "Neutrophils"))
if(macrophages ==T){
  plot.data$CellType <- factor(plot.data$CellType, levels="Macrophages")
  plot.data$Tissue <- c("Atelectasis: LUNG", "Emphysema: LUNG", "Fibrosis: LUNG", "Pneumonia: LUNG")
}

# plot.data$Tissue <- c(
#   "Atelectasis: LUNG",
#   "Atrophy: THYROID",
#   "Atrophy: MSCLSK",
#   "Congestion: ESPMCS", 
#   "Emphysema: LUNG",
#   "Fibrosis: LIVER",
#   "Fibrosis: LIVER",
#   "Fibrosis: LIVER",
#   "Fibrosis: LUNG",
#   "Gastritis: STMACH", 
#   "Gynecomastoid: BREAST",
#   "Gynecomastoid: BREAST",
#   "Hashimoto: THYROID",
#   "Hyperplasia: THYROID",
#   "Pneumonia: LUNG",
#   "Saponification: PNCREAS",
#   "Steatosis: LIVER",
#   "Steatosis: LIVER",
#   "Steatosis: LIVER",
#   "Type 1 diabetes: PNCREAS",
#   "Type 1 diabetes: ADPSBQ",
#   "Type 1 diabetes: ADPVSC",
#   "Type 2 diabetes: PNCREAS",
#   "Type 2 diabetes: ADPSBQ",
#   "Type 2 diabetes: ADPVSC"
# )

# plot.data$Trait <- factor(plot.data$Trait, levels = c("Age","Ancestry","Sex","BMI"),order = T)
# plot.data$Tissue <- factor(plot.data$Tissue, levels = unique(plot.data$Tissue))
# Plot ----
plot.data <- plot.data[!is.na(plot.data$CellType),]
g <- ggplot(plot.data) +
  geom_point(aes(x=CellType, y=Tissue, size = R2, col = FakeFDR, shape = Significant)) + 
  theme_minimal() +
  xlab("") + ylab("") +
  scale_color_gradient2(low = brewer.pal(11,"RdBu")[10],
                        mid = "white",
                        high = brewer.pal(11,"RdBu")[2]) +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1)) + #,
        #legend.position = "bottom") +
  scale_shape_manual(values =c(18,19)) +
  labs(colour = "-log10(FDR)",
       shape = "FDR < 0.05",
       size = "Variance explained (%)") #+
if(macrophages==T){
  ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell_Macrophages.png"), g, device="png",
         width = 4, height = 5)
}
ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell_0.01.png"), g, device="png",
       width = 6, height = 8)

ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell.png"), g, device="png",
       width = 6, height = 8)
ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell.pdf"), g, device="pdf",
       width = 6, height = 8)

ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell_all_cov.png"), g, device="png",
       width = 6, height = 8)
ggsave(paste0(first_dir, "Jose/03_Models/xCell/TissueComposition_xCell_all_cov.pdf"), g, device="pdf",
       width = 6, height = 8)
