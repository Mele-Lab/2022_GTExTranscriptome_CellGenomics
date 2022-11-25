#!/usr/bin/env Rscript

dse_sharing <- readRDS("~/GTEx_v8/Raquel/Draft/Cell_submission/data/events_DS.tissue_sharing.rds")
dse_sharing$Ancestry[grep("CYP3A5",dse_sharing$Ancestry$gene_name),][1,]
head(dse_sharing$Ancestry[,1:8])

# Input data ----
args = commandArgs(trailingOnly=TRUE)
tissue <- args[1] #"ArteryAorta" 
trait <- args[2] # "Age"  
event_id <- args [3] # "ENSG00000115414.18;SE:chr2:215391814-215392931:215393203-215394528:-"
outpath <- args[4] # "~/GTEx_v8/Raquel/Draft/Analysis/splicing/figures/splicing_events_examples/" # 
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

# Libraries  ----
library(ggplot2)
library(ggpubr)
library(ggtext)

# Demographic traits ----
traits_cols <- c("Age" = "#56B4E9",
                 "Ancestry" = "#E69F00",
                 "Sex" = "#009E73",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Gene annotation ----
# PCG and lincRNA genes with matched biotype annotation in gencode v38
# PAR genes excluded
#gene_annotation <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
gene_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#length(unique(gene_annotation$ensembl.id)) # 26,196 genes

# Transcript annotation ----
transcript_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(transcript_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "gene_biotype", "transcript.name", "transcript_biotype")

# Event annotation ----
#event_annotation <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")
event_annotation <- readRDS("~/GTEx_v8/Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")

# event label
event_label <- paste0(unlist(strsplit(event_annotation[event_id, "event.label"], split = ":"))[1],
                     " - ",
                     unlist(strsplit(event_annotation[event_id, "event.label"], split = ":"))[2],
                     " (", unlist(strsplit(event_annotation[event_id, "event.label"], split = ":"))[3], ")")

# Metadata ----
#metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/", tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues//",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")

# Event PSI values ----
#psi <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
psi <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))

# Event TPM values ----
#event.tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))
#event.tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))

# Transcript TPM ----
transcript.tpm <- read.delim(paste0("~/GTEx_v8/Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue, ".transcript_TPM.txt"))
colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))
# Subset to tissue samples --
transcript.tpm <- transcript.tpm[, metadata$Sample]
identical(colnames(transcript.tpm), metadata$Sample)

# Transcript splicing ratio ----
transcript.sratio <- read.delim(paste0("~/GTEx_v8/Raquel/Draft/SUPPA/PSI_isoforms/", tissue, "/", tissue, "_isoform.psi"))
colnames(transcript.sratio) <- gsub("\\.", "-", colnames(transcript.sratio))
rownames(transcript.sratio) <- sapply(rownames(transcript.sratio), function(i) unlist(strsplit(i, split = ";"))[[2]])
# Subset to tissue samples --
transcript.sratio <- transcript.sratio[, metadata$Sample]
identical(colnames(transcript.sratio), metadata$Sample)

# 1. Plot splicing event PSI ----
get_box_stats <- function(y){
  if(max(y) > 0.9){
    lower_limit = min(y) * 0.85
    return(data.frame(
      y = 0.95 * lower_limit,
      label = paste(
        "N =", length(y), "\n"#,
        #"Mean =", round(mean(y), 2), "\n",
        #"Median =", round(median(y), 2), "\n"
      )
    ))  
  }else{
    upper_limit = max(y) * 1.15
    return(data.frame(
      y = 0.95 * upper_limit,
      label = paste(
        "N =", length(y), "\n"#,
        #"Mean =", round(mean(y), 2), "\n",
        #"Median =", round(median(y), 2), "\n"
      )
    ))  
  }
  
}

# 1.1 get PSI data ----
pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}
psi_data <- metadata[,c("Donor",
                        "Sample",
                        "Ancestry",
                        "Sex",
                        "Age",
                        "BMI"),]
psi_data$Age_factor <- ifelse(psi_data$Age < 45, "[20-45)", "[45-70]")
psi_data$Age_factor <- factor(psi_data$Age_factor, levels = c("[20-45)", "[45-70]"), order = T)
psi_data$BMI_factor <- sapply(psi_data$BMI, function(bmi) pseudo_categorize_bmi(bmi))
psi_data$BMI_factor <- factor(psi_data$BMI_factor, levels = c("Normal", "Overweight", "Obese"), order = T)
colnames(psi_data)[colnames(psi_data)=="Age"] <- "Age_numeric"
colnames(psi_data)[colnames(psi_data)=="BMI"] <- "BMI_numeric"
colnames(psi_data)[colnames(psi_data)=="Age_factor"] <- "Age"
colnames(psi_data)[colnames(psi_data)=="BMI_factor"] <- "BMI"
psi_data$PSI <- sapply(psi_data$Sample, function(sample) psi[event_id, sample])
psi_data$event_label <- as.factor(event_label)

# 1.2 plot PSI data ----
event_label <- "CYP3A5 (A3) - LIVER"
psi_data$event_label <- event_label
p1 <- ggplot(data = psi_data,
       aes(x = eval(parse(text=trait)),
           y = PSI)
) +
  geom_violin(fill = traits_cols[trait],
              col = "black") +
  geom_boxplot(fill = "white",
               col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("PSI")+
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        #axis.ticks = element_line(size = 0.1),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        axis.text.x = element_text(angle = 0),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")  +
  facet_wrap(~event_label)
 
 
 #labs(title=event_label) 

# 1.3 plot % explained splicing variation ----
# read hier.part data  
hier_part <- readRDS("~/GTEx_v8/Raquel/Draft/Analysis/splicing/data/splicing.hier_parts.DSEs.rds")
100*hier_part[[tissue]][event_id,][,trait]

d <- data.frame("variable" = trait,
                "value"  = 100*hier_part[[tissue]][event_id,][,trait])
variable_col <- traits_cols[trait]
names(variable_col) <- trait

p2 <- ggplot(d, aes(x=variable, y=value)) + 
  geom_bar(aes(fill = as.factor(variable)), stat = "identity") +
  scale_fill_manual(values = variable_col) +
  coord_flip() +
  xlab("") + ylab("Splicing variation explained (%)") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks.x =  element_line(size = 0.5),
        axis.ticks.y =  element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = ggtext::element_markdown(size=10),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        axis.text = element_text(size = 12), 
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") 

  
# 1.4 PSI plot ----
PSI_plot <- ggarrange(p1, p2, nrow = 2, heights =  c(9,2))
pdf("~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/CYP3A5_A3_LIVER.pdf",
    width = 4, height = 5)
PSI_plot
dev.off()

# 2. Plot isoforms expression and abundance ratios ----
# for each event recover the biotype of the most abundant spliced-in isoform and spliced-out isoform --
get_isoform_info <- function(event.id){
  # Isoforms info --
  spliced.in.isoforms <- unlist(strsplit(events.info[event.id, "isoforms.spliced_in"], split = ","))
  spliced.out.isoforms <- unlist(strsplit(events.info[event.id, "isoforms.spliced_out"], split = ","))
  
  # Get isoform TPM data --
  tpm.data <- cbind.data.frame("Sample" = rep(metadata$Sample, length(c(spliced.in.isoforms, spliced.out.isoforms))),
                               "TPM" = unlist(lapply(c(spliced.in.isoforms, spliced.out.isoforms), function(isoform) 
                                 as.numeric(transcript.tpm[isoform,])))
  )
  tpm.data$Isoform <- unlist(lapply(c(spliced.in.isoforms, spliced.out.isoforms), function(isoform) rep(isoform, length(metadata$Sample))))
  tpm.data$Class <- ifelse(tpm.data$Isoform %in% spliced.in.isoforms,
                           "Spliced-in",
                           ifelse(tpm.data$Isoform %in% spliced.out.isoforms,
                                  "Spliced-out",       
                                  "Other"))
  
  # Order isoform by median TPM --
  spliced.in.isoform <- names(sort(sapply(spliced.in.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))[1]
  spliced.in.isoform.biotyoe <- transcript_annotation[transcript_annotation$transcript.id==spliced.in.isoform, "transcript_biotype"]
  spliced.out.isoform <- names(sort(sapply(spliced.out.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))[1]
  spliced.out.isoform.biotyoe <- transcript_annotation[transcript_annotation$transcript.id==spliced.out.isoform, "transcript_biotype"]
  return(c(spliced.in.isoform, spliced.in.isoform.biotyoe,
           spliced.out.isoform, spliced.out.isoform.biotyoe))  
}

# 2.1 Get isoform TPM and relative abundance data ----
# Isoforms info --
spliced.in.isoforms <- unlist(strsplit(event_annotation[event_id, "isoforms.spliced_in"], split = ","))
spliced.out.isoforms <- unlist(strsplit(event_annotation[event_id, "isoforms.spliced_out"], split = ","))
spliced.isoforms <- c(spliced.in.isoforms,
                      spliced.out.isoforms)
all.isoforms <- transcript_annotation[transcript_annotation$ensembl.id == event_annotation[event_id,"ensembl.id"], "transcript.id"]
other.isoforms <- all.isoforms[!all.isoforms %in% spliced.isoforms]

# Get isoform data --
tpm.data <- cbind.data.frame("Sample" = rep(metadata$Sample, length(all.isoforms)),
                             "TPM" = unlist(lapply(all.isoforms, function(isoform) 
                               as.numeric(transcript.tpm[isoform,]))),
                             "Ratio" = unlist(lapply(all.isoforms, function(isoform) 
                               as.numeric(transcript.sratio[isoform,])))
)
tpm.data$Isoform <- unlist(lapply(all.isoforms, function(isoform) rep(isoform, length(metadata$Sample))))
tpm.data$Class <- ifelse(tpm.data$Isoform %in% spliced.in.isoforms,
                         "Spliced-in",
                         ifelse(tpm.data$Isoform %in% spliced.out.isoforms,
                                "Spliced-out",       
                                "Other"))

# Order isoform by median TPM --
spliced.in.isoforms <- names(sort(sapply(spliced.in.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))
spliced.out.isoforms <- names(sort(sapply(spliced.out.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))
if(length(other.isoforms)>0){
  other.isoforms <- names(sort(sapply(other.isoforms, function(isoform) median(tpm.data[tpm.data$Isoform==isoform, "TPM"])), decreasing = T))  
}

# Prep dataframe --
tpm.data$Isoform <- factor(tpm.data$Isoform, levels = c(spliced.in.isoforms, spliced.out.isoforms, other.isoforms), order = T)
tpm.data$Class <- factor(tpm.data$Class, levels = c("Spliced-in", "Spliced-out","Other"), order = T)
tpm.data$Isoform2 <- sapply(as.character(tpm.data$Isoform), function(iso)
  paste0(transcript_annotation[transcript_annotation$transcript.id==iso, "transcript.name"], 
         "\n(", transcript_annotation[transcript_annotation$transcript.id==iso, "transcript_biotype"], ")")
)
tpm.data$Isoform2 <- factor(tpm.data$Isoform2,
                            levels = sapply(levels(tpm.data$Isoform), function(i)
                              unique(tpm.data[tpm.data$Isoform==i, "Isoform2"])
                            ),
                            order = T)


# 2.2 Subet most expressed spliced-in and spliced-out isoforms ----
tpm_data <- tpm.data[tpm.data$Isoform %in% c(spliced.in.isoforms[1], spliced.out.isoforms[1]),]
tpm_data$Isoform <- droplevels(tpm_data$Isoform)
tpm_data <- merge(tpm_data,
                  metadata[,c("Sample",
                              "Ancestry",
                              "Sex",
                              "Age",
                              "BMI"),],
                  by = "Sample")
tpm_data$Ancestry <- gsub("AFR", "AA", tpm_data$Ancestry)
tpm_data$Ancestry <- gsub("EUR", "EA", tpm_data$Ancestry)
tpm_data$Ancestry <- factor(tpm_data$Ancestry, levels = c("EA", "AA"), order = T)
tpm_data$Age_factor <- ifelse(tpm_data$Age < 45, "[20-45)", "[45-70]")
tpm_data$Age_factor <- factor(tpm_data$Age_factor, levels = c("[20-45)", "[45-70]"), order = T)
tpm_data$BMI_factor <- sapply(tpm_data$BMI, function(bmi) pseudo_categorize_bmi(bmi))
tpm_data$BMI_factor <- factor(tpm_data$BMI_factor, levels = c("Normal", "Overweight", "Obese"), order = T)
colnames(tpm_data)[colnames(tpm_data)=="Age"] <- "Age_numeric"
colnames(tpm_data)[colnames(tpm_data)=="BMI"] <- "BMI_numeric"
colnames(tpm_data)[colnames(tpm_data)=="Age_factor"] <- "Age"
colnames(tpm_data)[colnames(tpm_data)=="BMI_factor"] <- "BMI"

# 2.3 Plot isoforms TPM ----
p3 <- ggplot(data = tpm_data,
       aes(x = eval(parse(text=trait)),
           y = log2(1+TPM))
) +
  geom_violin(fill = traits_cols[trait],
              col = "black") +
  geom_boxplot(fill = "white",
               col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("log<sub>2</sub>(TPM)") +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.25),
        axis.line = element_line(colour = "black", size=0.25),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        axis.text.x = element_text(angle = 0),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")  +
  facet_wrap(~Isoform2)

# 2.4 Plot isoforms splicing ratios ----
p4 <- ggplot(data = tpm_data,
             aes(x = eval(parse(text=trait)),
                 y = Ratio)
) +
  geom_violin(fill = traits_cols[trait],
              col = "black") +
  geom_boxplot(fill = "white",
               col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("Splicing ratio") +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.25),
        axis.line = element_line(colour = "black", size=0.25),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        axis.text.x = element_text(angle = 0),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none")  +
  facet_wrap(~Isoform2)

# 3. Save plot ----
pdf_file <- paste0(outpath,
                   tissue, "_", trait, ".",
                   gsub(" ", "", gsub(" - ", "_", event_label)),
                   ".pdf")
pdf(pdf_file,
    width = 9, height = 4)
ggarrange(PSI_plot, p3, p4, nrow = 1)
dev.off()

