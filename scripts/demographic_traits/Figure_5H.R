#!/usr/bin/env Rscript

# Libraries ----
library(ggplot2)
library(ggpubr)
library(stringr)
library(RColorBrewer)

# DATA  ----
# Demograpchic traits ----
traits_cols <- c("Ancestry" = "#E69F00",
                 "Sex" = "#009E73",
                 "Age" = "#56B4E9",
                 "BMI" = "#CC79A7")
traits <- names(traits_cols)

# Tissues ----
tissue_info <- readRDS("~/GTEx_v8/Raquel/Draft/Data/Tissue_info.46_tissues.rds") # tissues ordered by sample size
tissues <- tissue_info$tissue_ID

# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
#gene_annotation  <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
gene_annotation  <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
# 26,196 genes

# Isoform annotation ----
isoform_annotation <- read.delim("~/GTEx_v8_data/GeneAnnotation/gencode.v26.GRCh38.transcripts.bed", header = F)
colnames(isoform_annotation) <- c("chr","start", "end","strand", "feature","ensembl.id", "transcript.id","gene.name", "gene_biotype", "transcript.name", "transcript_biotype")

# Event annotation ----
#events.info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")
event_anno <- readRDS("~/GTEx_v8/Raquel/Draft/SUPPA/gencode.v26.PC_lincRNA.biotype_matched_v38.splicing_events_coordinates.rds")

# Ribosomal proteins ----
# gene annotation
ribosomal_proteins <- gene_annotation[c(grep("^RPS", gene_annotation$gene.name),
                                        grep("^RPL", gene_annotation$gene.name)),] 
# isoform annotation
ribosomal_isoforms <- isoform_annotation[isoform_annotation$ensembl.id %in% ribosomal_proteins$ensembl.id,]
ribosomal_isoforms$transcript_id <- sapply(ribosomal_isoforms$transcript.id, function(i)  unlist(strsplit(i, split = "\\."))[[1]])

# splicing events
ribosomal_events <- event_anno[event_anno$ensembl.id %in% ribosomal_proteins$ensembl.id,]

# Differential splicing results ----
dsa_res <- readRDS("~/GTEx_v8/Raquel/Draft/Analysis/splicing/data/dsa_res.full.rds")

# DSE sharing ----
dse_sharing <- readRDS("~/GTEx_v8/Raquel/Draft/02.DiffSplic/Data/Events_DS.Tissue_sharing.rds")[["Ancestry"]]
dse_sharing <- dse_sharing[dse_sharing$ensembl.id %in% ribosomal_proteins$ensembl.id,]
#80/nrow(dse_sharing)
#sum(table(dse_sharing[dse_sharing$n.ds >1, "n.ds"]))
highly_shared_DSEs <- dse_sharing[dse_sharing$n.ds > 0, ]
#highly_shared_DSEs <- highly_shared_DSEs[!duplicated(highly_shared_DSEs$gene.name),]
events <- highly_shared_DSEs$event.id
length(events)
#events <- dsa_res[[tissue]][["Ancestry"]][dsa_res[[tissue]][["Ancestry"]]$adj.P.Val < 0.05 &
#                                            dsa_res[[tissue]][["Ancestry"]]$Ensembl_id %in% ribosomal_proteins$ensembl.id,]


# # Fst values ----
# Fst_values <- read.delim("~/GTEx_v8/Raquel/Draft/Data/Fst/GTEx_v8.MAF01.Diploid_sites.Fst.EUR_AFR_ancestry.tab", header = F)
# colnames(Fst_values) <- c("chr", "start", "Fst", "variant_id", "ref", " alt")

# Variants associated with DSEs in RPs ----
#infile <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/SUPPA/ancestry_DSEs_in_RPs.diploid_variants_MAF01_in_event_1Kb_window.tab", header = F)
infile <- read.delim("~/GTEx_v8/Raquel/Draft/SUPPA/ancestry_DSEs_in_RPs.diploid_variants_MAF01_in_event_1Kb_window.tab", header = F)
colnames(infile) <- c("event_label", "variants", "Fst")
rownames(infile) <- sapply(infile$event_label, function(event) rownames(event_anno[event_anno$event.label == event,]))
# remove duplicated variants if any
variants <- lapply(rownames(infile), function(event) unlist(strsplit(infile[event, "variants"], split = ",")))
names(variants) <- rownames(infile)
fst_values <- lapply(rownames(infile), function(event) unlist(strsplit(infile[event, "Fst"], split = ",")))
names(fst_values) <- rownames(infile)
variants_fst <- lapply(rownames(infile), function(event) unique(paste0(variants[[event]], ":", fst_values[[event]])))
names(variants_fst) <- rownames(infile)
infile$variants <- sapply(rownames(infile), function(event)
  paste(sapply(variants_fst[[event]], function(i) unlist(strsplit(i, split = ":"))[[1]]), collapse = ";")
)
infile$Fst <- sapply(rownames(infile), function(event)
  paste(sapply(variants_fst[[event]], function(i) unlist(strsplit(i, split = ":"))[[2]]), collapse = ";")
)

# GTEx v8: genotyped variants within a 1kb flank of ancestry-DSEs in RPs ----
#variants_genotype <- read.delim("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/SUPPA/GTEx_v8.Genotyped_Variants.MAF01.Variants_associated_with_DSEs_in_RPs.txt", header = F)
variants_genotype <- read.delim("~/GTEx_v8/Raquel/Draft/SUPPA/GTEx_v8.Genotyped_Variants.MAF01.Variants_associated_with_DSEs_in_RPs.txt", header = F)

# set proper column names --
colnames.variants <- c('chr','position','variant_id','ref','alt')
colnames.individuals <- str_split_fixed(as.character(variants_genotype[1,-c(1:5)]),'=',2)[,1]
colnames(variants_genotype) <- c(colnames.variants,colnames.individuals)

# customize the genotypes files --
variants_genotype_mod <- cbind(variants_genotype[1:5], apply(variants_genotype[6:ncol(variants_genotype)],2, FUN=function(x){str_split_fixed(x,'=',2)[,2]}))
variants_genotype_wide <- cbind(variants_genotype_mod[1:5], apply(variants_genotype_mod[6:ncol(variants_genotype_mod)],2, FUN=function(x){ifelse(x=='0|0','0',
                                                                                                                 ifelse(x=='0|1','1',
                                                                                                                        ifelse(x=='1|1','2',NA)))}))
variants_genotype_long <- reshape2::melt(variants_genotype_wide,
                                 id.vars = colnames.variants,
                                 variable.name = 'Individual_ID',
                                 value.name = 'genotype')
variants_genotype_long$Individual_ID <- as.character(variants_genotype_long$Individual_ID)

# ANALYSIS ----
get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  if(upper_limit < 1){
    return(data.frame(
      y = 0.95 * upper_limit,
      label = paste(
        "n =", length(y), "\n"#,
        #"Mean =", round(mean(y), 2), "\n",
        #"Median =", round(median(y), 2), "\n"
      )
    ))  
  }else{
    return(data.frame(
      y = min(y) * 0.9,
      label = paste(
        "n =", length(y), "\n"#,
        #"Mean =", round(mean(y), 2), "\n",
        #"Median =", round(median(y), 2), "\n"
      )
    ))  
  }
}

for(tissue in tissues[c(-3)]){
print(tissue)
#outpath <- paste0("~/GTEx_v8/Raquel/Draft/Analysis/splicing/ribosomal_proteins_have_large_population_differences/sQTLs/", tissue, "/")
#if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}

  # 1. tissue data ----
# Metadata ----
#metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/", tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")

# PSI values ----
#psi <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/02.DiffSplic/01.DSA/Tissues/", tissue, "/", tissue, ".Alternatively_spliced_events.PSI_values.rds"))
psi <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/02.DiffSplic/01.DSA/Tissues/", tissue, "/", tissue, ".Alternatively_spliced_events.PSI_values.rds"))

# sQTLs calling results ----
sqtls <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/02.DiffSplic/06.RPs_sQTLs/Tissues/", 
                        tissue, "/", tissue, ".ancestry_DSEs_in_RPs.candidate_sQTL.rds"))
sqtls <- sqtls[sqtls$adj_p_value < 0.05,]
sqtls <- sqtls[sqtls$class=="cis-driven",]
sqtls$event_label
head(sqtls)

#sqtls["ENSG00000197756.9;RI:chr2:216499949:216500031-216501341:216501467:+",]

# variants in events DS in RPs in tissue ----
variants <- infile[rownames(infile) %in% rownames(dsa_res[[tissue]][["Ancestry"]][dsa_res[[tissue]][["Ancestry"]]$adj.P.Val < 0.05,]),]
psi <- psi[rownames(dsa_res[[tissue]][["Ancestry"]][dsa_res[[tissue]][["Ancestry"]]$adj.P.Val < 0.05 &
                                                    dsa_res[[tissue]][["Ancestry"]]$Ensembl_id %in% ribosomal_proteins$ensembl.id,]),]
psi <- psi[rownames(psi) %in% rownames(sqtls),]
psi <- psi[rownames(psi) %in% events,]

if(nrow(psi)>0){
# 2. plot per event DS in RPs ----
for(event in rownames(psi)){
  print(event)
  #event <- rownames(psi)[1]
  #get_plot <- function(event){  
# variant with lowest adjusted p-value
variant_id <- sqtls[event, "variant"]
variant_Fst <- round(as.numeric(sqtls[event, "Fst"]), 2)

# genotype of 'sQTL' ----
gt_data <- variants_genotype_long[variants_genotype_long$variant_id == variant_id &
                                  variants_genotype_long$Individual_ID %in% metadata$Donor,]
ref <- unique(gt_data$ref)
alt <- unique(gt_data$alt)
head(gt_data)

# PSI data ----
psi_data <- metadata[,c("Donor",
                        "Sample",
                        "Ancestry"),]
psi_data$psi <- sapply(psi_data$Sample, function(sample) psi[event, sample])
psi_data$genotype <- unlist(lapply(psi_data$Donor, function(donor) gt_data[gt_data$Individual_ID==donor,"genotype"]))
psi_data$GT <- sapply(psi_data$genotype, function(gt) ifelse(is.na(gt),
                                                             NA,
                                                             ifelse(gt=="0", paste0(ref,"/",ref),
                                                                    ifelse(gt==1, paste0(ref,"/",alt), paste0(alt,"/",alt)))))
psi_data$GT <- factor(psi_data$GT, 
                      levels = c(paste0(ref,"/",ref), paste0(ref,"/",alt), paste0(alt,"/",alt)),
                      order = TRUE)
psi_data$Ancestry <- gsub("AFR", "AA", psi_data$Ancestry)
psi_data$Ancestry <- gsub("EUR", "EA", psi_data$Ancestry)
psi_data$Ancestry <- factor(psi_data$Ancestry, levels = c("EA", "AA"), order = T)
psi_data <- psi_data[!is.na(psi_data$GT),]

# compute allele frequency ----
# https://www.google.com/search?q=allele+frequency&sxsrf=APq-WBv1cf6eDgLmvyezRqA_OiBvdHo5fQ:1645636533039&tbm=isch&source=iu&ictx=1&vet=1&fir=IIPGH1APIleVxM%252CsKs7MBUhyapWEM%252C_%253B1_N8Dq2sZkthGM%252C8TfWmupJUNWcMM%252C_%253BhBa5UfDik_lfXM%252C4EIBzlaaOjSMIM%252C_%253B_f9Xb8Ocz5-nkM%252CQU5ZOglEFPkEyM%252C_%253BVUbGroKKrfm9mM%252CloOFseq9WBqglM%252C_%253BVZ7saU6AFrTgwM%252CGgxfzP9yWrzhDM%252C_%253BKXQBKiAazffEsM%252CjaHko7LDf5alVM%252C_&usg=AI4_-kQ1vvS-KU-c36NmwFvarkE5vlCdpA&sa=X&ved=2ahUKEwiMs-znqZb2AhUmxIUKHXegCi4Q_h16BAgqEAE#imgrc=1_N8Dq2sZkthGM&imgdii=VZ7saU6AFrTgwM

# calculate allele frequency ----
gt_freq <- t(apply(as.matrix(table(psi_data$Ancestry, psi_data$GT)), 1, function(x) x/sum(x)))
# c(gt_freq[1,1] + gt_freq[1,2]/2,
#   gt_freq[1,3] + gt_freq[2,2]/2)
d2 <- cbind.data.frame("Ancestry" = c(rep("EA", 2),
                                      rep("AA", 2)),
                       "allele" = c(ref, alt, ref, alt),
                       "proportion" = c(c(gt_freq[1,1] + gt_freq[1,2]/2,
                                          gt_freq[1,3] + gt_freq[1,2]/2),
                                        c(gt_freq[2,1] + gt_freq[2,2]/2,
                                          gt_freq[2,3] + gt_freq[2,2]/2)
                       ))
d2$Ancestry <- factor(d2$Ancestry, levels = c("EA", "AA"), order = T)                
d2$allele <- factor(d2$allele, levels = rev(c(ref, alt)), order = T)

# 2.1 plot allele frequency ----
allele_cols <- brewer.pal(3,"Greys")[c(2,3)]
names(allele_cols) <- c(ref, alt)

allele_frequency_plot <- ggplot(d2, mapping = aes(x = Ancestry, 
                               fill = allele, 
                               y = proportion,
                               label = allele)) + 
  geom_bar( stat = "identity") +
  facet_grid(~Ancestry, scales = "free_x") +
  theme_bw() +
  ylab("Allele frequency") +
  scale_fill_manual(values = allele_cols) +
  xlab(paste0(unlist(strsplit(variant_id, split = "_"))[[1]],
              " : ",
              unlist(strsplit(variant_id, split = "_"))[[2]],
              " (Fst = ", round(variant_Fst, 2), ")")) +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13)) +
  geom_text(size = 4, position = position_stack(vjust = 0.5))

# 2.2 plot PSi stratified vy population and variant genotype ----
event_label <- paste(unlist(strsplit(sqtls[event, "event_label"], split = ":"))[c(1:3)], collapse = ":")
outpath <- paste0("~/GTEx_v8/Raquel/Draft/Analysis/splicing/ribosomal_proteins_have_large_population_differences/sQTLs_by_event/", gsub(":", "_", event_label), "/")
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
psi_plot <- ggplot(data = psi_data,
             aes(x = Ancestry,
                 y = psi)) +
  geom_violin(aes(fill = Ancestry),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") + ylab("PSI")+
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  labs(title = paste0(paste(unlist(strsplit(event_label,split = ":"))[1:3], collapse = ": "),
                      " (",  tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"],")")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        #axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") +
  facet_grid(~GT) +
  scale_fill_manual(values = c("gray", "#E69F00")) 
  pdf("~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/Figure_5H.pdf",
      width = 10, height = 4)
  ggarrange(allele_frequency_plot, psi_plot, widths = c(1, 4), nrow = 1) 
  dev.off()
  #%>%
    #ggexport(filename = gsub(":", "_", paste0(outpath, tissue, ".", event_label, ".png")),
    #          width = 800, height = 400)
  #png(gsub(":", "_", paste0(outpath, tissue, ".", event_label, ".png")), width = 800, height = 400)  
  #p
  #dev.off()
  #return(psi_plot)
}
}
}

my_plots <- sapply(rownames(psi)[1:9], function(event) get_plot(event), simplify = T)
ggsave("~/Desktop/test2.png", my_plots[[1]])
ggarrange(plotlist = my_plots)
plot_sQTL(rownames(psi)[3])

# RPLP2 ----
event <- my_events[5]
tissue <- my_tissues[[event]][1]
event_label <- event_anno[event, "event.label"]
variant.df <- closest_variant[closest_variant$event.label == event_label,]

# Metadata ----
#metadata <- readRDS( paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/", tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata <- readRDS( paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.SampleMetadata.rds"))
metadata$Ancestry <- relevel(metadata$Ancestry, ref="EUR")

# Closest variant genotype ----
ref <- closest_variant[closest_variant$event.label == event_label, "ref"]
alt <- closest_variant[closest_variant$event.label == event_label, "alt"]
Fst.value <- closest_variant[closest_variant$event.label == event_label, "Fst"]
variant_id <- closest_variant[closest_variant$event.label == event_label, "variant.id"]
gt.data <- gt.long[gt.long$variant.id == variant_id &
                     gt.long$Donor %in% metadata$Donor,]
gt.data <- gt.data[!duplicated(gt.data$Donor),]


# PSI data ----
psi.data <- metadata[,c("Donor",
                        "Sample",
                        "Ancestry"),]
#psi.data$psi <- sapply(psi.data$Sample, function(sample) psi[event, sample])
psi.data$genotype <- unlist(lapply(psi.data$Donor, function(donor) gt.data[gt.data$Donor==donor,"Genotype"]))
psi.data$GT <- sapply(psi.data$genotype, function(gt) ifelse(is.na(gt),
                                                             NA,
                                                             ifelse(gt=="0", paste0(ref,"/",ref),
                                                                    ifelse(gt==1, paste0(ref,"/",alt), paste0(alt,"/",alt)))))
psi.data$GT <- factor(psi.data$GT, 
                      levels = c(paste0(ref,"/",ref), paste0(ref,"/",alt), paste0(alt,"/",alt)),
                      order = TRUE)
psi.data$Ancestry <- gsub("AFR", "AA", psi.data$Ancestry)
psi.data$Ancestry <- gsub("EUR", "EA", psi.data$Ancestry)
psi.data$Ancestry <- factor(psi.data$Ancestry, levels = c("EA", "AA"), order = T)
psi.data <- psi.data[!is.na(psi.data$GT),]
table(psi.data$Ancestry, psi.data$GT)
barplot(apply(as.matrix(table(psi.data$Ancestry, psi.data$GT)), 2, function(x) x/sum(x)))

# Event PSI values ----
#psi <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
psi <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.PSI_values.Splicing_events.PC_lincRNA.rds"))
psi.data$psi <- sapply(psi.data$Sample, function(sample) psi[event, sample])

# Event TPM values ----
#event.tpm <- readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Paper/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))
#event.tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/",tissue,"/",tissue,".SelectedSamples.TPM.Splicing_events.PC_lincRNA.rds"))

# Transcript TPM ----
#transcript.tpm <- read.delim(paste0("~/GTEx_v8/Raquel/Draft/SUPPA/TranscriptExpressionFiles/", tissue, ".transcript_TPM.txt"))
#colnames(transcript.tpm) <- gsub("\\.", "-", colnames(transcript.tpm))
# Subset to tissue samples --
#transcript.tpm <- transcript.tpm[, metadata$Sample]
#identical(colnames(transcript.tpm), metadata$Sample)

# Transcript splicing ratio ----
#transcript.sratio <- read.delim(paste0("~/GTEx_v8/Raquel/Draft/SUPPA/PSI_isoforms/", tissue, "/", tissue, "_isoform.psi"))
#colnames(transcript.sratio) <- gsub("\\.", "-", colnames(transcript.sratio))
#rownames(transcript.sratio) <- sapply(rownames(transcript.sratio), function(i) unlist(strsplit(i, split = ";"))[[2]])
# Subset to tissue samples --
#transcript.sratio <- transcript.sratio[, metadata$Sample]
#identical(colnames(transcript.sratio), metadata$Sample)


# 3.2 Even PSI stratified per population and genotype ----
ggplot(data = psi.data,
       aes(x = Ancestry,
           y = psi)
) +
  geom_violin(aes(fill = Ancestry),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") + ylab("PSI")+
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 2) +
  labs(title = paste0(paste(unlist(strsplit(event_label,split = ":"))[1:3], collapse = ": "),
                      " (",  tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"],")")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5,
                                  size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none") +
  facet_grid(~GT) +
  scale_fill_manual(values = c("gray", "#E69F00")) 


head(closest_variant)


d <- cbind.data.frame("Ancestry" = c(rep("EA", length(unique(psi_data$GT))),
                                     rep("AA", length(unique(psi_data$GT)))),
                      "GT" = c(names(table(psi_data$Ancestry, psi_data$GT)[1,]),
                               names(table(psi_data$Ancestry, psi_data$GT)[2,])),
                      "count" = c(table(psi_data$Ancestry, psi_data$GT)[1,],
                                  table(psi_data$Ancestry, psi_data$GT)[2,]),
                      "proportion" = c(apply(as.matrix(table(psi_data$Ancestry, psi_data$GT)), 2, function(x) x/sum(x))[1,],
                                       apply(as.matrix(table(psi_data$Ancestry, psi_data$GT)), 2, function(x) x/sum(x))[2,])
)
d$Ancestry <- factor(d$Ancestry, levels = rev(c("EA", "AA")), order = T)

# ggplot(d, mapping = aes(x = GT, 
#                               fill = Ancestry, 
#                               y = proportion)) + 
#   geom_bar( stat = "identity") +
#   facet_grid(~GT, scales = "free_x") +
#   theme_bw() +
#   ylab("Proportion of individuals") +
#   scale_fill_manual(values = rev(c("gray", "#E69F00"))) +
#   xlab(paste0(unlist(strsplit(variant_id, split = "_"))[[1]],
#               " : ",
#               unlist(strsplit(variant, split = "_"))[[2]],
#               " (Fst = ", variant_Fst, ")")) +
#   theme(legend.position = 'none') 
