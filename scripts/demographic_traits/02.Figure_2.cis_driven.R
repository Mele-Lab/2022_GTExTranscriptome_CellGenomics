#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

# metadata ----
load("00.metadata.RData")

# ancestry eGenes DE
ancestry_eGenes_DE <- readRDS("ancestry_eGenes_DE.rds")

# Figure 2C ----
d <- cbind.data.frame("tissue" = tissues,
                      "cis-driven" = 100*sapply(tissues, function(tissue) nrow(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "cis-driven",])/nrow(ancestry_eGenes_DE[[tissue]])))
d$tissue <- factor(d$tissue, levels = tissues, order = T)
d$dummy <- rep("ancestry-DEGs (eGenes)", nrow(d))


cor.test(sapply(tissues, function(tissue) nrow(mdata[[tissue]])),
         sapply(tissues, function(tissue) nrow(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "cis-driven",])),
                method = "spearman")
)
cor.test(sapply(tissues, function(tissue) nrow(mdata[[tissue]])),
         sapply(tissues, function(tissue) nrow(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "cis-driven",])/nrow(ancestry_eGenes_DE[[tissue]])),
         method = "spearman"
)

# plot
p1 <- ggplot(d,
       aes(x = 1, y = `cis-driven`, col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Cis-driven genes (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") +
  facet_grid(~dummy)

# Figure 2D ---
Fst <- lapply(tissues, function(tissue) readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Fst/Tissues/",tissue,"/", tissue, ".eVariants.Fst.rds")))
names(Fst) <- tissues
Fst_data <- cbind.data.frame("tissue" = tissues, 
                 "cis-driven" = sapply(tissues, function(tissue) 
  median(Fst[[tissue]][rownames(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "cis-driven",]), "avg.Fst"], na.rm = T)),
  "not cis-driven" = sapply(tissues, function(tissue) 
    median(Fst[[tissue]][rownames(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "not cis-driven",]), "avg.Fst"], na.rm = T))
)
Fst_data <- melt(Fst_data)
Fst_data$variable <- factor(Fst_data$variable, levels = c("cis-driven", "not cis-driven"), order = T)

p2 <- ggplot(Fst_data,
       aes(x = variable, 
           y = value, 
           col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  #geom_line(aes(group=tissue)) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Fst") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=brewer.pal(11, "PRGn")[c(3,9)]),
        legend.position = "none") +
  facet_grid(~variable, scale = "free")

# Figure 2D ----
tissue_sharing <- readRDS("genes_DE.tissue_sharing.rds")
tissue_sharing_data <- cbind.data.frame("tissue" = tissues, 
                             "cis-driven" = sapply(tissues, function(tissue) 
                              median(tissue_sharing$Ancestry[tissue_sharing$Ancestry$ensembl_id %in% rownames(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "cis-driven",]), "n_DE"])),
                             "not cis-driven" = sapply(tissues, function(tissue) 
                               median(tissue_sharing$Ancestry[tissue_sharing$Ancestry$ensembl_id %in% rownames(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type == "not cis-driven",]), "n_DE"]))
)
tissue_sharing_data <- melt(tissue_sharing_data)
tissue_sharing_data$variable <- factor(tissue_sharing_data$variable, levels = c("cis-driven", "not cis-driven"), order = T)

p3 <- ggplot(tissue_sharing_data,
             aes(x = variable, 
                 y = value, 
                 col = tissue)) +
  geom_violin(col = "black") +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  #geom_line(aes(group=tissue)) +
  scale_color_manual(values = tissue_info$colcodes) +
  ylab("Number of tissues") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=brewer.pal(11, "PRGn")[c(3,9)]),
        legend.position = "none") +
  facet_grid(~variable, scale = "free")


# Figure 2E ----
data <- cbind.data.frame("tissue" = rep(tissues, 3),
                         "value" = c(sapply(tissues, function(tissue) summary(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type=="cis-driven", "ieQTLs_R2"])[3]),
                                     sapply(tissues, function(tissue) summary(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type=="not cis-driven", "ieQTLs_R2"])[3]),
                                     sapply(tissues, function(tissue) summary(ancestry_eGenes_DE[[tissue]][ancestry_eGenes_DE[[tissue]]$type=="not cis-driven", "Ancestry_R2"])[3])),
                         "feature" = c(rep("eQTLs", length(tissues)*2), rep("Ancestry", length(tissues))),
                         "type" = c(rep("cis-driven", length(tissues)), rep("cis-independent", length(tissues)*2)))
data$tissue <- factor(data$tissue, levels = tissues, order = T)
data$feature <- factor(data$feature, levels = c("eQTLs", "Ancestry"), order = T)
data$type <- factor(data$type, levels = c("cis-driven", "cis-independent"), order = T)
data$class <- paste0(data$feature, "_", data$type)
data$class <- factor(data$class, levels = c("eQTLs_cis-driven", "eQTLs_cis-independent", "Ancestry_cis-independent"), order = T)
unique(data$class)

p4 <- ggplot(data,
       aes(x= feature, y =  value, col = feature, fill = feature)) +
  geom_violin(col = "black", fill = NA) +
  geom_boxplot(col = "black",
               fill = "white",
               outlier.shape = NA,
               notch = T,
               width = 0.25) + 
  geom_jitter(aes(col = tissue), 
              size = 2) +
  scale_color_manual(values = tissue_info$colcodes) +
  facet_grid(~class, scale = "free") +
  ylab("Gene expression variation explained (%)") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5,
                                  size = 15),
        strip.background = element_rect(fill=traits_cols["Ancestry"]),
        legend.position = "none") 
p4
pdf(paste0(plot_path, "Figure_2CDEF.cis_driven.NatureGenetics.pdf"),
    width = 9, height = 3)
ggarrange(p1, p2, p3, p4, widths = c(1,2,2,3), ncol = 4)
dev.off()

# Figure 2F ----
# The expression data  corresponds to gene-level TPM quantifications from the GTEx v8 main paper (GTEx Consortium. 2020), which are available on the GTEx portal
# From each tissue, we used a subset of samples with available metadata for the covariates included in the linear regression (STAR methods)

get_box_stats <- function(y, upper_limit = max(y) * 1.15) {
  return(data.frame(
    y = upper_limit,
    label = paste(
      "n =", length(y)
    )
  ))
}

pseudo_categorize_bmi <- function(bmi){
  if(bmi < 25){
    return("Normal")
  }else if(bmi < 30){
    return("Overweight")
  }else{
    return("Obese")
  }
}

get_gene_data <- function(tissue, tpm, gene_name){
  ensembl_id <- gene_annotation[gene_annotation$gene.name==gene_name, "ensembl.id"]
  df <- cbind.data.frame(tpm[ensembl_id,])
  colnames(df) <- "TPM"
  df$Ancestry <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Ancestry"])
  df$Ancestry <- gsub("EUR", "EA", df$Ancestry)
  df$Ancestry <- gsub("AFR", "AA", df$Ancestry)
  df$Ancestry <- factor(df$Ancestry, levels = c("EA", "AA"), order = T)
  df$Sex <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Sex"])
  df$Sex <- gsub("1", "Male", df$Sex)
  df$Sex <- gsub("2", "Female", df$Sex)
  df$Sex <- factor(df$Sex, levels = c("Male", "Female"), order = T)
  df$Age_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"Age"])
  df$Age <- sapply(df$Age_int, function(i) ifelse(i<45, "[20-45)", "[45-70]"))
  df$Age <- factor(df$Age, 
                   levels = c("[20-45)", "[45-70]"),
                   order = T)
  df$BMI_int <- sapply(rownames(df), function(i) mdata[[tissue]][mdata[[tissue]]$Sample==i,"BMI"])
  df$BMI <- sapply(df$BMI_int, function(bmi) pseudo_categorize_bmi(bmi))
  df$BMI <- factor(df$BMI,
                   levels = c("Normal", "Overweight", "Obese"),
                   order = T)
  df$Tissue <- rep(tissue_info[tissue_info$tissue_ID==tissue, "tissue_abbrv"], nrow(df))
  df$ensembl_id <- rep(ensembl_id, nrow(df))
  return(df)
}



# get gene data --
tissue <- "SkinSunExposedLowerleg"
gene_name <- "PWP2"
gene <- gene_annotation[gene_annotation$gene.name == gene_name, "ensembl.id"]

# gene TPM expression data (data accesbible from https://gtexportal.org/home/datasets)
tpm <- readRDS(paste0("~/GTEx_v8/Raquel/Draft/Data/Tissues/", tissue, "/", tissue, ".SelectedSamples.TPM.PC_lincRNA.rds"))

# tissue independent eQTLs (data accesbible from https://gtexportal.org/home/datasets)
#eqtls_dir <- "/gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_eQTL_independent/"
tissue_id <- tissue_info[tissue_info$tissue_ID==tissue, "tissue_id"]
eqtls_dir <- "~/GTEx_v8_data/cisQTLs/GTEx_Analysis_v8_eQTL_independent/"
ieqtls <-  as.data.frame(data.table::fread(paste0(eqtls_dir, tissue_id,".v8.independent_eqtls.txt.gz")))

# i-eQTL of the gene with the lowest p-value
eVariant <- ieqtls[ieqtls$gene_id==gene,][which.min(ieqtls[ieqtls$gene_id==gene,"pval_true_df"]),][, "variant_id"]

# recover the samples' genotype for the i-eQTL associated with the cis-driven eGene DE between population
# read the genotypes for the eVariants which are also independent cis-eQtls 
# Job run in cluster
#vcftools --gzvcf ${vcf_file} --snps <(zcat /gpfs/projects/bsc83/Data/GTEx_v8/cisQTLs/GTEx_Analysis_v8_sQTL_independent/${tissue_id}.v8.independent_sqtls.txt.gz  | cut -f7 | awk '{if(NR>1)print}')  --recode --stdout | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz
#vcf-query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE=%GTR]\n' ${pathOut}/${tissue}.isQTLs.Donor_genotypes.vcf.gz | gzip -c > ${pathOut}/${tissue}.isQTLs.Donor_genotypes.txt.gz
#gt_path <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Fst/ieQTL_GT/"
gt_path <- "~/GTEx_v8/Raquel/Draft/Data/Fst/ieQTL_GT/"
variants_genotype <- as.data.frame(data.table::fread(paste0(gt_path, tissue, ".ieQTLs.Donor_genotypes.txt.gz"))) 

# set proper colnames 
colnames.variants <- c('CHROM','POS','variant_id','REF','ALT')
library(stringr)
colnames.individuals <- str_split_fixed(as.character(variants_genotype[1,-c(1:5)]),'=',2)[,1]
colnames(variants_genotype) <- c(colnames.variants,colnames.individuals)

# Retain selected tissue samples 
variants_genotype <- variants_genotype[,c(colnames.variants, mdata[[tissue]]$Donor)]
# Donors in same order as in metadata
if(!identical(mdata[[tissue]]$Donor, colnames(variants_genotype)[6:ncol(variants_genotype)])){
  print("Donors not in the same order in metadata and genotype tables")
  quit()
}

# Customize the genotypes files 
vg.tissue_mod <- cbind(variants_genotype[1:5], apply(variants_genotype[6:ncol(variants_genotype)],2, FUN=function(x){str_split_fixed(x,'=',2)[,2]}))
vg.tissue_wide <- cbind(vg.tissue_mod[1:5], apply(vg.tissue_mod[6:ncol(vg.tissue_mod)],2, FUN=function(x){ifelse(x=='0|0','0',
                                                                                                                 ifelse(x=='0|1','1',
                                                                                                                        ifelse(x=='1|1','2',NA)))}))
vg.tissue_long <- reshape2::melt(vg.tissue_wide,
                                 id.vars = colnames.variants,
                                 variable.name = 'Individual_ID',
                                 value.name = 'genotype')

# Subset tissue samples ----
vg.tissue_long <- vg.tissue_long[vg.tissue_long$Individual_ID %in% mdata[[tissue]]$Donor,]


gt_data <- vg.tissue_long[vg.tissue_long$variant_id==eVariant,]
gt_data$Ancestry <- sapply(gt_data$Individual_ID, function(i) mdata[[tissue]][mdata[[tissue]]$Donor==i, "Ancestry"])
genotype_table <- table(gt_data$Ancestry, gt_data$genotype)
ref <- unique(gt_data$REF)
alt <- unique(gt_data$ALT)
colnames(genotype_table) <- c(paste0(ref,ref), paste0(ref, alt), paste0(alt, alt)) 

# calculate allele frequency ----
gt_freq <- t(apply(as.matrix(table(gt_data$Ancestry, gt_data$genotype)), 1, function(x) x/sum(x)))
allele_freq <- cbind.data.frame("Ancestry" = c(rep("EA", 2),
                                      rep("AA", 2)),
                       "allele" = c(ref, alt, ref, alt),
                       "proportion" = c(c(gt_freq[1,1] + gt_freq[1,2]/2,
                                          gt_freq[1,3] + gt_freq[1,2]/2),
                                        c(gt_freq[2,1] + gt_freq[2,2]/2,
                                          gt_freq[2,3] + gt_freq[2,2]/2)
                       ))
allele_freq$Ancestry <- factor(allele_freq$Ancestry, levels = c("EA", "AA"), order = T)                
allele_freq$allele <- factor(allele_freq$allele, levels = rev(c(ref, alt)), order = T)

# plot allele frequency 
allele_cols <- brewer.pal(3,"Greys")[c(2,3)]
names(allele_cols) <- c(ref, alt)

allele_frequency_plot <- ggplot(allele_freq, mapping = aes(x = Ancestry, 
                                                  fill = allele, 
                                                  y = proportion, 
                                                  label = allele)) + 
  geom_bar( stat = "identity") +
  facet_grid(~Ancestry, scales = "free_x") +
  theme_bw() +
  ylab("Allele frequency") +
  scale_fill_manual(values = allele_cols) +
  xlab(paste0(unlist(strsplit(eVariant, split = "_"))[[1]],
              ": ",
              unlist(strsplit(eVariant, split = "_"))[[2]],
              "\n(Fst = 0.46)")) +
  theme(legend.position = 'none',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_text(size = 4, position = position_stack(vjust = 0.5))

data <- get_gene_data(tissue, tpm, gene_name)
data$individual <- sapply(rownames(data), function(i) mdata[[tissue]][mdata[[tissue]]$Sample == i, "Donor"])
data$GT <- sapply(data$individual, function(i) 
  ifelse(is.na(gt_data[gt_data$Individual_ID==i, "genotype"]), 
         NA,
         ifelse(gt_data[gt_data$Individual_ID==i, "genotype"]==0,
                paste0(gt_data[gt_data$Individual_ID==i, "REF"],"/",gt_data[gt_data$Individual_ID==i, "REF"]),
                ifelse(gt_data[gt_data$Individual_ID==i, "genotype"]==1,
                       paste0(gt_data[gt_data$Individual_ID==i, "REF"],"/",gt_data[gt_data$Individual_ID==i, "ALT"]),
                       paste0(gt_data[gt_data$Individual_ID==i, "ALT"],"/",gt_data[gt_data$Individual_ID==i, "ALT"])
                ))))         
data$GT <- factor(data$GT, 
                  levels = c(paste0(ref,"/",ref), paste0(ref,"/",alt), paste0(alt,"/",alt)),
                  order = TRUE)
data <- data[!is.na(data$GT),]

tpm_plot <- ggplot(data = data,
       aes(x = Ancestry,
           y = log2(TPM))) +
  geom_violin(aes(fill = Ancestry),
              col = "black") +
  geom_boxplot(col = "black",
               outlier.shape = NA,
               notch = T,
               width = 0.25) +
  geom_jitter(col = "black", 
              alpha = 0.1,
              size = 0.8) +
  xlab("") +
  ylab("TPM (log<sub>2</sub>)") +
  scale_fill_manual(values = c("#c9dc87",	"#87c9dc")) +
  stat_summary(fun.data = get_box_stats, geom = "text",
               hjust = 0.5, vjust = 0.9, size = 3) +
  #theme_bw() + 
  ggtitle(gene_name) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill="#B3B3B3"),
        legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  facet_grid(~GT)

pdf(paste0(plot_path, "Figure_2F.cis_driven_example.pdf"),
    width = 9, height = 4)
ggarrange(allele_frequency_plot, tpm_plot, widths = c(1,4))
dev.off()

pdf(paste0(plot_path, "Figure_2F.cis_driven_example.legend.pdf"),
    width = 4, height = 4)
plot.new()
legend("center",
       c("EA", "AA"),
       col = c("#c9dc87",	"#87c9dc"),
       pch = 15, bty = 'n', ncol = 1)
dev.off()
