#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggtext)
library(ggpubr)

# metadata ----
load("00.metadata.RData")

# differential expression analysis ----
dea <- readRDS("differential_expression_analysis.rds")

# list of DEGs de per tissue and trait
get_deg <- function(tissue,trait){
  if(length(dea[[tissue]][[trait]])==1){
    return(NA)
  }else{
    de_genes <- rownames(dea[[tissue]][[trait]])[ dea[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(de_genes)
  }
}

# genes DE with each trait in each tissue
DEGs <- lapply(tissues, function(tissue)
  lapply(traits, function(trait) 
    get_deg(tissue,trait)
  )
)
names(DEGs) <- tissues
for(tissue in tissues){names(DEGs[[tissue]]) <- traits}

# combination of 2 traits 
traits_two <- c("Ancestry-Sex",
                "Ancestry-Age",
                "Ancestry-BMI",
                "Sex-Age",
                "Sex-BMI",
                "Age-BMI")

# list of genes with additive effect per tissue and combination of 2 traits
DEGs_2_traits <- lapply(traits_two, function(i)
  lapply(tissues, function(tissue) 
    intersect(DEGs[[tissue]][[unlist(strsplit(i, split = "-"))[[1]]]], DEGs[[tissue]][[unlist(strsplit(i, split = "-"))[[2]]]])
  )
)
names(DEGs_2_traits) <- traits_two
for(pw in traits_two){
  names(DEGs_2_traits[[pw]]) <- tissues
}


# Fisher's exact test to test if there are more DEGs with 2 traits than expected ----
overlap_enrichment_fun <- function(tissue, trait.1, trait.2){
  if(tissue %in% sex_tissues & (trait.1 == "Sex" | trait.2 == "Sex")){
    return(list("overlap" = NA,
                "odds.ratio" = NA,
                "p.value" = NA))
  }else{
    #                   de.trait.2  not.de.trait.2
    # de.trait.1
    # not.de.trait.1
    x11 <- length(intersect(DEGs[[tissue]][[trait.1]], DEGs[[tissue]][[trait.2]]))
    x12 <- length(DEGs[[tissue]][[trait.1]][! DEGs[[tissue]][[trait.1]] %in% DEGs[[tissue]][[trait.2]]])
    x21 <- length(DEGs[[tissue]][[trait.2]][! DEGs[[tissue]][[trait.2]] %in% DEGs[[tissue]][[trait.1]]])
    x22 <- length(exprs_genes[[tissue]][! exprs_genes[[tissue]] %in% unique(c(DEGs[[tissue]][[trait.1]], DEGs[[tissue]][[trait.2]]))])
    m <- matrix(c(x11,x12,x21,x22),2,2,byrow = T)
    rownames(m) <- c(paste0(trait.1,".deg"), paste0(trait.1, ".not_deg"))
    colnames(m) <- c(paste0(trait.2,".deg"), paste0(trait.2, ".not_deg"))
    #sum(m) == length(exprs.genes[[tissue]])
    f <- fisher.test(m)
    f$estimate
    return(list("overlap" = x11,
                "odds.ratio" = f$estimate,
                "p.value" = f$p.value,
                "counts.matrix" = m,
                "lower_CI" = f$conf.int[1],
                "upper_CI" = f$conf.int[2]))
  }
}

# Genes expressed per tissue --
exprs_genes <- lapply(tissues, function(tissue) rownames(dea[[tissue]][["Age"]]))
names(exprs_genes) <- tissues

# Enrichment analysis  --
fisher_results <- lapply(traits_two, function(i)
  lapply(tissues, function(tissue) 
    overlap_enrichment_fun(tissue, unlist(strsplit(i, split = "-"))[[1]], unlist(strsplit(i, split = "-"))[[2]])
  ))
names(fisher_results) <- traits_two
for(i in traits_two){names(fisher_results[[i]]) <- tissues}

# P-values
p_values <- lapply(traits_two, function(i)
  sapply(tissues, function(tissue)
    fisher_results[[i]][[tissue]][["p.value"]]
  )
)
names(p_values) <- traits_two
fisher_pvalues <- p_values

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p_values), method = "BH")
adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(tissues), ncol = length(traits_two),
                         byrow = F)

fdr <- -log10(adjPVal_matrix)
colnames(fdr) <- traits_two
rownames(fdr) <- tissues
fdr[fdr <= -log10(0.05)] <- NA
apply(fdr, 2, function(x) sum(!is.na(x)))

sum(apply(fdr, 2, function(x) sum(!is.na(x))))
fisher_tissues <- apply(fdr, 2, function(x) tissues[which(!is.na(x))])

# Table S3A ----
Table_S3A <- t(sapply(tissues, function(tissue) sapply(traits_two, function(traits) fisher_results[[traits]][[tissue]]$overlap)))

# Table S3B ----
my_tissues <- lapply(traits_two, function(i) tissues[!is.na(fdr[,i])])
names(my_tissues) <- traits_two
parse_table <- function(tissue, traits){
  return(c(tissue,
           traits,
           fisher_results[[traits]][[tissue]]$overlap,
           fisher_results[[traits]][[tissue]]$odds.ratio,
           fisher_results[[traits]][[tissue]]$p.value,
           fdr[tissue, traits]))
}

Table_S3B <- do.call(rbind.data.frame, lapply(names(my_tissues), function(traits) t(sapply(my_tissues[[traits]], function(tissue) parse_table(tissue, traits)))))
rownames(Table_S3B) <- NULL
colnames(Table_S3B) <- c("tissue", "traits", "overlap", "odds ratio (log2)", "P-value", "FDR (-log10)")
Table_S3B$`odds ratio (log2)` <- log2(as.numeric(Table_S3B$`odds ratio (log2)`))


# Is there a bias towards a particular direction of change ----
Xsq_fun <- function(tissue, trait1, trait2){
  #print(paste0(tissue, ": ", trait1, "-", trait2))
  if(tissue %in% sex_tissues & (trait1 == "Sex" | trait2 == "Sex")){
    return(list("P-value" = NA, 
                "O/E" = rep(NA, 4),
                "counts" = rep(NA, 4)))
  }else{
    # Do we observe a higher than expected overlap of DEGs in a particular direction of change?
    # a numeric vector representing the observed proportions
    # a vector of probabilities (of the same length of the observed proportions) representing the expected proportions
    trait1.up <- rownames(dea[[tissue]][[trait1]][dea[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                        dea[[tissue]][[trait1]]$logFC > 0,])
    trait1.down <- rownames(dea[[tissue]][[trait1]][dea[[tissue]][[trait1]]$adj.P.Val < 0.05 &
                                                          dea[[tissue]][[trait1]]$logFC < 0,])
    
    trait2.up <- rownames(dea[[tissue]][[trait2]][dea[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                        dea[[tissue]][[trait2]]$logFC > 0,])
    trait2.down <- rownames(dea[[tissue]][[trait2]][dea[[tissue]][[trait2]]$adj.P.Val < 0.05 &
                                                          dea[[tissue]][[trait2]]$logFC < 0,])
    
    # Observed counts
    counts <- c(sum(trait1.up %in% trait2.up), # upup
                sum(trait1.down %in% trait2.up), # downup
                sum(trait1.up %in% trait2.down), # updown
                sum(trait1.down %in% trait2.down) # downdown
    )
    # Expected proportions
    trait1.up.p <- length(trait1.up)/length(c(trait1.up, trait1.down))
    trait1.down.p <- length(trait1.down)/length(c(trait1.up, trait1.down))
    trait2.up.p <- length(trait2.up)/length(c(trait2.up, trait2.down))
    trait2.down.p <- length(trait2.down)/length(c(trait2.up, trait2.down))
    expected_prob <- c(trait1.up.p * trait2.up.p, trait1.down.p * trait2.up.p, trait1.up.p * trait2.down.p, trait1.down.p * trait2.down.p)
    expected_counts <- round(sum(counts) * expected_prob)
    
    # Return results
    if(sum(counts) < 20){
      print(paste0(tissue, ": Fewer than 20 genes DE with ", trait1, " and ", trait2))
      return(list("P-value" = NA,
                  "O/E" = rep(NA, 4),
                  "counts" = counts))
      break
    }else{
      if(min(expected_counts) < 5){
        print(paste0(tissue, ": Number of observations is not enough for Chi-Square Test\nUsing Monte Carlo simulations"))
        Xsq <- chisq.test(counts, 
                          p = expected_prob,
                          simulate.p.value = T)
      }else{
        Xsq <- chisq.test(counts,
                          p = expected_prob)
      }
      oe <- Xsq$observed/round(Xsq$expected)
      return(list("P-value" = Xsq$p.value,
                  "O/E" = oe,
                  "counts" = counts))
    }
  }
}

# Enrichment analysis ----
Xsq_results <- lapply(traits_two, function(i)
  lapply(tissues, function(tissue) 
    Xsq_fun(tissue,
            unlist(strsplit(i, split = "-"))[[1]],
            unlist(strsplit(i, split = "-"))[[2]])
  ))
names(Xsq_results) <- traits_two
for(i in traits_two){names(Xsq_results[[i]]) <- tissues}

# Enrichment statistics --
# P-values
p_values <- lapply(traits_two, function(i)
  sapply(tissues, function(tissue)
    Xsq_results[[i]][[tissue]][["P-value"]]
  )
)
names(p_values) <- traits_two
Xsq_pvalues <- p_values

# multiple testing correction across tissues and number of pairwise combinations of traits
adj_p_values <- p.adjust(unlist(p_values), method = "BH")
adjPVal_matrix <- matrix(adj_p_values, 
                         nrow = length(tissues), ncol = length(traits_two),
                         byrow = F)
fdr <- -log10(adjPVal_matrix)
colnames(fdr) <- traits_two
rownames(fdr) <- tissues
fdr[fdr <= -log10(0.05)] <- NA # if not tested NA (grey in plot); if FDR >= 0.05, 0 (white in plot)
sum(apply(fdr, 2, function(x) sum(!is.na(x))))
Xsq_tissues <- apply(fdr, 2, function(x) tissues[which(!is.na(x))])

m <- t(cbind.data.frame(
  sapply(traits_two, function(i) sum(sapply(tissues, function(tissue) length(DEGs_2_traits[[i]][[tissue]])) > 20) ), 
  sapply(traits_two, function(i) length(fisher_tissues[[i]])),
  sapply(traits_two, function(i) length(intersect(fisher_tissues[[i]], Xsq_tissues[[i]])))
  ))
m2 <- rbind.data.frame(
  m[1,] - m[2,],
  m[2,] - m[3,],
  m[3,]
)
colnames(m2) <- traits_two
rownames(m2) <- c("> 20 DEGs", "enriched","enriched & bias")
m2$type <- rownames(m2)
library(reshape2)
library(ggplot2)
d <- melt(m2)
d$type <- factor(d$type, levels = c("> 20 DEGs", "enriched","enriched & bias"), order = T)
p1 <- ggplot(d,
        aes(x = variable,
            y = value,
            col = type,
            fill = type)) +
  geom_bar(stat="identity") +
  ylab("Number of tissues") + xlab("") +
  scale_color_manual(values = brewer.pal(3, "Reds")) +
  scale_fill_manual(values = brewer.pal(3, "Reds")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
pdf(paste0(plot_path, "Figure_3A.pdf"),
    width = 6, height = 5)
p1
dev.off()

i <- "Sex-Age"
intersect(fisher_tissues[[i]], Xsq_tissues[[i]])
data1 <- as.data.frame(sapply(intersect(fisher_tissues[[i]], Xsq_tissues[[i]]), function(tissue) Xsq_results[[i]][[tissue]][["counts"]]))
rownames(data1) <- c("female - old", "male - old", "female - young", "male - young")
data1$type <- rownames(data1)
data2 <- as.data.frame(sapply(intersect(fisher_tissues[[i]], Xsq_tissues[[i]]), function(tissue) Xsq_results[[i]][[tissue]][["counts"]]/sum(Xsq_results[[i]][[tissue]][["counts"]])))
rownames(data2) <- c("female - old", "male - old", "female - young", "male - young")
data2$type <- rownames(data2)

data <- melt(data1)
data$y <- 100*melt(data2)[,3]
data$tissue <- sapply(data$variable, function(i) tissue_info[tissue_info$tissue_ID==i, "tissue_abbrv"])
data$tissue <- factor(data$tissue, levels = tissue_info$tissue_abbrv, order = T)
data$tissue <- droplevels(data$tissue)
data$type <- factor(data$type, levels = rev(c("male - young", "male - old", "female - young", "female - old")), order = T)
 ggplot(data,
       aes(x = tissue, 
           y = y, 
           fill = type,
           label = value) ) +
  geom_bar(stat = "identity") +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = rev(c("#00b159", "#00aedb", "#f37735", "#ffc425"))) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("DEGs (%)") + xlab("")

i <- "Sex-BMI"
intersect(fisher_tissues[[i]], Xsq_tissues[[i]])
data1 <- as.data.frame(sapply(intersect(fisher_tissues[[i]], Xsq_tissues[[i]]), function(tissue) Xsq_results[[i]][[tissue]][["counts"]]))
rownames(data1) <- c("female - high BMI", "male - high BMI", "female - low BMI", "male - low BMI")
data1$type <- rownames(data1)
data2 <- as.data.frame(sapply(intersect(fisher_tissues[[i]], Xsq_tissues[[i]]), function(tissue) Xsq_results[[i]][[tissue]][["counts"]]/sum(Xsq_results[[i]][[tissue]][["counts"]])))
rownames(data2) <- c("female - high BMI", "male - high BMI", "female - low BMI", "male - low BMI")
data2$type <- rownames(data2)

data <- melt(data1)
data$y <- 100*melt(data2)[,3]
data$tissue <- sapply(data$variable, function(i) tissue_info[tissue_info$tissue_ID==i, "tissue_abbrv"])
data$tissue <- factor(data$tissue, levels = tissue_info$tissue_abbrv, order = T)
data$tissue <- droplevels(data$tissue)
data$type <- factor(data$type, levels = rev(c("male - low BMI", "male - high BMI", "female - low BMI", "female - high BMI")), order = T)
p3 <- ggplot(data,
             aes(x = tissue, 
                 y = y, 
                 fill = type,
                 label = value) ) +
  geom_bar(stat = "identity") +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = brewer.pal(11, "RdBu")[c(2,3,10,9)]) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("DEGs (%)") + xlab("")

pdf(paste0(plot_path, "Figure_3b.pdf"),
    width = 9, height = 5)
ggarrange(p2, p3, widths = c(2.7,2))
dev.off()
