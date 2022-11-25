#!/usr/bin/env Rscript

# **************************************************************************** #
Sys.setenv(TZ="Europe/Madrid")
# ---------------------- #
start_time <- Sys.time()
# ---------------------- #
# **************************************************************************** #

# ---- Data ----
# Traits ----
traits <- c("Age","Ancestry","Sex","BMI")
# Tissues ----
tissue_info <- readRDS("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/Data/Tissue_info.46_tissues.rds") # tissues ordered by sample size
#tissue_info <- readRDS("~/GTEx_v8/Raquel/Draft/Data/Tissue_info.47_tissues.rds")
tissues <- tissue_info$tissue_ID

# Gene annotation ----
# PCG and lincRNA genes with mathched biotype annotation in gencode v38
# PAR genes excluded
gene_annotation  <- read.delim("/gpfs/projects/bsc83/Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
#gene_annotation  <- read.delim("~/bsc83_Data/GTEx_v8/GeneAnnotation/gencode.v26.GRCh38.genes.biotype_matched_v38.bed")
# 26,196 genes
Y_genes <- gene_annotation[gene_annotation$chr=="chrY",]$ensembl.id

# Differential expression analyses: results tables ----
# DEG ----
vl_results <- lapply(tissues, function(tissue)
  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/01.DEA/Tissues/",tissue,"/", tissue,".voom_limma.covariates_and_traits.results.rds")))
#readRDS(paste0("~/GTEx_v8/Raquel/03_DEA/Tissues/",tissue,"/", tissue,".voom_limma.covariates_and_traits.results.rds")))
names(vl_results) <- tissues

# Expressed genes ----
exprs_genes <- lapply(tissues, function(tissue) 
  rownames(vl_results[[tissue]][["Age"]])
)
n_exprs_tissues <- as.data.frame(table(unlist(exprs_genes)))
colnames(n_exprs_tissues) <- c("ensembl_id","n_tissues_expressed")
n_exprs_tissues$ensembl_id <- as.character(n_exprs_tissues$ensembl_id)

# ---- Code ----
# 1. Create tables with tissue sharing of genes DE ----
# list of ensembl.id of genes de per tissue and trait
get_deg <- function(tissue,trait){
  if(length(vl_results[[tissue]][[trait]])==1){
    return(NA)
  }else{
    de_genes <- rownames(vl_results[[tissue]][[trait]])[ vl_results[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(de_genes)
  }
}
# One list per tissue with the ensembl id of genes DE with each trait
genes_de <- lapply(tissues, function(tissue)
  lapply(traits, function(trait) 
    get_deg(tissue,trait)
  )
)
names(genes_de) <- tissues
for(tissue in tissues){names(genes_de[[tissue]]) <- traits}

# list of all genes DE per trait 
all_de_genes <- lapply(traits, function(trait) 
  unique(unlist( lapply(tissues, function(tissue) 
    genes_de[[tissue]][[trait]]))))
names(all_de_genes) <- traits
# remove NA
for(trait in traits){all_de_genes[[trait]] <- all_de_genes[[trait]][ !is.na( all_de_genes[[trait]] ) ]}

# Function  to calculate the tissue sharing of each gene DE with each trait
de_gene_sharing <- function(gene,trait){
  n <- sum(sapply(tissues, function(tissue) ifelse(gene %in% genes_de[[tissue]][[trait]],1,0)))
  avg_fc <- mean(sapply(tissues, function(tissue) ifelse(gene %in% genes_de[[tissue]][[trait]],
                                                         vl_results[[tissue]][[trait]][gene,"logFC"],
                                                         NA)), na.rm = T)
  up <- sum(sapply(tissues, function(tissue) ifelse(gene %in% genes_de[[tissue]][[trait]],
                                                    ifelse(vl_results[[tissue]][[trait]][gene,"logFC"] > 0, 1, 0),
                                                    0)))
  down <- sum(sapply(tissues, function(tissue) ifelse(gene %in% genes_de[[tissue]][[trait]],
                                                      ifelse(vl_results[[tissue]][[trait]][gene,"logFC"] < 0, 1, 0),
                                                      0)))
  sign <- ifelse(n == up, 
                 1,
                 ifelse( n == down,
                         -1,
                         0))
  de_tissues <- sapply(tissues, function(tissue) ifelse(gene %in% genes_de[[tissue]][[trait]],
                                                        tissue,
                                                        NA))
  de_tissues <-  de_tissues[ !is.na(de_tissues)]      
  de_tissues <- paste(de_tissues, collapse = ";")
  n_tissues_exprs <- n_exprs_tissues[n_exprs_tissues$ensembl_id==gene, "n_tissues_expressed"]
  exprs_tissues <- sapply(tissues, function(tissue) 
    ifelse(gene %in% rownames(vl_results[[tissue]][[trait]]),
           tissue,
           NA)
  )
  exprs_tissues <- exprs_tissues[!is.na(exprs_tissues)]
  exprs_tissues <- paste(exprs_tissues, collapse = ";")
  l <- list("ensembl.id" = gene,
            "n.de" = n,
            "avg_FC" = avg_fc,
            "up" = up,
            "down" = down,
            "sign" = sign,
            "which.de" = de_tissues,
            "n.exprs" = n_tissues_exprs,
            "which.exprs" = exprs_tissues)
  return(l)
}

# DEG sharing  ----
# one table per trait with 
deg_sharing <- lapply(traits, function(trait)
  do.call(rbind.data.frame,
          lapply(all_de_genes[[trait]], function(gene)
            de_gene_sharing(gene,trait) 
          )
  )
)
names(deg_sharing) <- traits
for(trait in traits){
  deg_sharing[[trait]] <- deg_sharing[[trait]][ order(deg_sharing[[trait]]$n.de, decreasing = T), ]
}

# Add gene name ----
for(trait in traits){
  deg_sharing[[trait]]$gene.name <- sapply(deg_sharing[[trait]]$ensembl.id, function(gene) gene_annotation[gene_annotation$ensembl.id==gene,"gene.name"])
}

# Reorganize
for(trait in traits){
  deg_sharing[[trait]] <- deg_sharing[[trait]][,c(1,10,2,3,4,5,6,8,7,9)]
}

# Save DE data ----
outpath <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/01.DiffExprs/Data/"
#outpath <- "~/GTEx_v8/Raquel/03_DEA/Data/Tissue_Sharing/"
dir.create(outpath, recursive = T)
saveRDS(genes_de, paste0(outpath,"Genes_DE.rds"))
saveRDS(all_de_genes, paste0(outpath,"All_genes_DE_per_trait.rds"))
saveRDS(deg_sharing, paste0(outpath, "Genes_DE.Tissue_sharing.rds"))

# **************************************************************************** #
# ---------------------- #
end_time <- Sys.time()
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Elapsed time ---- #")
end_time - start_time
print("# ---------------------- #\n")
# ---------------------- #

# ---------------------- #
print("")
print("# ---- Session info ---- #")
sessionInfo()
print("# ---------------------- #\n")
print("")
# ---------------------- #
# **************************************************************************** #
q()
# ------------------------------------ #

# Plots ####
trait <- "Age"
m <- sapply(1:length(tissues), function(i)
  sapply(1:length(tissues), function(j)
    nrow(deg_sharing[[trait]][deg_sharing[[trait]]$N==i &
                                deg_sharing[[trait]]$N_tissues_exprs==j,])
  )
)
colSums(m)[2]
nrow(deg_sharing[[trait]][deg_sharing[[trait]]$N==2,])
# Column -> tissues DE
# Row -> tissue exprs
nrow(deg_sharing[[trait]][deg_sharing[[trait]]$N==1 &
                            deg_sharing[[trait]]$N_tissues_exprs==2,])
head(m[,1:4])
col_ha <- HeatmapAnnotation("DEG (log10)" = anno_barplot(log10(colSums(m)+1),
                                                         gp = gpar(fill = traits_cols[trait],
                                                                   col = traits_cols[trait]),
                                                         border= F),
                            "Tissues" = anno_barplot(1:46,
                                                     gp = gpar(fill = "light grey",
                                                               col = "light grey"),
                                                     border = F),
                            annotation_name_gp = gpar(fontsize = 6))

row_ha <- HeatmapAnnotation("EG (log10)" = anno_barplot(log10(rowSums(m)+1),
                                                        gp = gpar(fill = "dark grey",
                                                                  col = "dark grey"),
                                                        border = F ),
                            "Tissues" = anno_barplot(1:46,
                                                     gp = gpar(fill = "light grey",
                                                               col = "light grey"),
                                                     border = F),
                            which = "row",
                            annotation_name_gp = gpar(fontsize = 6),
                            annotation_name_rot = 90)
m_anno <- m
m_anno[m_anno==0] <- ""
Heatmap(as.matrix(log10(m+1)),
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = col_ha,
        left_annotation = row_ha,
        col = brewer.pal(9,"Blues"),
        name = "DEG (log10)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", m_anno[i, j]), x, y, gp = gpar(fontsize = 6))
        })

plot(cumsum(m[,1]))
plot(cumsum(unlist(lapply(1:length(m[,1]), function(i) 
  rep(i,m[,1][i])))/colSums(m)[1]))
