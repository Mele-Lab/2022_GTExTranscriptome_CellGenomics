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

# Differential splicing analyses: results tables ----
# DSE 
fr_results <- lapply(tissues, function(tissue)
  readRDS(paste0("/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/02.DiffSplic/01.DSA/Tissues/",tissue,"/", tissue,".fractional_regression.covariates_and_traits.results.rds")))
names(fr_results) <- tissues

# Alternatively spliced events
AS_events <- lapply(tissues, function(tissue) 
  rownames(fr_results[[tissue]][["Age"]])
)
n_AS_tissues <- as.data.frame(table(unlist(AS_events)))
colnames(n_AS_tissues) <- c("event_id","n_tissues_AS")
n_AS_tissues$event_id <- as.character(n_AS_tissues$event_id)

##### Code #####
# 1. Create tables with tissue sharing of events DS ####
# list of evemt id of events ds per tissue and trait
get_dse <- function(tissue,trait){
  if(length(fr_results[[tissue]][[trait]])==1){
    return(NA)
  }else{
    ds_events <- rownames(fr_results[[tissue]][[trait]])[ fr_results[[tissue]][[trait]]$adj.P.Val < 0.05]
    return(ds_events)
  }
}
# One list per tissue with the event id of events DS with each trait
events_ds <- lapply(tissues, function(tissue)
  lapply(traits, function(trait) 
    get_dse(tissue,trait)
  )
)
names(events_ds) <- tissues
for(tissue in tissues){names(events_ds[[tissue]]) <- traits}

# list of all genes DE per trait ####
all_ds_events <- lapply(traits, function(trait) 
  unique(unlist( lapply(tissues, function(tissue) 
    events_ds[[tissue]][[trait]]))))
names(all_ds_events) <- traits
# remove NA
for(trait in traits){all_ds_events[[trait]] <- all_ds_events[[trait]][ !is.na( all_ds_events[[trait]] ) ]}

# Function  to calculate the tissue sharing of each gene DE with each trait
ds_event_sharing <- function(event,trait){
  n <- sum(sapply(tissues, function(tissue) ifelse(event %in% events_ds[[tissue]][[trait]],1,0)))
  avg_estimate <- mean(sapply(tissues, function(tissue) ifelse(event %in% events_ds[[tissue]][[trait]],
                                                               fr_results[[tissue]][[trait]][event,"Estimate"],
                                                               NA)), na.rm = T)
  up <- sum(sapply(tissues, function(tissue) ifelse(event %in% events_ds[[tissue]][[trait]],
                                                    ifelse(fr_results[[tissue]][[trait]][event,"Estimate"] > 0, 1, 0),
                                                    0)))
  down <- sum(sapply(tissues, function(tissue) ifelse(event %in% events_ds[[tissue]][[trait]],
                                                      ifelse(fr_results[[tissue]][[trait]][event,"Estimate"] < 0, 1, 0),
                                                      0)))
  sign <- ifelse(n == up, 
                 1,
                 ifelse( n == down,
                         -1,
                         0))
  ds_tissues <- sapply(tissues, function(tissue) ifelse(event %in% events_ds[[tissue]][[trait]],
                                                        tissue,
                                                        NA))
  ds_tissues <-  ds_tissues[ !is.na(ds_tissues)]      
  ds_tissues <- paste(ds_tissues, collapse = ";")
  
  n_tissues_AS <- n_AS_tissues[n_AS_tissues$event_id==event, "n_tissues_AS"]
  AS_tissues <- sapply(tissues, function(tissue) 
    ifelse(event %in% rownames(fr_results[[tissue]][[trait]]),
           tissue,
           NA)
  )
  AS_tissues <- AS_tissues[!is.na(AS_tissues)]
  AS_tissues <- paste(AS_tissues, collapse = ";")
  l <- list("event.id" = event,
            "n.ds" = n,
            "avg_Estimate" = avg_estimate,
            "spliced.in" = up,
            "spliced.out" = down,
            "sign" = sign,
            "which.ds" = ds_tissues,
            "n.as" = n_tissues_AS,
            "which.as" = AS_tissues)
  return(l)
}

# DSE sharing
# one table per trait with 
# Event id, N (number tissues), Average Estimate Up Down Sign and TissuesDE
dse_sharing <- lapply(traits, function(trait)
  do.call(rbind.data.frame,
          lapply(all_ds_events[[trait]], function(event)
            ds_event_sharing(event,trait) 
          )
  )
)
names(dse_sharing) <- traits
for(trait in traits){
  dse_sharing[[trait]] <- dse_sharing[[trait]][ order(dse_sharing[[trait]]$n.ds, decreasing = T), ]
}

# Add gene name ####
for(trait in traits){
  dse_sharing[[trait]]$ensembl.id <- sapply(dse_sharing[[trait]]$event.id, function(event_id)
    unlist(strsplit(event_id, split = ";"))[[1]]
  )
  dse_sharing[[trait]]$gene.name <- sapply(dse_sharing[[trait]]$ensembl.id, function(gene) gene_annotation[gene_annotation$ensembl.id==gene,"gene.name"])
  dse_sharing[[trait]]$Type <- sapply(dse_sharing[[trait]]$event.id, function(event_id)
    unlist(strsplit(unlist(strsplit(event_id, split = ";"))[[2]],split = ":"))[[1]]
  )
}

for(trait in traits){
  dse_sharing[[trait]] <- dse_sharing[[trait]][,c(1,12,10,11,2,3,4,5,6,8,7,9)] 
}

# Save DSE sharing data ----
outpath <- "/gpfs/projects/bsc83/Projects/GTEx_v8/Raquel/Draft/02.DiffSplic/Data/"
if(!dir.exists(outpath)){dir.create(outpath, recursive = T)}
saveRDS(events_ds, paste0(outpath,"Events_DS.rds"))
saveRDS(all_ds_events, paste0(outpath,"All_events_DS_per_trait.rds"))
saveRDS(dse_sharing, paste0(outpath, "Events_DS.Tissue_sharing.rds"))


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
