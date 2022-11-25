#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
table_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/tables/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}
if(!dir.exists(table_path)){dir.create(table_path, recursive = T)}

# metadata ----
load("00.metadata.RData")

# differential splicing tables ----
dsa <- readRDS("00.differential_splicing_analysis.rds")
for(tissue in tissues){
  if(tissue %in% sex_tissues){
    for(trait in traits[-2]){
      dsa[[tissue]][[trait]]$type <- sapply(rownames(dsa[[tissue]][[trait]]), function(e) unlist(strsplit(unlist(strsplit(e, split = ";"))[[2]], split = ":"))[[1]])
    } 
  }else{
    for(trait in traits){
      dsa[[tissue]][[trait]]$type <- sapply(rownames(dsa[[tissue]][[trait]]), function(e) unlist(strsplit(unlist(strsplit(e, split = ";"))[[2]], split = ":"))[[1]])
    }
  }
}


# get tables ----
my_fun <- function(trait, tissue){
  if(tissue %in% sex_tissues & trait == "Sex"){
    return(c(NA, NA, NA, NA))
  }else{
    return(c(sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "biotype"] == "NC-NC"),
      sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "biotype"] %in% c("NC-PC", "PC-NC")),
      sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "biotype"] == "PC-PC"),
      sum(apply(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05 &
                                 dsa[[tissue]][[trait]]$biotype == "PC-PC",c("spliced_in_domains", "spliced_out_domains")], 
            1,
            function(x) sum(!is.na(x))) > 0)
    ))
  }
}
my_tables <- lapply(traits, function(trait)
  t(sapply(tissues, function(tissue) 
    my_fun(trait, tissue)
    ))
  )
names(my_tables) <- traits
for(trait in traits){
  colnames(my_tables[[trait]]) <- c("NC-NC", "NC-PC/PC-NC", "PC-PC", "PC-PC:protein domain")
}

# save Table S4E ----
for(trait in traits){
  write.table(my_tables[[trait]],
              paste0(table_path, "Table_S4E.", trait, ".tab"),
              col.names = T, row.names = T,
              quote = F,
              sep = "\t")
  
}
