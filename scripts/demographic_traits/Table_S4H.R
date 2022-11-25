#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
table_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/tables/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}
if(!dir.exists(table_path)){dir.create(table_path, recursive = T)}

# libraries ----
library(rcompanion)
library(gtools)

# metadata ----
load("00.metadata.RData")

# differential splicing tables ----
dsa <- readRDS("differential_splicing_analysis.rds")
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

# comparison of the amount of alternative splicing variations explained by the different types of splicing events ----
# pairs of events
event_combinations <- as.data.frame(combinations(splicing_events, n = length(splicing_events), r = 2))
colnames(event_combinations) <- c("type2", "type1")  
event_combinations <- event_combinations[,c(2,1)]
event_combinations$type1 <- factor(event_combinations$type1, levels = splicing_events, order = T)
event_combinations$type2 <- factor(event_combinations$type2, levels = splicing_events, order = T)
event_combinations <- event_combinations[order(event_combinations$type1, event_combinations$type2),]
m <- matrix(NA, 7,7)
rownames(m) <- splicing_events
colnames(m) <- splicing_events
for(i in 1:nrow(m)){
  for(j in 1:ncol(m)){
    m[i,j] <- paste0(splicing_events[i], "-",  splicing_events[j])
  }
}
m[upper.tri(m)]

effect_size_fun <- function(event_1, event_2, r2_values){
  if(event_1 == event_2){
    return(NA)
  }else{
    x <-  unlist(r2_values[c(event_1, event_2)])
    g <- unlist(lapply(c(event_1, event_2), function(event) rep(event, length(r2_values[[event]]))))
    g <- factor(g, levels = c(event_1, event_2), order = T)
    p <- tryCatch(
      {
        e <- wilcoxonRG(x = x,
                        g = g)
      },
      warning=function(cond) {
        print(paste0(event_1, "-", event_2))
        e <- wilcoxonRG(x = x,
                        g = g,)
      }
    )  
    return(e)  
  }
}

run_analysis <- function(tissue, trait){
  print(tissue)
  # r2 values associatd with each type of event
  r2_values <- lapply(splicing_events, function(event) 
    dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$type == event, "R2"][!is.na(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$type == event, "R2"])]
  )
  names(r2_values) <- splicing_events
  # kruskal-wallis test
  kruskal_p_value <- kruskal.test(r2_values)$p.value
  if(kruskal_p_value < 0.05){
    # pairwise wilcoxon test
    g <- unlist(lapply(splicing_events, function(event) rep(event, length(r2_values[[event]]))))
    g <- factor(g, levels = splicing_events, order = T)
    mu_p_values <- t(pairwise.wilcox.test(x = unlist(r2_values),
                                          g = g,
                                          p.adjust.method = "BH"
    )$p.value)
    mu_p_values <- rbind.data.frame(cbind.data.frame(rep(NA, 6), mu_p_values),
                                    rep(NA, 7))
    colnames(mu_p_values)[1] <- "SE"
    rownames(mu_p_values)[7] <- "AL"
    mu_p_values[mu_p_values >= 0.05] <- NA
    p_values <- mu_p_values[upper.tri(mu_p_values)]
    names(p_values) <- m[upper.tri(m)]
  }else{
    p_values <- rep(NA, length(m[upper.tri(m)]))
    names(p_values) <- m[upper.tri(m)]
  }
  
  # pairwise effect size
  effect_sizes <-  sapply(splicing_events, function(event_1)
    sapply(splicing_events, function(event_2)
      effect_size_fun(event_1, event_2, r2_values)
    ))
  effect_sizes <- round(effect_sizes, 2)
  effect_sizes[lower.tri(effect_sizes)] <- NA
  effect_sizes <- effect_sizes*(-1) # i changed the sign cause the matrix was fill by column, but i want to read it by row
  es <- effect_sizes[upper.tri(effect_sizes)]
  names(es) <-  m[upper.tri(m)]
  return(list("kruskal_p_value" = kruskal_p_value,
              "fdr" = p_values,
              "effect_sizes" = es))
}

# Ancestry ----
trait <- "Ancestry"
# number of DSEs per type of splicing events
number_of_DSEs <- t(sapply(tissues, function(tissue) 
  sapply(splicing_events, function(event) 
    sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "type"]==event)
  )))
# tissues with at least 3 DSEs of each type of splicing event
my_tissues <- tissues[apply(number_of_DSEs, 1, function(x) sum(x>2))==7]
# run analaysis
results <- lapply(my_tissues, function(tissue) run_analysis(tissue, trait))
names(results) <- my_tissues
# parse and create results table
x <- grep("SE", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- grep("MX", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- y[!y%in%x]
z <- grep("A5", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
z <- z[!z%in%c(x,y)]
w <- grep("A3", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
w <- w[!w%in%c(x,y,z)]
a <- grep("RI", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
a <- a[!a%in%c(x,y,z,w)]
b <- grep("AF", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
b <- b[!b%in%c(x,y,z,w,a)]
c <- grep("AL", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
c <- c[!c%in%c(x,y,z,w,a,b)]
column_order <- c(x,y,z,w,a,b,c)

df1 <- cbind.data.frame(sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                        t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))[, column_order]
) 
df2 <- cbind.data.frame( sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                       t(sapply(my_tissues, function(tissue) results[[tissue]]$effect_sizes))[, column_order]
)
df2[is.na(df1)] <- NA
df <- cbind.data.frame(df1, df2[,-1])
df <- df[df[,1] < 0.05, apply(df, 2, function(x) sum(!is.na(x))) > 0]

write.table(df,
            paste0(table_path, "Table_S4H.", trait, ".tab"),
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")

# Sex ----
trait <- "Sex"
# number of DSEs per type of splicing events
number_of_DSEs <- t(sapply(tissues, function(tissue) 
  sapply(splicing_events, function(event) 
    sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "type"]==event)
  )))
# tissues with at least 3 DSEs of each type of splicing event
my_tissues <- tissues[apply(number_of_DSEs, 1, function(x) sum(x>2))==7]
# run analaysis
results <- lapply(my_tissues, function(tissue) run_analysis(tissue, trait))
names(results) <- my_tissues
# parse and create results table
x <- grep("SE", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- grep("MX", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- y[!y%in%x]
z <- grep("A5", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
z <- z[!z%in%c(x,y)]
w <- grep("A3", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
w <- w[!w%in%c(x,y,z)]
a <- grep("RI", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
a <- a[!a%in%c(x,y,z,w)]
b <- grep("AF", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
b <- b[!b%in%c(x,y,z,w,a)]
c <- grep("AL", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
c <- c[!c%in%c(x,y,z,w,a,b)]
column_order <- c(x,y,z,w,a,b,c)

df1 <- cbind.data.frame(sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                        t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))[, column_order]
) 
df2 <- cbind.data.frame( sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                         t(sapply(my_tissues, function(tissue) results[[tissue]]$effect_sizes))[, column_order]
)
df2[is.na(df1)] <- NA
df <- cbind.data.frame(df1, df2[,-1])
df <- df[df[,1] < 0.05, apply(df, 2, function(x) sum(!is.na(x))) > 0]
write.table(df,
            paste0(table_path, "Table_S4H.", trait, ".tab"),
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")

# Age ----
trait <- "Age"
# number of DSEs per type of splicing events
number_of_DSEs <- t(sapply(tissues, function(tissue) 
  sapply(splicing_events, function(event) 
    sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "type"]==event)
  )))
# tissues with at least 3 DSEs of each type of splicing event
my_tissues <- tissues[apply(number_of_DSEs, 1, function(x) sum(x>2))==7]
# run analaysis
results <- lapply(my_tissues, function(tissue) run_analysis(tissue, trait))
names(results) <- my_tissues
# parse and create results table
x <- grep("SE", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- grep("MX", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
y <- y[!y%in%x]
z <- grep("A5", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
z <- z[!z%in%c(x,y)]
w <- grep("A3", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
w <- w[!w%in%c(x,y,z)]
a <- grep("RI", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
a <- a[!a%in%c(x,y,z,w)]
b <- grep("AF", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
b <- b[!b%in%c(x,y,z,w,a)]
c <- grep("AL", colnames(t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))), value = T)
c <- c[!c%in%c(x,y,z,w,a,b)]
column_order <- c(x,y,z,w,a,b,c)

df1 <- cbind.data.frame(sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                        t(sapply(my_tissues, function(tissue) results[[tissue]]$fdr))[, column_order]
) 
df2 <- cbind.data.frame( sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                         t(sapply(my_tissues, function(tissue) results[[tissue]]$effect_sizes))[, column_order]
)
df2[is.na(df1)] <- NA
df <- cbind.data.frame(df1, df2[,-1])
df <- df[df[,1] < 0.05, apply(df, 2, function(x) sum(!is.na(x))) > 0]
write.table(df,
            paste0(table_path, "Table_S4H.", trait, ".tab"),
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")

# BMI ----
trait <- "BMI"
# number of DSEs per type of splicing events
number_of_DSEs <- t(sapply(tissues, function(tissue) 
  sapply(splicing_events, function(event) 
    sum(dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05, "type"]==event)
  )))
# tissues with at least 3 DSEs of each type of splicing event
my_tissues <- tissues[apply(number_of_DSEs, 1, function(x) sum(x>2))==7]
# run analaysis
results <- lapply(my_tissues, function(tissue) run_analysis(tissue, trait))
names(results) <- my_tissues
# parse and create results table
column_order <- c(x,y,z,w,a,b,c)

df <- t(as.data.frame(c(sapply(my_tissues, function(tissue) results[[tissue]]$effect_sizes)[column_order,],
                       sapply(my_tissues, function(tissue) results[[tissue]]$kruskal_p_value),
                       sapply(my_tissues, function(tissue) results[[tissue]]$fdr)[column_order,]
)))
rownames(df) <- my_tissues
colnames(df) <- c(column_order, "Kruskal-Wallis p-value", column_order)
write.table(df,
            paste0(table_path, "Table_S4H.", trait, ".tab"),
            col.names = T, row.names = T,
            quote = F,
            sep = "\t")

