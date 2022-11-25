#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

load("~/GTEx_v8/Raquel/Draft/Analysis/expression/expression_breadth_and_sharing/scripts/05.sup_figure_exprs_ranking.RData")

for(trait in traits){
 names(traits.rnk[[trait]]) <- as.character(1:length(tissues)) 
}
trait <- "Age"
parse_data <- function(trait){
  df <- do.call(rbind.data.frame,lapply(1:46, function(i)  cbind.data.frame(rep(i, length(traits.rnk[[trait]][[i]])),
                                             traits.rnk[[trait]][[i]])
  ))
  colnames(df) <- c("number_of_tissues", "rank")
  df$trait <- rep(trait, nrow(df))
  return(df)
} 

data <- do.call(rbind.data.frame, 
                lapply(traits, function(trait) parse_data(trait)))
data$trait <- factor(data$trait, levels = traits[c(2,3,1,4)], order = T)
data$number_of_tissues <- as.character(data$number_of_tissues)
data$number_of_tissues <- factor(data$number_of_tissues, levels = rev(as.character(1:length(tissues))), order = T)
head(data)

get_box_stats <- function(y, upper_limit = quantile(y)[3]  ){
  return(data.frame(
    y = 0.95*upper_limit,
    label = paste(
     "n = ", length(y), "\n"#,
      #"Mean =", round(mean(y), 2), "\n",
      #"Median =", round(median(y), 2), "\n"
    )
  ))
}
p <- ggplot(data,
       aes(x = rank,
           y = number_of_tissues,
           fill = trait,
           col = trait)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(~trait, scale = "free_x") +
  theme_bw() +
  scale_colour_manual(values = traits_cols[c(2,3,1,4)]) +
  scale_fill_manual(values = sapply(traits_cols[c(2,3,1,4)], function(i) alpha(i, 0.5))) +
  ylab("Number of tissues in which a gene is expressed") + xlab("Expression ranking of the tissue with the lowest FDR") +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  stat_summary(fun.data = get_box_stats, geom = "text",
             hjust = 0.5, vjust = 0.8, size = 2, col = "black") 

pdf("~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/Figure_S2A.expression_ranking_genes_DE.pdf",
    width = 9, height = 9)
p
dev.off()


