#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)
library(ggtext)

# metadata ----
load("00.metadata.RData")


# Figure 4A ----
# DSEs tissue sharing --
tissue_sharing <- readRDS("events_DS.tissue_sharing.rds")
for(trait in traits){
  tissue_sharing[[trait]]$trait <- rep(trait, nrow(tissue_sharing[[trait]]))
}

# Table S3E and S3F --
ribosomal_proteins <- gene_annotation[c(grep("^RPS", gene_annotation$gene.name, perl = T),
                                        grep("^RPL", gene_annotation$gene.name, perl = T)),]

# tissue-shared genes distribution --
traits_cols <- c(traits_cols, "black")
names(traits_cols)[5] <- "gene"

d0 <- do.call(rbind.data.frame, lapply(traits, function(trait) tissue_sharing[[trait]][, c(1,2,3,4,9)]))
color_variable <- function(gene, trait){
  if(gene %in% ribosomal_proteins$gene.name & trait == "Ancestry"){
    "gene"
  }else{
    trait
  }
}
d0$col <- sapply(1:nrow(d0), function(row) color_variable(d0[row, "gene_name"], d0[row, "trait"]))
d0$trait <- factor(d0$trait, levels = c(traits, "gene"), order = T)
d0 <- rbind.data.frame(d0[d0$col != "gene",], d0[d0$col == "gene",])
# so the plot is not so heavy
df <- rbind.data.frame(do.call(rbind.data.frame, 
                               lapply(traits[-4], function(trait) 
                                 do.call(rbind.data.frame, lapply(1:2, function(i) d0[d0$trait == trait & d0$n_DE == i,][sample(nrow(d0[d0$trait == trait & d0$n_DE == i,]), 500),]))
                               )),
                       d0[d0$n_DE > 4 & d0$trait %in% traits[c(1:3)],],
                       d0[d0$trait == "BMI",])
df <- d0
# plot
p1 <- ggplot(df,
             aes(x = trait,
                 y = n_DS,
                 col = col)) +
  geom_jitter() +
  scale_color_manual(values = traits_cols) +
  theme_bw() +
  ylab("Number of tissues") + xlab("") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = 10, lty = 2, size = 0.25) +
  geom_hline(yintercept = 1.5, lty = 2, size = 0.25) +
  geom_hline(yintercept = 5.5, lty = 2, size = 0.25) +
  geom_text_repel(data = rbind.data.frame(df[df$col == "gene" & df$trait == "Ancestry",]),#,
                                          #df[df$trait == "Sex",][1:3,],
                                          #df[df$trait == "Age",][1:3,],
                                          #df[df$trait == "BMI",][1:3,]), 
                  aes(x = trait, y = n_DS, label = gene_name), 
                  cex = 3.5, col = "black")

pdf(paste0(plot_path, "Figure_4A.tissue_sharing.pdf"),
    width = 6, height = 6)
p1
dev.off()

# proportion of events with a given tissue-sharing per trait --
total_DSEs <- sapply(traits, function(trait) nrow(tissue_sharing[[trait]]))
d <- rbind.data.frame(sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DS==1)/total_DSEs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DS>=2 & tissue_sharing[[trait]]$n_DS<=5)/total_DSEs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DS>=6 & tissue_sharing[[trait]]$n_DS<=9)/total_DSEs[trait])
),
sapply(traits, function(trait) 
  100*(sum(tissue_sharing[[trait]]$n_DS>=10 )/total_DSEs[trait])
))
colnames(d) <- traits
rownames(d) <- c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")


df2 <- melt(d)
df2$type <- rep(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing"), 4)
df2$type <- factor(df2$type, levels = rev(c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")), order = T) 
cols <- brewer.pal(4, "Greys")
names(cols) <- c("Tissue-specific", "Low sharing", "Moderate sharing", "High sharing")
p2 <- ggplot(df2, aes(x = 1,
                      y = value,
                      fill = type)) +
  geom_bar(stat= "identity") + 
  coord_flip() +
  theme_bw() +
  facet_grid(~variable) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  xlab("") + ylab("DSEs (%)")


pdf(paste0(plot_path, "Figure_4A.tissue_sharing.bar_plot.pdf"),
    width = 12, height = 1.25)
p2
dev.off()


