#!/usr/bin/env Rscript

path_to_data <- "~/GTEx_v8/Raquel/Draft/Cell_submission/data/"
setwd(path_to_data)
plot_path <- "~/GTEx_v8/Raquel/Draft/Cell_submission/Figures/"
if(!dir.exists(plot_path)){dir.create(plot_path, recursive = T)}

# libraries ----
library(ggplot2)


# metadata ----
load("00.metadata.RData")

# differential splicing analysis ----
dsa <- readRDS("00.differential_splicing_analysis.rds")

# Figure 4C ----
ASEs_isoforms_biotype <- lapply(tissues, function(tissue) paste0(rownames(dsa[[tissue]][["Age"]]),
                                                                 ":",
                                                                 dsa[[tissue]][["Age"]]$spliced_in, 
                                                                 "",
                                                                 dsa[[tissue]][["Age"]]$spliced_out, 
                                                                 "_",
                                                                 dsa[[tissue]][["Age"]]$biotype))
names(ASEs_isoforms_biotype) <- tissues
table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]]))
sum(table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]])))

ASEs_with_domain <- lapply(tissues, function(tissue)
  unique(c(rownames(dsa[[tissue]][["Age"]][!is.na(dsa[[tissue]][["Age"]]$spliced_in_domains),]),
                  rownames(dsa[[tissue]][["Age"]][!is.na(dsa[[tissue]][["Age"]]$spliced_out_domains),]))))
names(ASEs_with_domain) <- tissues
events_domain <- lapply(tissues, function(tissue) 
  paste0(rownames(dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]],]),
         ":",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]], "spliced_in"], 
         "",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]], "spliced_out"], 
         "_",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]],"biotype"])
  )
names(events_domain) < tissues
table(sapply(unique(unlist(events_domain)), function(i) unlist(strsplit(i, split = "_"))[[2]]))

m <- matrix(round(c(c(0,10067),
                    c(18013,13466),
                    c(35656-7849, 7849))), 2,3)



pdf(paste0(plot_path, "Figure_4C.pdf"),
    width = 3, height = 6)
barplot(m, yaxt = 'n',
        ylab = "Alternatively spliced events",
        border = NA,
        horiz = F)
axis(2, at = axTicks(2), labels = prettyNum(axTicks(2), big.mark = ","), las = 2)
dev.off()

pdf(paste0(plot_path, "Figure_4C.legend.pdf"),
    width = 6, height = 4)
plot.new()
legend("center",
       title = "switch between isoforms",
       c("non-coding isoforms",
         "non-coding (spliced-in) and protein-coding (spliced-out)",
         "protein-coding (spliced-in) and non-coding (spliced-out)",
         "protein-coding isoforms",
         "AS events affecs a PFAM domain"),
       bty = 'n',
       pch = 15,
       col = c("#e6e6e6ff", "#91C483", "#FF6464", "#FFE162", "#ccb44e")
)
dev.off()


biotype_cols <- rev(c("#EEEEEE", "#91C483", "#FF6464", "#FFE162"))
names(biotype_cols) <- rev(c("NC-NC", "NC-PC", "PC-NC", "PC-PC"))
sum(table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]])))
table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]]))["PC-PC"]/
  length(unique(unlist(ASEs_isoforms_biotype)))
table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]]))[c("NC-PC", "PC-NC")]/
  length(unique(unlist(ASEs_isoforms_biotype)))
sum(table(sapply(unique(unlist(ASEs_isoforms_biotype)), function(i) unlist(strsplit(i, split = "_"))[[2]]))[c("NC-PC", "PC-NC")]/
      length(unique(unlist(ASEs_isoforms_biotype))))

# Figure 4D ----
for(tissue in tissues){
  dsa[[tissue]]$Age$Type <- sapply(rownames(dsa[[tissue]]$Age), function(i) unlist(strsplit(unlist(strsplit(i, split = ";"))[[2]], split = ":"))[[1]] )
}

df <- rbind.data.frame("Tissue" = rep(tissues, 2),
                       "variable" = c(rep("NC-PC", length(tissues)),
                                      rep("PC-NC", length(tissues))), 
                       "type" = unlist(lapply(splicing_events, function(e) rep(e, 2*length(tissues)))),
                       "value" = unlist(lapply(splicing_events, function(e)
                         c(sapply(tissues, function(tissue) table(dsa[[tissue]][["Age"]][dsa[[tissue]][["Age"]]$Type==e, "biotype"])["NC-PC"]/
                                    nrow(dsa[[tissue]][["Age"]][dsa[[tissue]][["Age"]]$Type==e &
                                                                      dsa[[tissue]][["Age"]]$biotype %in% c("NC-PC", "PC-NC"), ]) ),
                           sapply(tissues, function(tissue) table(dsa[[tissue]][["Age"]][dsa[[tissue]][["Age"]]$Type==e, "biotype"])["PC-NC"]/
                                    nrow(dsa[[tissue]][["Age"]][dsa[[tissue]][["Age"]]$Type==e &
                                                                      dsa[[tissue]][["Age"]]$biotype %in% c("NC-PC", "PC-NC"), ]) ))))
)

df <- as.data.frame(t(df))                        
colnames(df) <- c("tissue", "variable", "type","value")
df$type <- factor(df$type, levels = splicing_events, order = T)
df$variable <- factor(df$variable, levels = c("NC-PC", "PC-NC"), order = T)
df$value <- as.numeric(df$value)

pdf(paste0(plot_path, "Figure_4D.pdf"),
    width = 4, height = 4)
ggplot(df,
       aes( x = type,
            y = value,
            col = variable)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_manual(values = c("#91C483", "#FF6464")) +
  theme_bw()+
  ylab("Alternatively spliced events (%)") +
  xlab("") +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

for(tissue in tissues){
  for(trait in traits){
    dsa[[tissue]][[trait]] <- dsa[[tissue]][[trait]][dsa[[tissue]][[trait]]$adj.P.Val < 0.05,]
  }
}

ASEs_with_domain <- lapply(tissues, function(tissue)
  unique(unlist(lapply(traits, function(trait) 
  unique(c(rownames(dsa[[tissue]][[trait]][!is.na(dsa[[tissue]][[trait]]$spliced_in_domains),]),
           rownames(dsa[[tissue]][[trait]][!is.na(dsa[[tissue]][[trait]]$spliced_out_domains),])))))))
names(ASEs_with_domain) <- tissues
events_domain <- lapply(tissues, function(tissue) 
  paste0(rownames(dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]],]),
         ":",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]], "spliced_in"], 
         "",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]], "spliced_out"], 
         "_",
         dsa[[tissue]][["Age"]][ASEs_with_domain[[tissue]],"biotype"])
)
names(events_domain) < tissues
table(sapply(unique(unlist(events_domain)), function(i) unlist(strsplit(i, split = "_"))[[2]]))
