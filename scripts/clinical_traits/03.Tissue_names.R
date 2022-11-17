data <- readRDS("~/Documents/mn4/Jose/00_Data/Tissue_info.rds")
tissues <- data$tissue_ID
write.table(tissues, "~/Documents/mn4/Jose/00_Data/Tissues.txt", row.names = F, col.names = F, quote = F)
