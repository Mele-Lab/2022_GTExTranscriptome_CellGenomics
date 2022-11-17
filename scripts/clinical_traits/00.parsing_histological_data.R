#Loading libraries
library(optparse)
 
#Reading arguments
# Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-f", "--file"), type="character",
                     dest="file",
                     help="File downloaded from the GTEx Portal: histological data, 'GTEX Portal.csv'")
options=parse_args(parser)

file=options$file # file <- "GTEx Portal.csv"

# Parsing data into the format we want: From one variable with all the diseases (comma separated) to one variable per disease (levels 0=healthy and 1=disease)
data <- read.csv(file)
df1 <- na.omit(stack(setNames(strsplit(data$Pathology.Categories, ","), seq_len(nrow(data))))[, 2:1])
df1$values <- gsub(" ","",df1$values)
tab <- as.data.frame.matrix(table(df1))
final <- cbind(data, tab)

#Updating the variable spermatogenesis to our own definition: Reduced spermatogenesis
spermatogenesis <- final[final$spermatogenesis==1,]
spermatogenesis$spermatogenesis <- 0
spermatogenesis[grep("reduce|redueced|reduction|decrease|diminish|arrest|absent|atrophy", spermatogenesis$Pathology.Notes),]$spermatogenesis <- 1
#I annotated atrophy as reduced spermatogenesis as 90 % of the samples annotated as atrophy are also annotated as reduced, and by direct visual inspection of the histology images, they look as reduced spermatogenesis as well
#Autolysis slightly correlates with reduced spermatogenesis too, but I am leaving it as it is as by image inspection I don't see it clearly
final <- final[,-62]

final_2 <- merge(final, spermatogenesis, by=colnames(final), all.x = TRUE)
final_2[final_2$monckeberg==1,"sclerotic"] <- 1  #Monckeberg samples are annotated as sclerotic

dir <- gsub("/[^/]+$","",file)
write.csv(final_2, paste0(dir,"/histological_data.csv"), na="NA")

