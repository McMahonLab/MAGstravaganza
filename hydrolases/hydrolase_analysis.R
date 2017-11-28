# Analyze carbohydrate active enzymes in the MAGs

# I have 2 heatmaps to make - one of diff types of CAZy enzymes by order by lake, one of glucoside hydrolase families by order by lake

# All my data is in a file by MAG - read in, make table, append to a larger dataset, and delete the object

# Read in list of MAGs

mag_data <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)

order <- sapply(strsplit(mag_data$Taxonomy, ";"), '[', 3)

files <- mag_data$IMG_OID
start <- read.table(paste("C:/Users/Alex/Desktop/MAGstravaganza/hydrolases/", files[1], ".txt", sep = ""))
start$V2 <- sub(".hmm", "", start$V2)
types <- substr(start$V2, start = 1, stop = 2)
type_data <- table(types)
gh_data <- table(start$V2[which(types == "GH")])

start_MAG_cazy <- rep(files[1], length(type_data))
lake_cazy <- rep(mag_data$Lake[1], length(type_data))
tax_cazy <- rep(order[1], length(type_data))
cazy_dataset <- data.frame(start_MAG_cazy, lake_cazy, tax_cazy, type_data)
colnames(cazy_dataset) <- c("MAG", "Lake", "Order", "Enzyme", "Freq")

start_MAG_gh <- rep(files[1], length(gh_data))
lake_gh <- rep(mag_data$Lake[1], length(gh_data))
tax_gh <- rep(order[1], length(gh_data))
gh_dataset <- data.frame(start_MAG_gh, lake_gh, tax_gh, gh_data)
colnames(gh_dataset) <- c("MAG", "Lake", "Order", "GH_Family", "Freq")

for(i in 2:dim(mag_data)[1]){
  start <- read.table(paste("C:/Users/Alex/Desktop/MAGstravaganza/hydrolases/", files[i], ".txt", sep = ""))
  start$V2 <- sub(".hmm", "", start$V2)
  types <- substr(start$V2, start = 1, stop = 2)
  type_data <- table(types)
  gh_data <- table(start$V2[which(types == "GH")])
  
  start_MAG_cazy <- rep(files[i], length(type_data))
  lake_cazy <- rep(mag_data$Lake[i], length(type_data))
  tax_cazy <- rep(order[i], length(type_data))
  partial_dataset <- data.frame(start_MAG_cazy, lake_cazy, tax_cazy, type_data)
  colnames(partial_dataset) <- c("MAG", "Lake", "Order", "Enzyme", "Freq")
  cazy_dataset <- rbind(partial_dataset, cazy_dataset)
  
  start_MAG_gh <- rep(files[i], length(gh_data))
  lake_gh <- rep(mag_data$Lake[i], length(gh_data))
  tax_gh <- rep(order[i], length(gh_data))
  partial_dataset <- data.frame(start_MAG_gh, lake_gh, tax_gh, gh_data)
  colnames(partial_dataset) <- c("MAG", "Lake", "Order", "GH_Family", "Freq")
  gh_dataset <- rbind(partial_dataset, gh_dataset)
}

# Plot heatmaps

library(ggplot2)
library(cowplot)


ggplot(cazy_dataset[which(cazy_dataset$Lake == "Mendota"),], aes(x = Order, y = Enzyme, fill = Freq)) + geom_tile() + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 70) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
ggplot(cazy_dataset[which(cazy_dataset$Lake == "Trout Bog Epilimnion"),], aes(x = Order, y = Enzyme, fill = Freq)) + geom_tile() + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 50) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
ggplot(cazy_dataset[which(cazy_dataset$Lake == "Trout Bog Hypolimnion"),], aes(x = Order, y = Enzyme, fill = Freq)) + geom_tile() + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 175) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

ggplot(gh_dataset[which(gh_dataset$Lake == "Mendota" & gh_dataset$Order != "[Blank]"  & gh_dataset$Freq > 3), ], aes(x = Order, y = GH_Family, fill = Freq)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 25)

ggplot(gh_dataset[which(gh_dataset$Lake == "Trout Bog Epilimnion" & gh_dataset$Order != "[Blank]"  & gh_dataset$Freq > 3), ], aes(x = Order, y = GH_Family, fill = Freq)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 20)

ggplot(gh_dataset[which(gh_dataset$Lake == "Trout Bog Hypolimnion" & gh_dataset$Order != "[Blank]"  & gh_dataset$Freq > 3), ], aes(x = Order, y = GH_Family, fill = Freq)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_gradient2(low = "white", mid = "#6baed6", high = "#08306b", midpoint = 60)
