# Make panels of nitrogen and sulfur pathways by lake
library(ggplot2)
library(cowplot)
# Use data from Excel pathway analysis

table <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T, row.names = 2)

colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 2:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}

# Read in metadata info to find out which MAGs came from which lake

metadata <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

# Mendota Nitrogen
MEnitrogen <- table[which(table$Category == "Nitrogen Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
MEnitrogen_totals <- rowSums(MEnitrogen > 0)
MEnitrogen_df <- data.frame(names(MEnitrogen_totals), MEnitrogen_totals)
colnames(MEnitrogen_df) <- c("Pathway", "MAGs")
ggplot(data = MEnitrogen_df, aes(x = Pathway, y = MAGs)) + geom_bar(stat = "identity", fill = "#FF0000", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(19, 0))

#Trout Bog totals
TEnitrogen <- table[which(table$Category == "Nitrogen Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake != "Mendota")])) == F)]
TEnitrogen_totals <- rowSums(TB > 0)
TEnitrogen_df <- data.frame(names(TEnitrogen_totals), TEnitrogen_totals)
colnames(TEnitrogen_df) <- c("Pathway", "MAGs")
ggplot(data = TEnitrogen_df, aes(x = Pathway, y = MAGs)) + geom_bar(stat = "identity", fill = "#00A287", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(19, 0))

#8EEB00