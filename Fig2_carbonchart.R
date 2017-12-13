# Make figure from carbon chart

table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv")
library(ggplot2)
library(cowplot)
library(reshape2)

#make a long format table - for each phylum and lake, count how many times each pathway appears.

phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)
class<- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 2)
order <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 3)

colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 3:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}
rTCA <- table[2,3:196]
rTCA[which(rTCA < 100)] <- 0
table[2,3:196] <- rTCA

table$Pathway <- gsub("(*)", "", table$Pathway)
ME.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
ME.MAGs$Pathways <- make.unique(table$Pathway)

ME.long <- melt(ME.MAGs)
ME.long$Pathways <- factor(ME.long$Pathways, levels = table$Pathway)
ME.long$variable <- factor(ME.long$variable, levels = colnames(table))
ME.long$Taxonomy <- phylum[match(ME.long$variable, metadata$IMG_OID)]
ME.long$Taxonomy[which(ME.long$Taxonomy == "[Blank]")] <- "Unclassified"
ME.long <- ME.long[which(is.na(ME.long$Pathways) == F), ]
ME.agg <- aggregate(value ~ Pathways + Taxonomy, ME.long, mean)

ME.heatmap <- ggplot(data = ME.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#FF0000", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/ME.heatmap.pdf", ME.heatmap, base_height = 7, base_aspect_ratio = 0.75)

TE.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TE.MAGs$Pathways <- make.unique(table$Pathway)

TE.long <- melt(TE.MAGs)
TE.long$Pathways <- factor(TE.long$Pathways, levels = table$Pathway)
TE.long$variable <- factor(TE.long$variable, levels = colnames(table))
TE.long$Taxonomy <- phylum[match(TE.long$variable, metadata$IMG_OID)]
TE.long$Taxonomy[which(TE.long$Taxonomy == "[Blank]")] <- "Unclassified"
TE.long <- TE.long[which(is.na(TE.long$Pathways) == F), ]
TE.agg <- aggregate(value ~ Pathways + Taxonomy, TE.long, mean)

TE.heatmap <- ggplot(data = TE.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#8EEB00", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/TE.heatmap.pdf", TE.heatmap, base_height = 7, base_aspect_ratio = 0.67)

TH.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
TH.MAGs$Pathways <- make.unique(table$Pathway)

TH.long <- melt(TH.MAGs)
TH.long$Pathways <- factor(TH.long$Pathways, levels = table$Pathway)
TH.long$variable <- factor(TH.long$variable, levels = colnames(table))
TH.long$Taxonomy <- phylum[match(TH.long$variable, metadata$IMG_OID)]
TH.long$Taxonomy[which(TH.long$Taxonomy == "[Blank]")] <- "Unclassified"
TH.long <- TH.long[which(is.na(TH.long$Pathways) == F), ]
TH.agg <- aggregate(value ~ Pathways + Taxonomy, TH.long, mean)

TH.heatmap <- ggplot(data = TH.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#00A287", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/TH.heatmap.pdf", TH.heatmap, base_height = 7, base_aspect_ratio = .73)
