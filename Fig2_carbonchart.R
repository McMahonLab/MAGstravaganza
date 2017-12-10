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


ggplot(data = ME.long, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#5ab4ac", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

TE.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TE.MAGs$Pathways <- make.unique(table$Pathway)

TE.long <- melt(TE.MAGs)
TE.long$Pathways <- factor(TE.long$Pathways, levels = table$Pathway)
TE.long$variable <- factor(TE.long$variable, levels = colnames(table))
TE.long$Taxonomy <- phylum[match(TE.long$variable, metadata$IMG_OID)]


ggplot(data = TE.long, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#67a9cf", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

TH.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
TH.MAGs$Pathways <- make.unique(table$Pathway)

TH.long <- melt(TH.MAGs)
TH.long$Pathways <- factor(TH.long$Pathways, levels = table$Pathway)
TH.long$variable <- factor(TH.long$variable, levels = colnames(table))
TH.long$Taxonomy <- phylum[match(TH.long$variable, metadata$IMG_OID)]


ggplot(data = TH.long, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#d8b365", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
##############
lake <- c(rep("ME", 26), rep("TB", 52))
label <- c(rep("Mendota", 26), rep("TroutEpi", 26), rep("TroutHypo", 26))
carbon <- rep(colnames(table)[4:29], 3)
counts <- c(ME.total, TE.total, TH.total)

to.plot <- data.frame(lake, label, carbon, counts)
to.plot$carbon <- factor(to.plot$carbon, levels = rev(colnames(table)[4:29]))

ggplot(data = to.plot, aes(x = label, y = carbon, fill = counts)) + geom_tile() + labs(x = "", y = "") + scale_fill_gradient2(low = "#f7fcfd", mid = "#9ebcda", high = "#88419d", midpoint = 0.25) + theme(axis.text.x = element_text(angle = 270))

ME.cyano <- colSums(table[which(table$Lake == "ME" & table$phylum == "Cyanobacteria"), 4:29])/length(which(table$Lake == "ME" & table$phylum == "Cyanobacteria"))
TB.chlorobi <- colSums(table[which(table$Lake == "TE" & table$phylum == "Chlorobi" | table$Lake == "TH" & table$phylum == "Chlorobi"), 4:29])/length(which(table$Lake == "TE" & table$phylum == "Chlorobi" | table$Lake == "TH" & table$phylum == "Chlorobi")) * -1
ME.burk <- colSums(table[which(table$Lake == "ME" & table$order == "Burkholderiales"), 4:29])/length(which(table$Lake == "ME" & table$order == "Burkholderiales"))
TB.burk <- colSums(table[which(table$Lake == "TE" & table$order == "Burkholderiales" | table$Lake == "TH" & table$order == "Burkholderiales"), 4:29])/length(which(table$Lake == "TE" & table$order == "Burkholderiales" | table$Lake == "TH" & table$order == "Burkholderiales")) * -1
ME.bact <- colSums(table[which(table$Lake == "ME" & table$phylum == "Bacteroidetes"), 4:29])/length(which(table$Lake == "ME" & table$phylum == "Bacteroidetes"))
TB.bact <- colSums(table[which(table$Lake == "TE" & table$phylum == "Bacteroidetes" | table$Lake == "TH" & table$phylum == "Bacteroidetes"), 4:29])/length(which(table$Lake == "TE" & table$phylum == "Bacteroidetes" | table$Lake == "TH" & table$phylum == "Bacteroidetes")) * -1
ME.ver <- colSums(table[which(table$Lake == "ME" & table$phylum == "Verrucomicrobia"), 4:29])/length(which(table$Lake == "ME" & table$phylum == "Verrucomicrobia"))
TB.ver <- colSums(table[which(table$Lake == "TE" & table$phylum == "Verrucomicrobia" | table$Lake == "TH" & table$phylum == "Verrucomicrobia"), 4:29])/length(which(table$Lake == "TE" & table$phylum == "Verrucomicrobia" | table$Lake == "TH" & table$phylum == "Verrucomicrobia"))* -1
ME.act <- colSums(table[which(table$Lake == "ME" & table$phylum == "Actinobacteria"), 4:29])/length(which(table$Lake == "ME" & table$phylum == "Actinobacteria"))
TB.act <- colSums(table[which(table$Lake == "TE" & table$phylum == "Actinobacteria" | table$Lake == "TH" & table$phylum == "Actinobacteria"), 4:29])/length(which(table$Lake == "TE" & table$phylum == "Actinobacteria" | table$Lake == "TH" & table$phylum == "Actinobacteria")) * -1
ME.methT <- colSums(table[which(table$Lake == "ME" & table$order == "Methylophilales"), 4:29])/length(which(table$Lake == "ME" & table$order == "Methylophilales"))
TB.methT <- colSums(table[which(table$Lake == "TE" & table$order == "Methylophilales" | table$Lake == "TH" & table$order == "Methylophilales"), 4:29])/length(which(table$Lake == "TE" & table$order == "Methylophilales" | table$Lake == "TH" & table$order == "Methylophilales")) * -1
ME.methC <- colSums(table[which(table$Lake == "ME" & table$order == "Methylococcales"), 4:29])/length(which(table$Lake == "ME" & table$order == "Methylococcales"))
TB.methC <- colSums(table[which(table$Lake == "TE" & table$order == "Methylococcales" | table$Lake == "TH" & table$order == "Methylococcales"), 4:29])/length(which(table$Lake == "TE" & table$order == "Methylococcales" | table$Lake == "TH" & table$order == "Methylococcales")) * -1
ME.plan <- colSums(table[which(table$Lake == "ME" & table$phylum == "Planctomycetes"), 4:29])/length(which(table$Lake == "ME" & table$phylum == "Planctomycetes"))
TB.acid <- colSums(table[which(table$Lake == "TE" & table$phylum == "Acidobacteria" | table$Lake == "TH" & table$phylum == "Acidbacteria"), 4:29])/length(which(table$Lake == "TE"  & table$phylum == "Acidobacteria"| table$Lake == "TH" & table$phylum == "Acidobacteria")) * -1



lake <- c(rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26), rep("ME", 26), rep("TB", 26))
label <- c(rep("Cyanobacteria", 26), rep("Burkholderiales-ME", 26), rep("Bacteroidetes-ME", 26), rep("Verrucomicrobia-ME", 26), rep("Actinobacteria-ME", 26), rep("Methylophilales-ME", 26), rep("Methylococcales-ME", 26), rep("Planctomycetes", 26), rep("Chlorobi", 26), rep("Burkholderiales-TB", 26), rep("Bacteroidetes-TB", 26), rep("Verrucomicrobia-TB", 26), rep("Actinobacteria-TB", 26), rep("Methylophilales-TB", 26), rep("Methylococcales-TB", 26), rep("Acidobacteria", 26))
carbon <- rep(colnames(table)[4:29], 16)
counts <- c(ME.cyano, ME.burk, ME.bact, ME.ver, ME.act, ME.methT, ME.methC, ME.plan, TB.chlorobi, TB.burk, TB.bact, TB.ver, TB.act, TB.methT, TB.methC, TB.acid)

to.plot <- data.frame(lake, label, carbon, counts)
to.plot$carbon <- factor(to.plot$carbon, levels = rev(colnames(table)[4:29]))
to.plot$label <- factor(to.plot$label, levels = c("Cyanobacteria", "Chlorobi", "Burkholderiales-ME", "Burkholderiales-TB", "Bacteroidetes-ME", "Bacteroidetes-TB", "Verrucomicrobia-ME", "Verrucomicrobia-TB", "Actinobacteria-ME", "Actinobacteria-TB", "Methylophilales-ME", "Methylophilales-TB", "Methylococcales-ME", "Methylococcales-TB", "Planctomycetes", "Acidobacteria"))

v1 <- ggplot(data = to.plot, aes(x = label, y = carbon, fill = counts)) + geom_tile() + labs(x = "", y = "") + scale_fill_gradient2(low = "white", mid = "#9ebcda", high = "#88419d", midpoint = 0.4) + theme(axis.text.x = element_text(angle = 270))

v2 <- ggplot(data = to.plot, aes(x = label, y = carbon, fill = counts)) + geom_tile() + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#5ab4ac", midpoint = 0.4) + theme(axis.text.x = element_text(angle = 270))

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/fig2_me.pdf", v2, base_height = 10, base_width = 7)
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/fig2_tb.pdf", v1, base_height = 10, base_width = 7)

v3 <- ggplot(data = to.plot, aes(x = label, y = carbon, fill = counts)) + geom_tile() + labs(x = "", y = "") + scale_fill_gradient2(low = "#d8b365", mid = "white", high = "#5ab4ac", midpoint = 0) + theme(axis.text.x = element_text(angle = 270))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/fig2.pdf", v3, base_height = 10, base_width = 7)
