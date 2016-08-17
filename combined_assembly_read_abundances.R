#Goal: Compare phylogeny differences in TB vs ME by read annotations

library(ggplot2)
library(scales)
path2repo <- "C:/Users/Alex/Desktop/MAGstravaganza/"

TBH <- read.csv(paste(path2repo, "TBhypo_IMGblast_phylogeny.csv", sep = ""), header = T)
TBE <- read.csv(paste(path2repo, "TBepi_IMGblast_phylogeny.csv", sep = ""), header = T)
ME <- read.csv(paste(path2repo, "ME_IMGblast_phylogeny.csv", sep = ""), header = T)

#Plot the top 10 most abundant phyla at 60% identity for each lake

TBHabun <- TBH[,c(2,5)]
TBHabun <- TBHabun[order(TBHabun[,2], decreasing = T),]
TBHtop <- TBHabun[c(2:11),] #Not including "unassigned" in this measure
colnames(TBHtop) <- c("Phylum", "Hits")
TBHtop$Phylum <- factor(TBHtop$Phylum, levels = TBHtop$Phylum)

TBEabun <- TBE[,c(2,5)]
TBEabun <- TBEabun[order(TBEabun[,2], decreasing = T),]
TBEtop <- TBEabun[c(2:11),] #Not including "unassigned" in this measure
colnames(TBEtop) <- c("Phylum", "Hits")
TBEtop$Phylum <- factor(TBEtop$Phylum, levels = TBEtop$Phylum)

TBHtop$Layer <- c(rep("TBH", dim(TBHtop)[1]))
TBEtop$Layer <- c(rep("TBE", dim(TBEtop)[1]))
TBtop <- rbind(TBEtop, TBHtop)

pdf(file = paste(path2repo, "TBabundance.pdf"), height = 6, width = 12)
ggplot(data = TBtop, aes(x = Phylum, y = Hits, fill = Layer)) + geom_bar(stat = "identity", position = "dodge" ) + theme_bw() + labs(x = NULL, y = "Hits at 60%+ BLAST Identity", title = "Trout Bog Read Abundances") + theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5, size = 28), axis.text.y = element_text(size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 32)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 100000), labels = comma) + scale_fill_manual(values = c("darkred", "dodgerblue"))
dev.off()

MEabun <- ME[,c(2,5)]
MEabun <- MEabun[order(MEabun[,2], decreasing = T),]
MEtop <- MEabun[c(2:11),] #Not including "unassigned" in this measure
colnames(MEtop) <- c("Phylum", "Hits")
MEtop$Phylum <- factor(MEtop$Phylum, levels = MEtop$Phylum)

pdf(file = paste(path2repo, "MEabundance.pdf"), height = 6, width = 12)
ggplot(data = MEtop, aes(x = Phylum, y = Hits)) + geom_bar(stat = "identity", fill = c("darkred")) + theme_bw() + labs(x = NULL, y = "Hits at 60%+ BLAST Identity", title = "Lake Mendota Read Abundances") + theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust= 0.5, size = 24), axis.text.y = element_text(size = 20), axis.title = element_text(size = 20), plot.title = element_text(size = 28)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 575000), labels = comma)
dev.off()

