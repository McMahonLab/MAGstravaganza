# CAn I use genes like OTUs?
# We'll start in the Trout Bog Hypo, where I think we'll see the most intersting trends

hypo.data <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Previous_MAG_analyses/hypo_metagenome_genesurvey.csv", row.names=1)
hypo.dates <- as.Date(substr(colnames(hypo.data), start=2, stop=12), format = "%Y.%m.%d")

groups <- factor(substr(colnames(hypo.data), start = 4, stop = 5), levels = c("07", "08", "09"))

# How uneven are gene counts?
hist(rowSums(hypo.data))
# Pretty uneven
# What are the most abundant genes?
abundances <- rowSums(hypo.data)
abundances <- abundances[order(abundances, decreasing = T)]
abundances[1:15]
# Lots of central metabolites - acetyl-Coa , succinate, ATP synthase, glutamate synthase, aconitase(TCA), inorganic pyrophosphatase (fatty acid synthesis), phosphoglyeracte kinase, G3P dehydrogenase, NADH dehydrogenase, pyruvate:ferredoxin, pyruvate phosphate dikinase, cytochrome c oxidase

# Let's run a PCoA using genes as "taxa"

library(vegan)
library(ggplot2)
x <- vegdist(t(hypo.data), method = "bray")
pcoa <- betadisper(x, groups)
scores <- scores(pcoa)

plot.pcoa <- data.frame(scores$sites, groups)
colnames(plot.pcoa) <- c("PCoA1", "PCoA2", "Year")

axis1 <- round(pcoa$eig[1]/sum(pcoa$eig), digits = 2)
axis2 <- round(pcoa$eig[2]/sum(pcoa$eig), digits = 2)

ggplot(data=plot.pcoa, aes(x = PCoA1, y = PCoA2, color = Year)) + geom_point() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + labs(title = "TBH Genes", x = paste("PCoA1 (", axis1, ")", sep = ""), y = paste("PCoA2 (", axis2, ")"))

# That's one wicked arching pattern. PcoA probably isn't a good bet here. Try DCA instead
dca <- decorana(x)
coords <- dca$rproj
plot.dca <- data.frame(dca$rproj[,1:2], groups)
ggplot(data=plot.dca, aes(x = DCA1, y = DCA2, color = groups)) + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1))
# slight clustering by years, but not definitive. Add Julian data instead
plot.dca$Date <- as.numeric(format(hypo.dates, format = "%j"))
ggplot(data=plot.dca, aes(x = DCA1, y = DCA2, color = Date, shape = groups)) + geom_point(size = 3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black")) + theme(axis.text.x = element_text(hjust = 0.5, size = 10, colour = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 10, hjust = 0.5, vjust = 0.1), panel.border = element_rect(colour = "black", fill=NA, size=1)) + scale_color_gradient2(low = "green", mid = "yellow", high = "red", midpoint = 225)

#Not really clustering by date, either

heatmap(as.matrix(x), Rowv = T, Colv = T)

#Nope, definitely no clustering by year or Julian date

#Is this because changes in gene abundances are miniscule?

# Plot the top 100 genes over time
library(reshape2)
top100 <- hypo.data[match(names(abundances)[1:100], rownames(hypo.data)), ]
top100 <- prop.table(as.matrix(top100), margin = 2)
top100 <- melt(top100)
ggplot(top100, aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity") + theme(legend.position = "none")

# very constant over time, but admittedly the most abundant genes are also very common central metabolite genes.
# Can I break these down by category?

genekey <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Previous_MAG_analyses/genekey2.csv", row.names=1)
long.hypo <- melt(prop.table(as.matrix(hypo.data), margin = 2))
colnames(long.hypo) <- c("Gene", "Sample", "Hits")
long.hypo$Category <- genekey$Kegg.Category[match(long.hypo$Gene, genekey$Product.Name)]

#First a sanity check - do the broad categories stay constant over time?
ggplot(long.hypo, aes(x = Sample, y = Hits, fill = Category)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Yep, doing good. Some slight variation but not much

#Now break it down by categories
ggplot(long.hypo[which(long.hypo$Category == "Aromatic degradation"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Aromatic Degradataion")

ggplot(long.hypo[which(long.hypo$Category == "Carbon fixation - photosynthetic organisms"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Carbon fixation - photosynthetic organisms")

ggplot(long.hypo[which(long.hypo$Category == "Carbon fixation - prokaryotic"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Carbon fixation - prokaryotic")

ggplot(long.hypo[which(long.hypo$Category == "Methane metabolism"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Methane metabolism")

ggplot(long.hypo[which(long.hypo$Category == "Nitrogen metabolism"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Nitrogen metabolism")

ggplot(long.hypo[which(long.hypo$Category == "Oxidative Phosphorylation"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Oxidative Phosphorylation")

ggplot(long.hypo[which(long.hypo$Category == "Sulfur metabolism"), ], aes(x = Sample, y = Hits, fill = Gene)) + geom_bar(stat = "identity", position = "fill") + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(title = "Sulfur metabolism")

# Most constant, although that 2007 November date really throws off prokaryotic carbon fixation and nitrogen metabolism. Nitrogen fixation appears to decrease steadily from November 2007

long.hypo[which(long.hypo$Category == "Nitrogen metabolism" & long.hypo$Sample == "X2007.11.14"), ]
