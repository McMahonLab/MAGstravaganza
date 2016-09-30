# Make test phylogenetic tree

library(ape)
MEfile <- read.tree(file.choose())

#Get trait info
table <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Data_files/Metapathways_output/master_table.csv")
keep <- match(substr(MEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]


trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(file, type = "fan", show.tip.label = F, main = "Mendota Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

TEfile <- read.tree(file.choose())
keep <- match(substr(TEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]
trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(TEfile, type = "fan", show.tip.label = F, main = "Trout Bog Epi Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

THfile <- read.tree(file.choose())
keep <- match(substr(THfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]
trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(THfile, type = "fan", show.tip.label = F, main = "Trout Bog Hypo Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

