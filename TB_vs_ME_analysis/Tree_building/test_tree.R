# Make test phylogenetic tree

library(ape)
MEfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/ME_MAGS.nwk")

#Get trait info
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/Metapathways_output/master_table.csv")
keep <- match(substr(MEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]


trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(MEfile, type = "fan", show.tip.label = F, main = "Mendota Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

TEfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/TBE_MAGS.nwk")
keep <- match(substr(TEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]
trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(TEfile, type = "fan", show.tip.label = F, main = "Trout Bog Epi Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

THfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/TBH_MAGS.nwk")
keep <- match(substr(THfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
nfix <- table[keep, grep("nitrogen.fixation", colnames(table))]
trait.colors <- c()
trait.colors[which(nfix == "0")] <- "red"
trait.colors[which(nfix == "1")] <- "green"
plot(THfile, type = "fan", show.tip.label = F, main = "Trout Bog Hypo Nitrogen Fixation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

# Repeat with sulfur oxidation

MEfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/ME_MAGS.nwk")

#Get trait info
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/Metapathways_output/master_table.csv")
keep <- match(substr(MEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
sox <- table[keep, grep("sulfi.e.oxidation", colnames(table))]
sox <- sox != "0"
sox <- rowSums(sox)


trait.colors <- c()
trait.colors[which(sox == "0")] <- "red"
trait.colors[which(sox == "1")] <- "green"
plot(MEfile, type = "fan", show.tip.label = F, main = "Mendota Sulfur Oxidation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

TEfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/TBE_MAGS.nwk")
keep <- match(substr(TEfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
sox <- table[keep, grep("sulfi.e.oxidation", colnames(table))]
sox <- sox != "0"
sox <- rowSums(sox)

trait.colors <- c()
trait.colors[which(sox == "0")] <- "red"
trait.colors[which(sox == "1")] <- "green"
plot(TEfile, type = "fan", show.tip.label = F, main = "Trout Bog Epi Sulfur Oxidation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

THfile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/TBH_MAGS.nwk")
keep <- match(substr(THfile$tip.label, start = 0, stop = 10), substr(table$X, start = 2, stop = 11))
#keep <- keep[which(is.na(keep) == F)]
sox <- table[keep, grep("sulfi.e.oxidation", colnames(table))]
sox <- sox != "0"
sox <- rowSums(sox)

trait.colors <- c()
trait.colors[which(sox == "0")] <- "red"
trait.colors[which(sox == "1")] <- "green"
plot(THfile, type = "fan", show.tip.label = F, main = "Trout Bog Hypo Sulfur Oxidation", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 2)

