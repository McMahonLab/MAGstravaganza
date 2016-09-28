# Make test phylogenetic tree

library(ape)
file <- read.tree(file.choose())

trait <- runif(9, 0, 1)
trait.colors <- c()
trait.colors[which(trait >= 0.5)] <- "green"
trait.colors[which(trait < 0.5)] <- "red"


plot(file, type = "fan", show.tip.label = F, main = "TBH Test Tree", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = trait.colors, cex = 3)
