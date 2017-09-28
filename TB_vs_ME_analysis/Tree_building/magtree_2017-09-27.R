library(ape)
treefile <- read.tree("C:/Users/Alex/Dropbox/allMAGS.nwk")

#Get metadata

metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/Genome_info/MAG_phylogeny.csv", header = T)

tips <- treefile$tip.label
lakekey <- substr(metadata$GenomeName, start = 1, stop = 2)
tips <- gsub(".genes.fna", "", tips)
matching <- match(tips, metadata$IMGOID)
tip.lakes <- lakekey[matching]
lake.colors <- c()
lake.colors[which(tip.lakes == "ME")] <- "#8dd3c7"
lake.colors[which(tip.lakes == "TE")] <- "#ffffb3"
lake.colors[which(tip.lakes == "TH")] <- "#bebada"

# Plot the tree with colors by lake

plot(treefile, type = "fan", show.tip.label = F, main = "", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = lake.colors, cex = 2)

#Add another with tip labels so I can add phylum colors
plot(treefile, type = "fan", show.tip.label = T, main = "", edge.width = 3)
