library(ape)
treefile <- read.tree("C:/Users/Alex/Dropbox/allMAGS.nwk")

#Get metadata

metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T, colClasses = c("character", "character", "factor", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))


#load n2 fixer data
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T)
colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 3:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}
hasN2 <- table[6, ] > 0


tips <- treefile$tip.label
tips <- gsub(".genes.fna", "", tips)

lakes <- metadata$Lake[match(tips, metadata$IMG_OID)]
Nfix <- hasN2[match(tips, colnames(table))]
key <- data.frame(tips, lakes, Nfix)

lake.colors <- c()
lake.colors[which(key$lakes == "Mendota" & key$Nfix == TRUE)] <- "#FF0000"
lake.colors[which(key$lakes == "Trout Bog Epilimnion" & key$Nfix == TRUE)] <- "#8EEB00"
lake.colors[which(key$lakes == "Trout Bog Hypolimnion" & key$Nfix == TRUE)] <- "#00A287"
lake.colors[which(is.na(key$lakes) == T | key$Nfix == FALSE)] <- "white"


m1 <- match(colnames(table)[3:196], metadata$IMG_OID)
sorted_metadata <- metadata[m1[which(is.na(m1) == F)],]

key <- paste(sorted_metadata$Lake, hasN2[3:196], sep = ":")


# Plot the tree with colors by lake for N2 fixers

plot(treefile, type = "fan", show.tip.label = F, main = "", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = lake.colors, cex = 2)

#Add another with tip labels so I can add phylum colors
plot(treefile, type = "fan", show.tip.label = T, main = "", edge.width = 3)
