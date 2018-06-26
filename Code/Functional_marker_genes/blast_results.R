# Does BLASTing reveal any assembly bias?

library(reshape2)
library(ggplot2)
library(cowplot)
# To make sure that the pathway content I'm seeing in the MAGs is representative of the ecosystem, I'm checking for genes in the metagenomes. If I get a lot of hits across lakes, then the gene is probably present in both and didn't assemble in one of them for some reason.

# I'm using BLASTx (nucleotide vs protein) on marker genes I've selected for pathways mentioned in the manuscript. I used IMG's function search to look for these marker genes in the MAGs and downloaded any gene with those functions from both lakes. This resulted in about 1000 genes representing ~30 pathways. I then BLASTed all metagenomes against this list of genes.

# The search terms used are in MAGstravaganza/Data_files/marker_genes_to_blast.csv

# Load BLAST results. I'm only looking at things with a percent ID greater than 80%.

blast_results <- read.table(file = "D:/pid80_results.txt", header = F)
colnames(blast_results) <- c("GeneID", "PID", "metagenome")
blast_results$metagenome <- gsub(".len150", "", blast_results$metagenome)
blast_results$metagenome <- gsub("IHSC1", "IHSC2", blast_results$metagenome)

# Load the gene names that go with the IDs
fasta_headers <- read.csv(file = "D:/fastaheaders_to_blast.csv", header = T)

# Load metagenome size info
mendota <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
tbe <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
tbh <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")

all_lengths <- rbind(mendota[,c(1,2)], tbe[,c(1,2)], tbh[,c(1,2)])

# Add column of lake for the metagenome
lakekey <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metagenome_metadata.csv", header = T)
blast_results$Lake <- lakekey$site[match(blast_results$metagenome, lakekey$sample)]

counts <- table(GeneID = blast_results$GeneID, blast_results$metagenome)
counts <- melt(counts)
counts <- counts[which(counts$value > 0),]

#I need to normalize counts by metagenome size and gene size.
counts$GeneLength <- fasta_headers$Length[match(counts$GeneID, fasta_headers$GeneID)]
counts$metaSize <- all_lengths$reads[match(counts$Var2, all_lengths$sample)]

# Convert to RPKM - general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
counts$RPKM <- counts$value / (counts$GeneID/1000 * counts$metaSize/1000000)

# Add the lake info for each metagenome
counts$Lake <- lakekey$site[match(counts$Var2, lakekey$sample)]

# Add product info
counts$Product <- fasta_headers$Product[match(counts$GeneID, fasta_headers$GeneID)]

# The hard part - linking product back to pathway

pathway <- rep(NA, dim(counts)[1])
hits <- grep("carboxylase", counts$Product)
pathway[hits] <- "CBB Cycle"
hits <- grep("lyase", counts$Product)
pathway[hits] <- "rTCA Cycle"
hits <- grep("rhodopsin", counts$Product)
pathway[hits] <- "Rhodopsin"
hits <- grep("aa3", counts$Product)
pathway[hits] <- "Cytochrome Oxidase aa3-type"
hits <- grep("cbb3", counts$Product)
pathway[hits] <- "Cytochrome Oxidase cbb3-type"
hits <- grep("Cbb3", counts$Product)
pathway[hits] <- "Cytochrome Oxidase cbb3-type"
hits <- grep("oxidase d", counts$Product)
pathway[hits] <- "Cytochrome Oxidase d"
hits <- grep("Fe-hydrogenase III", counts$Product)
pathway[hits] <- "Ni,Fe-hydrogenase III"
hits <- grep("Fe-hydrogenase I ", counts$Product)
pathway[hits] <- "Ni,Fe-hydrogenase I or II"
hits <- grep("Fe-hydrogenase 1 ", counts$Product)
pathway[hits] <- "Ni,Fe-hydrogenase I or II"
hits <- grep("Fe-hydrogenase 2 ", counts$Product)
pathway[hits] <- "Ni,Fe-hydrogenase I or II"
hits <- grep("rhamnulose", counts$Product)
pathway[hits] <- "Rhamnose Degradation"
hits <- grep("fuculose", counts$Product)
pathway[hits] <- "Fucose Degradation"
hits <- grep("Fuculose", counts$Product)
pathway[hits] <- "Fucose Degradation"
hits <- grep("ribulokinase", counts$Product)
pathway[hits] <- "Arabinose Degradation"
hits <- grep("aldose 1-epimerase", counts$Product)
pathway[hits] <- "Galactose Degradation"
hits <- grep("Aldose 1-epimerase", counts$Product)
pathway[hits] <- "Galactose Degradation"
hits <- grep("mannose-6-phosphate", counts$Product)
pathway[hits] <- "Mannose Degradation"
hits <- grep("xylose isomerase", counts$Product)
pathway[hits] <- "Xylose Degradation"
hits <- grep("methanol", counts$Product)
pathway[hits] <- "Methanol Degradation"
hits <- grep("nitrogenase", counts$Product)
pathway[hits] <- "Nitrogen Fixation"
hits <- grep("Nitrogenase", counts$Product)
pathway[hits] <- "Nitrogen Fixation"
hits <- grep("putrescine", counts$Product)
pathway[hits] <- "Putrescine Degradation"
hits <- grep("nitrite", counts$Product)
pathway[hits] <- "Dissimilatory Nitrate Reduction"
hits <- grep("nitric oxide", counts$Product)
pathway[hits] <- "Denitrification"
hits <- grep("quinone", counts$Product)
pathway[hits] <- "Sulfide Oxidation"
hits <- grep("flavocytochrome", counts$Product)
pathway[hits] <- "Sulfide Oxidation"
hits <- grep("Flavocytochrome", counts$Product)
pathway[hits] <- "Sulfide Oxidation"
hits <- grep("sulfite reductase", counts$Product)
pathway[hits] <- "Sulfide Reduction"
hits <- grep("Sox", counts$Product)
pathway[hits] <- "Sulfur Oxidation"
hits <- grep("sulfate adenylyltransferase", counts$Product)
pathway[hits] <- "Sulfite Oxidation"
hits <- grep("phosphosulfate reductase", counts$Product)
pathway[hits] <- "Sulfite Oxidation"
hits <- grep("urea", counts$Product)
pathway[hits] <- "Urea Degradation"
hits <- grep("Urea", counts$Product)
pathway[hits] <- "Urea Degradation"
hits <- grep("chrome c", counts$Product)
pathway[hits] <- "Cytochrome C"
hits <- grep("chrome bd", counts$Product)
pathway[hits] <- "Cytochrome bd"

counts$Pathway <- pathway

# Sum RPKM by pathway and lake

pathway_rpkm <- aggregate(RPKM ~ Lake + Pathway, data = counts, mean)

ggplot(pathway_rpkm, aes(x = Pathway, y = RPKM, fill = Lake)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()

# Cool, but what I'm really interested in is the differences between lakes. Can I pick a lake and always normalize to one?

pathway_rpkm2 <- pathway_rpkm[c(1:6, 18:32, 40:72),]
bars <- unique(pathway_rpkm2$Pathway)
for(i in 1:length(bars)){
  me_amount <- pathway_rpkm2$RPKM[which(pathway_rpkm2$Lake == "ME_epi" & pathway_rpkm2$Pathway == bars[i])]
  ratio <- 1/me_amount
  pathway_rpkm2$RPKM[which(pathway_rpkm2$Pathway == bars[i])] <- pathway_rpkm2$RPKM[which(pathway_rpkm2$Pathway == bars[i])] * ratio
}

x <- ggplot(pathway_rpkm2[which(pathway_rpkm2$Pathway != "Sulfide Reduction"),], aes(x = Pathway, y = RPKM, fill = Lake)) + geom_bar(stat = "identity", position = "dodge") + labs(y = "Normalized RPKM") + scale_y_log10() + coord_flip()

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/blast_results.pdf", x, base_height = 8, base_aspect_ratio = 1)


# Aggregate by product and reshape counts back to wide to input into Lefse
product_rpkm <- aggregate(RPKM ~ Var2 + Product, data = counts, mean)
pathway_rpkm <- aggregate(RPKM ~ Var2 + Pathway, data = counts, mean)
lefse_input <- dcast(pathway_rpkm, Pathway ~ Var2)

# Replace NAs with 0

for (i in 1:dim(lefse_input)[1]){
  row <- lefse_input[i,]
  row[which(is.na(row) == T)] <- 0
  lefse_input[i,] <- row
}

# Add a 1st colum of the lake for each metagenome
# Make products a row. since output is tabular, replace spaces with underscores
rownames(lefse_input) <- lefse_input[,1]
rownames(lefse_input) <- gsub(" ", "_", rownames(lefse_input))
lefse_input <- lefse_input[,2:165]
lakerow <- as.character(lakekey$site[match(colnames(lefse_input), as.character(lakekey$sample))])
lefse_input <- rbind(lakerow, lefse_input)

# Output
write.table(lefse_input, file = "C:/Users/Alex/Desktop/MAGstravaganza/Data_files/lefse_input.txt", quote = F, row.names = T, col.names = T, sep = "\t")
