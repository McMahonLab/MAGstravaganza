############
#MAGstravaganza_Manuscript_plots_2018-01-23

### Make all Manuscript_plots used in the manuscript

#Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
library(ape)
library(zoo)
dadjoke::groan()
############
### Fig 1. 16S vs metagenomic read classifications



############
### Fig 2. Metabolism chart

#Load data
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv")

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

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/ME.heatmap.pdf", ME.heatmap, base_height = 7.5, base_aspect_ratio = 0.75)

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
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/TE.heatmap.pdf", TE.heatmap, base_height = 7.5, base_aspect_ratio = 0.67)

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
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/TH.heatmap.pdf", TH.heatmap, base_height = 7.5, base_aspect_ratio = .73)

################
# Fig 3. Glycoside hydrolases

# Load data
mag_data <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)
mag_data$Lake <- as.character(mag_data$Lake)

# Set up data frame to store results
gh_df <- mag_data[,c(1,3,4,10)]
gh_density <- c()
gh_diversity <- c()

# Read in each MAG's CAZy results and store GH information
for(i in 1:dim(mag_data)[1]){
  datafile <- read.table(paste("C:/Users/Alex/Desktop/MAGstravaganza/hydrolases/", mag_data$IMG_OID[i], ".txt", sep = ""))
  hits <- sub(".hmm", "", datafile$V2)
  types <- substr(hits, start = 1, stop = 2)
  gh_hits <- hits[which(types == "GH")]
  gh_density[i] <- length(gh_hits)/mag_data$Gene_Count[i]*100
  gh_diversity[i] <- length(unique(gh_hits))
}

gh_df$Diversity <- gh_diversity
gh_df$Density <- gh_density

gh_df$Order <- sapply(strsplit(mag_data$Taxonomy, ";"), "[", 3)
gh_df$Order <- factor(gh_df$Order, levels = rev(c("Actinomycetales", "Chlamydiales", "Cytophagales", "Mycoplasmatales", "Planctomycetales", "Puniceicoccales", "Rhodocyclales", "Sphingomonadales", "Xanthomonadales", "Solibacterales", "Bdellovibrionales", "Chlorobiales", "Holophagales", "Rickettsiales", "Bacteroidales", "Campylobacterales", "Desulfobacterales", "Desulfuromonadales", "Gallionellales", "Ignavibacteriales", "Legionellales", "Nitrosomonadales", "Pseudomonadales", "Rhizobiales", "Solirubrobacterales", "Flavobacteriales", "Acidimicrobiales", "Burkholderiales", "Methylococcales", "Methylophilales", "Sphingobacteriales", "Verrucomicrobiales")))

# Plot the heatmap
gh.agg <- aggregate(Density ~ Lake + Order, data = gh_df, mean)
gh_density_plot <- ggplot(data = gh.agg[which(gh.agg$Order != "[Blank]" & is.na(gh.agg$Order) == F), ], aes(x = Lake, y = Order, fill = Density)) + geom_tile(color = "black") + scale_fill_gradient2(low = "white", mid = "#8c96c6", high = "#810f7c", midpoint = 4) + labs(x = "", y = "") + background_grid(major = "xy") + theme(axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/GH.density.heatmap.pdf", gh_density_plot, base_height = 5, base_aspect_ratio = 1)

# Correlation of density and diversity
cor.test(gh_df$Gene_Count, gh_df$Density)

#################
# Fig 4 - Nitrogen and Sulfur Metabolisms
# Panels A-C
# Use data from Excel pathway analysis

table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T, row.names = 2)

colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 2:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}

# Read in metadata info to find out which MAGs came from which lake

metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

# Mendota Nitrogen

MEnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Polyamine Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
# Convert to presence/absence and sum by pathway
MEnitrogen_totals <- rowSums(MEnitrogen > 0)
# Select pathways to plot
MEnitrogen_totals <- MEnitrogen_totals[c(1, 2, 4, 5, 7, 11)]
MEnitrogen_df <- data.frame(names(MEnitrogen_totals), MEnitrogen_totals)
colnames(MEnitrogen_df) <- c("Pathway", "MAGs")

#Trout Bog Nitrogen
TEnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Polyamine Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TEnitrogen_totals <- rowSums(TEnitrogen > 0)
TEnitrogen_totals <- TEnitrogen_totals[c(1, 2, 4, 5, 7, 11)]
TEnitrogen_df <- data.frame(names(TEnitrogen_totals), TEnitrogen_totals)
colnames(TEnitrogen_df) <- c("Pathway", "MAGs")

THnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Polyamine Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
THnitrogen_totals <- rowSums(THnitrogen > 0)
THnitrogen_totals <- THnitrogen_totals[c(1, 2, 4, 5, 7, 11)]
THnitrogen_df <- data.frame(names(THnitrogen_totals), THnitrogen_totals)
colnames(THnitrogen_df) <- c("Pathway", "MAGs")

#Repeat with sulfur
MEsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
MEsulfur_totals <- rowSums(MEsulfur > 0)
MEsulfur_totals <- MEsulfur_totals[c(1,3,4)]
MEsulfur_df <- data.frame(names(MEsulfur_totals), MEsulfur_totals)
colnames(MEsulfur_df) <- c("Pathway", "MAGs")

#Trout Bog totals
TEsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TEsulfur_totals <- rowSums(TEsulfur > 0)
TEsulfur_totals <- TEsulfur_totals[c(1,3,4)]
TEsulfur_df <- data.frame(names(TEsulfur_totals), TEsulfur_totals)
colnames(TEsulfur_df) <- c("Pathway", "MAGs")

THsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
THsulfur_totals <- rowSums(THsulfur > 0)
THsulfur_totals <- THsulfur_totals[c(1,3,4)]
THsulfur_df <- data.frame(names(THsulfur_totals), THsulfur_totals)
colnames(THsulfur_df) <- c("Pathway", "MAGs")

# Combine sulfur and nitrogen tables into one dataframe per lake

MEdata <- rbind(MEnitrogen_df, MEsulfur_df)
TEdata <- rbind(TEnitrogen_df, TEsulfur_df)
THdata <- rbind(THnitrogen_df, THsulfur_df)

ME_plot <- ggplot(data = MEdata, aes(x = Pathway, y = MAGs/102*100)) + geom_bar(stat = "identity", fill = "#FF0000", color = "black") + labs(x = "", y = "") + scale_y_continuous(limits = c(0, 95)) + coord_flip() + labs(y = "% MAGs with pathway", title = "Lake Mendota")
TE_plot <- ggplot(data = TEdata, aes(x = Pathway, y = MAGs/31*100)) + geom_bar(stat = "identity", fill = "#8EEB00", color = "black") + labs(x = "", y = "") + scale_y_continuous(limits = c(0, 95)) + coord_flip() + labs(y = "% MAGs with pathway", title = "Trout Bog Epiliminion")
TH_plot <- ggplot(data = THdata, aes(x = Pathway, y = MAGs/63*100)) + geom_bar(stat = "identity", fill = "#00A287", color = "black") + labs(x = "", y = "") + scale_y_continuous(limits = c(0, 95)) + coord_flip()  + labs(y = "% MAGs with pathway", title = "Trout Bog Hypoliminion")

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/fig4A_mendota.pdf", ME_plot, base_height = 3, base_aspect_ratio = 1.7)
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/fig4B_TBE.pdf", TE_plot, base_height = 3, base_aspect_ratio = 1.7)
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/fig4C_TBH.pdf", TH_plot, base_height = 3, base_aspect_ratio = 1.7)

# Panel D (or B) - the tree
# Load data

treefile <- read.tree("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/allMAGS.nwk")
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T, colClasses = c("character", "character", "factor", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T)

# These genomes don't meet the stricter quality cutoffs I'm using in the pathway analysis (> 50% complete)
treefile <- drop.tip(treefile, "2582580570.genes.fna")
treefile <- drop.tip(treefile, "2582580574.genes.fna")
treefile <- drop.tip(treefile, "2593339180.genes.fna")
treefile <- drop.tip(treefile, "2593339186.genes.fna")
treefile <- drop.tip(treefile, "2593339193.genes.fna")
treefile <- drop.tip(treefile, "2588253500.genes.fna")
treefile <- drop.tip(treefile, "2582580534.genes.fna")
treefile <- drop.tip(treefile, "2582580704.genes.fna")
treefile <- drop.tip(treefile, "2582580707.genes.fna")
treefile <- drop.tip(treefile, "2582580709.genes.fna")
treefile <- drop.tip(treefile, "2582580714.genes.fna")
treefile <- drop.tip(treefile, "2582580604.genes.fna")

# Fix column names and grab N2 fixation data
colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 3:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}
hasN2 <- table[3, ] > 0

tips <- treefile$tip.label
tips <- gsub(".genes.fna", "", tips)

# Make keys for coloring tips
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

# Plot with colors by lake
key <- data.frame(tips, lakes, Nfix)
lake.colors2 <- c()
lake.colors2[which(key$lakes == "Mendota")] <- "#FF0000"
lake.colors2[which(key$lakes == "Trout Bog Epilimnion")] <- "#8EEB00"
lake.colors2[which(key$lakes == "Trout Bog Hypolimnion")] <- "#00A287"

plot(treefile, type = "fan", show.tip.label = F, main = "", edge.width = 3)
tiplabels(pch = 21, col = "black", bg = lake.colors2, cex = 2)

# And finally, the base tree for outputting.
pdf("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/tree_base.pdf", height = 10, width = 10)
plot(treefile, type = "fan", show.tip.label = F, main = "", edge.width = 3)
dev.off()

# Output the tree as a pdf file, import into Inkscape, and manually do the following:
# - color branches by phylum
# - add different shaped tips for each lake
# - color tips that are nitrogen fixing MAGs
# It's a pain but I can't figure out how to get the tree exactly the way I want it with code!

##################
# Figure 5
# Trace methylotrophs over time

Mendota_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/Mendota_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
Mendota_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
Mendota_dates$sample <- as.character(Mendota_dates$sample)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

TE_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/TE_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TE_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
TE_dates$sample <- as.character(TE_dates$sample)

TH_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/TH_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TH_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")
TH_dates$sample <- as.character(TH_dates$sample)

# Convert to RPKM (reads per kilobase per million reads)
# general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
size_vector <- metadata$Genome_Size[match(Mendota_raw_input$V1, metadata$IMG_OID)]
Mendota_RPKM <- Mendota_raw_input$V3/(size_vector/1000 * Mendota_raw_input$V4/1000000)
ME_results <- data.frame(Mendota_raw_input[,1:2], Mendota_RPKM)
colnames(ME_results) <- c("MAG", "metaG", "RPKM")
ME_results$Date <- as.Date(Mendota_dates$date[match(ME_results$metaG, Mendota_dates$sample)], format = "%m/%d/%y")
ME_results$Julian_Date <- as.numeric(format(ME_results$Date, "%j"))
ME_results <- ME_results[which(is.na(ME_results$Date) == F),]

size_vector <- metadata$Genome_Size[match(TE_raw_input$V1, metadata$IMG_OID)]
TE_RPKM <- TE_raw_input$V3/(size_vector/1000 * TE_raw_input$V4/1000000)
TE_results <- data.frame(TE_raw_input[,1:2], TE_RPKM)
colnames(TE_results) <- c("MAG", "metaG", "RPKM")
TE_results$Date <- as.Date(TE_dates$date[match(TE_results$metaG, TE_dates$sample)], format = "%m/%d/%y")
TE_results$Julian_Date <- as.numeric(format(TE_results$Date, "%j"))
TE_results <- TE_results[which(is.na(TE_results$Date) == F),]

size_vector <- metadata$Genome_Size[match(TH_raw_input$V1, metadata$IMG_OID)]
TH_RPKM <- TH_raw_input$V3/(size_vector/1000 * TH_raw_input$V4/1000000)
TH_results <- data.frame(TH_raw_input[,1:2], TH_RPKM)
colnames(TH_results) <- c("MAG", "metaG", "RPKM")
TH_results$Date <- as.Date(TH_dates$date[match(TH_results$metaG, TH_dates$sample)], format = "%m/%d/%y")
TH_results$Julian_Date <- as.numeric(format(TH_results$Date, "%j"))
TH_results <- TH_results[which(is.na(TH_results$Date) == F),]

# Get taxonomy
metadata$phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)
metadata$class<- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 2)
metadata$order <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 3)

ME_results$phylum <- metadata$phylum[match(ME_results$MAG, metadata$IMG_OID)]
ME_results$class <- metadata$class[match(ME_results$MAG, metadata$IMG_OID)]
ME_results$order <- metadata$order[match(ME_results$MAG, metadata$IMG_OID)]

TE_results$phylum <- metadata$phylum[match(TE_results$MAG, metadata$IMG_OID)]
TE_results$class <- metadata$class[match(TE_results$MAG, metadata$IMG_OID)]
TE_results$order <- metadata$order[match(TE_results$MAG, metadata$IMG_OID)]

TH_results$phylum <- metadata$phylum[match(TH_results$MAG, metadata$IMG_OID)]
TH_results$class <- metadata$class[match(TH_results$MAG, metadata$IMG_OID)]
TH_results$order <- metadata$order[match(TH_results$MAG, metadata$IMG_OID)]


seasons <- data.frame(x1 = c(75, 152, 244), x2 = c(151, 243, 334), y1 = c(0, 0, 0), y2 = c(Inf, Inf, Inf))
trace <- data.frame(input = ME_results$RPKM[which(ME_results$order == "Methylococcales")], dates = ME_results$Julian_Date[which(ME_results$order == "Methylococcales")])
trace <- trace[order(trace$dates),]
trace$line <- c(trace$input[1], mean(trace$input[1:2]), mean(trace$input[1:3]), mean(trace$input[1:4]), mean(trace$input[1:5]), mean(trace$input[1:6]), mean(trace$input[1:7]), mean(trace$input[1:8]), rollmean(trace$input, 9))
trace$MAG <- rep("2582580566", length(trace$line))

p1 <- ggplot(data = ME_results[which(ME_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = ME_results$RPKM[which(ME_results$order == "Methylophilales")], dates = ME_results$Julian_Date[which(ME_results$order == "Methylophilales")], MAG =  ME_results$MAG[which(ME_results$order == "Methylophilales")])
trace1 <- trace[which(trace$MAG == "2582580561"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580598"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace <- rbind(trace1, trace2)

p2 <- ggplot(data = ME_results[which(ME_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TE_results$RPKM[which(TE_results$order == "Methylococcales")], dates = TE_results$Julian_Date[which(TE_results$order == "Methylococcales")], MAG =  TE_results$MAG[which(TE_results$order == "Methylococcales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580615"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580639"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2582580646"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace <- rbind(trace1, trace2, trace3)

p3 <- ggplot(data = TE_results[which(TE_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TE_results$RPKM[which(TE_results$order == "Methylophilales")], dates = TE_results$Julian_Date[which(TE_results$order == "Methylophilales")], MAG =  TE_results$MAG[which(TE_results$order == "Methylophilales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580623"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580649"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace <- rbind(trace1, trace2)

p4 <- ggplot(data = TE_results[which(TE_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TH_results$RPKM[which(TH_results$order == "Methylococcales")], dates = TH_results$Julian_Date[which(TH_results$order == "Methylococcales")], MAG =  TH_results$MAG[which(TH_results$order == "Methylococcales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580677"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2593339178"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2593339191"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace <- rbind(trace1, trace2, trace3)

p5 <- ggplot(data = TH_results[which(TH_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TH_results$RPKM[which(TH_results$order == "Methylophilales")], dates = TH_results$Julian_Date[which(TH_results$order == "Methylophilales")], MAG =  TH_results$MAG[which(TH_results$order == "Methylophilales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580712"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2593339176"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2593339182"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace4 <- trace[which(trace$MAG == "2593339192"), ]
trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], mean(trace4$input[1:2]), mean(trace4$input[1:3]), mean(trace4$input[1:4]), mean(trace4$input[1:5]), mean(trace4$input[1:6]), mean(trace4$input[1:7]), mean(trace4$input[1:8]), rollmean(trace4$input, 9))
trace <- rbind(trace1, trace2, trace3, trace4)

p6 <- ggplot(data = TH_results[which(TH_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

meth <- plot_grid(p1, p3, p5, p2, p4, p6, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/Fig5.pdf", meth, base_height = 5, base_aspect_ratio = 3)
