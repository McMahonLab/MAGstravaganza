############
#MAGstravaganza_Manuscript_plots_2018-01-23

### Make all Manuscript_plots used in the manuscript

#Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
library(ape)
library(zoo)
############
### Fig 1. 16S vs metagenomic read classifications

# Read in tag data. Note that TB and ME were amplified with different primers (ugh)

mendota_tags <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/mendota_taxonomy_phylum_bacteria_archaea_L2.csv", header = T, row.names = 1)
trout_tags <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/Trout_Bog_220_10_pyroclust.phylum_classification_percent.csv", header = T, row.names = 1, fill = TRUE)
trout_tags[is.na(trout_tags)] <- 0
#Split Trout Bog data into epi and hypo
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

tbe_tags <- trout_tags[,which(substrRight(colnames(trout_tags), 1) == "E")]
tbh_tags <- trout_tags[,which(substrRight(colnames(trout_tags), 1) == "H")]

# Get average % of community for each phylum, and switch to long format

mendota_tag_totals <- rowSums(mendota_tags)/dim(mendota_tags)[2]
tbe_tag_totals <- rowSums(tbe_tags)/dim(tbe_tags)[2]
tbh_tag_totals <- rowSums(tbh_tags)/dim(tbh_tags)[2]

phylum_names <- c(rownames(mendota_tags), rownames(tbe_tags), rownames(tbh_tags))
phylum_percents <- c(mendota_tag_totals, tbe_tag_totals, tbh_tag_totals)
phylum_lake <- c(rep("Mendota", length(mendota_tag_totals)), rep("Trout Bog Epilimnion", length(tbe_tag_totals)), rep("Trout Bog Hypolimnion", length(tbh_tag_totals)))

tag_data <- data.frame(phylum_names, phylum_percents, phylum_lake)
tag_data <- tag_data[which(tag_data$phylum_percents >= 1), ]

# Keep phyla and match order to the MAGs plot, then add a few of the most abundant or so
# Which are most abundant?
head(tag_data[order(tag_data$phylum_percents, decreasing = T), ], 50)
tag_data <- tag_data[which(tag_data$phylum_names == "Acidobacteria" | tag_data$phylum_names == "Actinobacteria" | tag_data$phylum_names == "Bacteroidetes" | tag_data$phylum_names == "Chlamydiae" | tag_data$phylum_names == "Chlorobi" | tag_data$phylum_names == "Cyanobacteria" | tag_data$phylum_names == "Ignavibacteria" | tag_data$phylum_names == "Planctomycetes" | tag_data$phylum_names == "Proteobacteria" | tag_data$phylum_names == "Tenericutes" | tag_data$phylum_names == "Verrucomicrobia" | tag_data$phylum_names == "Cryptophyta" | tag_data$phylum_names == "stramenopiles" | tag_data$phylum_names == "Firmicutes"),]
# Order the same as the MAGs plot and change names or add zeroes as necessary
tag_data$phylum_names <- as.character(tag_data$phylum_names)
tag_data$phylum_names[which(tag_data$phylum_names == "stramenopiles")] <- "Heterokonta"

zeroes <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
zphyla <- c("Chlamydiae", "Chlamydiae", "Chlamydiae", "Tenericutes", "Tenericutes", "Tenericutes", "Ignavibacteria", "Ignavibacteria", "Ignavibacteria")
zlakes <- c("Mendota", "Trout Bog Epilimnion", "Trout Bog Hypolimnion", "Mendota", "Trout Bog Epilimnion", "Trout Bog Hypolimnion", "Mendota", "Trout Bog Epilimnion", "Trout Bog Hypolimnion")
zdata <- data.frame(zphyla, zeroes, zlakes)
colnames(zdata) <- colnames(tag_data)
tag_data <- rbind(tag_data, zdata)

tag_data$phylum_names <- factor(tag_data$phylum_names, levels = c("Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chlamydiae", "Chlorobi", "Cyanobacteria", "Ignavibacteria", "Planctomycetes", "Proteobacteria", "Tenericutes", "Verrucomicrobia", "Cryptophyta", "Heterokonta", "Firmicutes"))

panelA <- ggplot(data = tag_data, aes(x = phylum_names, y = phylum_percents, fill = phylum_lake)) + geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Mean Percent of Reads", x = NULL, title = "Phylum Assignments of 16S Ribosomal Reads") + scale_fill_manual(values = c("#FF0000", "#8EEB00", "#00A287"))

# Panel B

# Get # of reads assigned to each phylum (mapped MAG) by lake-layer
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

# Normalize (reads per kilobase MAG per million reads)
# general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
# but technically we're not calling it RPKM because that number is meant for genomes, not metagenomes
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

ME_results$phylum <- metadata$phylum[match(ME_results$MAG, metadata$IMG_OID)]
TE_results$phylum <- metadata$phylum[match(TE_results$MAG, metadata$IMG_OID)]
TH_results$phylum <- metadata$phylum[match(TH_results$MAG, metadata$IMG_OID)]

ME_results$Taxonomy <- metadata$Taxonomy[match(ME_results$MAG, metadata$IMG_OID)]
TE_results$Taxonomy <- metadata$Taxonomy[match(TE_results$MAG, metadata$IMG_OID)]
TH_results$Taxonomy <- metadata$Taxonomy[match(TH_results$MAG, metadata$IMG_OID)]

# Aggregate by phylum

agg_ME_phyla <- aggregate(RPKM ~ phylum, data = ME_results, mean)
agg_TE_phyla <- aggregate(RPKM ~ phylum, data = TE_results, mean)
agg_TH_phyla <- aggregate(RPKM ~ phylum, data = TH_results, mean)

agg_phyla <- rbind(agg_ME_phyla, agg_TE_phyla, agg_TH_phyla)
agg_phyla$Lake <- c(rep("Mendota", dim(agg_ME_phyla)[1]), rep("Trout Bog Epilimnion", dim(agg_TE_phyla)[1]), rep("Trout Bog Hypolimnion", dim(agg_TH_phyla)[1]))

# Remove unclassified
agg_phyla <- agg_phyla[grep("Blank", agg_phyla$phylum, invert = T), ]

panelB <- ggplot(data = agg_phyla, aes(x = phylum, y = RPKM, fill = Lake)) + geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "Mean Normalized Reads", x = NULL, title = "Metagenomic reads recruited by time series MAGs") + scale_fill_manual(values = c("#FF0000", "#8EEB00", "#00A287"))

panelA <- panelA + theme(legend.position = "none")
fig1 <- plot_grid(panelA, panelB, ncol = 1, labels = c("A", "B"))
save_plot("C:/Users/Alex//Desktop/MAGstravaganza/Manuscript_plots/fig1.pdf", fig1, base_height = 8, base_aspect_ratio = 1.6)


################
# Fig 2. Glycoside hydrolases

# Load data
mag_data <- read.csv("C:/Users/Alex//Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)
mag_data$Lake <- as.character(mag_data$Lake)

# Set up data frame to store results
gh_df <- mag_data[,c(1,3,4,10)]
gh_density <- c()
gh_diversity <- c()
tophit <- c()
tophit_counts <- c()

# Read in each MAG's CAZy results and store GH information
for(i in 1:dim(mag_data)[1]){
  datafile <- read.table(paste("C:/Users/Alex/Desktop/MAGstravaganza/dbCAN_results/", mag_data$IMG_OID[i], ".txt", sep = ""))
  hits <- sub(".hmm", "", datafile$V2)
  types <- substr(hits, start = 1, stop = 2)
  gh_hits <- hits[which(types == "GH")]
  gh_density[i] <- length(gh_hits)/mag_data$Gene_Count[i]*100
  gh_diversity[i] <- length(unique(gh_hits))
  gh_counts <- table(gh_hits)
  tophit[i] <- names(gh_counts)[which(gh_counts == max(gh_counts))]
  tophit_counts[i] <- as.numeric(gh_counts[which(gh_counts == max(gh_counts))])
}

gh_df$Diversity <- gh_diversity
gh_df$Density <- gh_density
gh_df$TopGHHit <- tophit
gh_df$TopGHHit_Count <- tophit_counts

gh_df$Order <- sapply(strsplit(mag_data$Taxonomy, ";"), "[", 3)
gh_df$Order <- factor(gh_df$Order, levels = rev(c("Actinomycetales", "Chlamydiales", "Cytophagales", "Mycoplasmatales", "Planctomycetales", "Puniceicoccales", "Rhodocyclales", "Sphingomonadales", "Xanthomonadales", "Solibacterales", "Bdellovibrionales", "Chlorobiales", "Holophagales", "Rickettsiales", "Bacteroidales", "Campylobacterales", "Desulfobacterales", "Desulfuromonadales", "Gallionellales", "Ignavibacteriales", "Legionellales", "Nitrosomonadales", "Pseudomonadales", "Rhizobiales", "Solirubrobacterales", "Flavobacteriales", "Acidimicrobiales", "Burkholderiales", "Methylococcales", "Methylophilales", "Sphingobacteriales", "Verrucomicrobiales")))

# Output data at this stage for supplemental file
write.csv(gh_df, file = "C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/dbCANN_results.csv", quote = F, row.names = F)

# Plot the heatmap
gh.agg <- aggregate(Density ~ Lake + Order, data = gh_df, mean)
gh_density_plot <- ggplot(data = gh.agg[which(gh.agg$Order != "[Blank]" & is.na(gh.agg$Order) == F), ], aes(x = Lake, y = Order, fill = Density)) + geom_tile(color = "black") + scale_fill_gradient2(low = "white", mid = "#8c96c6", high = "#810f7c", midpoint = 4) + labs(x = "", y = "") + background_grid(major = "xy") + theme(axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/GH.density.heatmap.pdf", gh_density_plot, base_height = 5, base_aspect_ratio = 1)

# Correlation of density and diversity
cor.test(gh_df$Gene_Count, gh_df$Density)

############
### Fig 3. Metabolism chart

#Load data
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv")
# Fake taxonomy data so that Proteobacterial classes function as phyla
metadata$Taxonomy <- gsub("Proteobacteria;", "", metadata$Taxonomy)

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
table <- table[which(table$Pathway != "Nitrogen fixation, Fe-only"), ]
table$Pathway[which(table$Pathway == "Nitrogen fixation, Mo-Fe")] <- "Nitrogen fixation"
ME.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
ME.MAGs$Pathways <- make.unique(table$Pathway)

ME.long <- melt(ME.MAGs)
ME.long$Pathways <- factor(ME.long$Pathways, levels = table$Pathway)
ME.long$variable <- factor(ME.long$variable, levels = colnames(table))
ME.long$Taxonomy <- phylum[match(ME.long$variable, metadata$IMG_OID)]
ME.long$Taxonomy[which(ME.long$Taxonomy == "[Blank]")] <- "Unclassified"
ME.long <- ME.long[which(is.na(ME.long$Pathways) == F), ]
ME.agg <- aggregate(value ~ Pathways + Taxonomy, ME.long, mean)
# Don't plot unclassified groups
ME.agg <- ME.agg[which(ME.agg$Taxonomy != "Unclassified"),]

ME.heatmap <- ggplot(data = ME.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#FF0000", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/ME.heatmap.pdf", ME.heatmap, base_height = 7.75, base_aspect_ratio = 0.73)

TE.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TE.MAGs$Pathways <- make.unique(table$Pathway)

TE.long <- melt(TE.MAGs)
TE.long$Pathways <- factor(TE.long$Pathways, levels = table$Pathway)
TE.long$variable <- factor(TE.long$variable, levels = colnames(table))
TE.long$Taxonomy <- phylum[match(TE.long$variable, metadata$IMG_OID)]
TE.long$Taxonomy[which(TE.long$Taxonomy == "[Blank]")] <- "Unclassified"
TE.long <- TE.long[which(is.na(TE.long$Pathways) == F), ]
TE.agg <- aggregate(value ~ Pathways + Taxonomy, TE.long, mean)

TE.agg$value[which(TE.agg$Pathways == "Calvin Cycle" & TE.agg$Taxonomy == "Chlorobi")] <- 0
TE.agg$value[which(TE.agg$Pathways == "Dissimilatory sulfate reduction" & TE.agg$Taxonomy == "Chlorobi")] <- 0
TE.agg <- TE.agg[which(TE.agg$Taxonomy != "Unclassified"),]

TE.heatmap <- ggplot(data = TE.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#8EEB00", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/TE.heatmap.pdf", TE.heatmap, base_height = 7.75, base_aspect_ratio = 0.69)

TH.MAGs <- table[,which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
TH.MAGs$Pathways <- make.unique(table$Pathway)

TH.long <- melt(TH.MAGs)
TH.long$Pathways <- factor(TH.long$Pathways, levels = table$Pathway)
TH.long$variable <- factor(TH.long$variable, levels = colnames(table))
TH.long$Taxonomy <- phylum[match(TH.long$variable, metadata$IMG_OID)]
TH.long$Taxonomy[which(TH.long$Taxonomy == "[Blank]")] <- "Unclassified"
TH.long <- TH.long[which(is.na(TH.long$Pathways) == F), ]
TH.agg <- aggregate(value ~ Pathways + Taxonomy, TH.long, mean)

TH.agg$value[which(TH.agg$Pathways == "Dissimilatory sulfate reduction" & TH.agg$Taxonomy == "Chlorobi")] <- 0
TH.agg$value[which(TH.agg$Pathways == "Calvin Cycle" & TH.agg$Taxonomy == "Chlorobi")] <- 0

TH.agg <- TH.agg[which(TH.agg$Taxonomy != "Unclassified"),]

TH.heatmap <- ggplot(data = TH.agg, aes(y = Pathways, x = Taxonomy, fill = value)) + geom_tile(color = "black") + labs(x = "", y = "") + scale_fill_gradient2(low = "white", high = "#00A287", midpoint = 25) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/TH.heatmap.pdf", TH.heatmap, base_height = 7.75, base_aspect_ratio = .73)

#################
# Fig 4 - Nitrogen and Sulfur Metabolisms
# Panel A
# Use marker gene data

# Read in datasets - marker gene counts + LEfSe significance

table <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/marker_gene_table.txt", header = T, row.names = 1)
lakekey <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metagenome_metadata.csv", header = T)
lakerow <- as.character(lakekey$site[match(colnames(table), as.character(lakekey$sample))])

epilefse <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/marker_genes_LDA_significance.csv", header = T)
TBlefse <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/TBmarker_genes_LDA_significance.csv", header = T)

genedata <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metabolic_gene_info.txt", fill = TRUE)
# for my selected marker genes, get the number of hits of each one by lake-layer

genes2keep <- genedata$V1[c(grep("rubisco", genedata$V3), grep("citrate_lyase", genedata$V3), grep("nitrogenase", genedata$V3), grep("Sox", genedata$V3), grep("dissimilatory_sulfite_reductase", genedata$V3), grep("urease", genedata$V3), grep("^sulfite_reductase", genedata$V3), grep("sulfate_adenylyltransferase", genedata$V3), grep("nitrate_reductase", genedata$V3), grep("nitrite_reductase", genedata$V3), grep("nitric_oxide_reductase", genedata$V3), grep("nitrous_oxide_reductase", genedata$V3), grep("sulfide_quinone_reductase", genedata$V3))]

table2keep <- table[match(genes2keep, rownames(table)),]
table2keep$genes <- rownames(table2keep)
table2keep <- melt(table2keep)
table2keep$annotations <- genedata$V3[match(table2keep$genes, genedata$V1)]
table2keep$lake <- as.character(lakekey$site[match(table2keep$variable, as.character(lakekey$sample))])
table2keep <- table2keep[which(is.na(table2keep$lake) == F),]

#Sum by function and lake, then divide by number of metagenomes from each lake
#First combine categories that are variations of the same thing
table2keep$annotations <- gsub("thiosulfohydrolase_SoxB", "SOX", table2keep$annotations)
table2keep$annotations <- gsub("thiosulfate_oxidation_carrier_SoxY", "SOX", table2keep$annotations)
table2keep$annotations <- gsub("sulfite_dehydrogenase_SoxC", "SOX", table2keep$annotations)
table2keep$annotations <- gsub("dissimilatory_sulfite_reductase", "sulfite_reductase", table2keep$annotations)
table2keep$annotations <- gsub("periplasmic_nitrate_reductase", "nitrate_reductase", table2keep$annotations)
table2keep$annotations <- gsub("nitrate_reductase_cytochrome_c-type", "nitrate_reductase", table2keep$annotations)
table2keep$annotations <- gsub("formate_dependent_cytochrome_c_nitrite_reductase", "nitrite_reductase", table2keep$annotations)
table2keep$annotations <- gsub("cytochrome_c_nitrite_reductase", "nitrite_reductase", table2keep$annotations)
table2keep$annotations <- gsub("nitrogenase_VaFe", "nitrogenase", table2keep$annotations)
table2keep$annotations <- gsub("nitrogenase_MoFe", "nitrogenase", table2keep$annotations)
table2keep$annotations <- gsub("Fe-only_nitrogenase", "nitrogenase", table2keep$annotations)
table2keep$annotations <- gsub("nitrogenase_iron_protein", "nitrogenase", table2keep$annotations)
table2keep$annotations <- gsub("_", " ", table2keep$annotations)
table2keep$annotations <- factor(table2keep$annotations, levels = c("rubisco", "citrate lyase", "methane ammonia monooxygenase", "methanol dehydrogenase", "urease", "nitrogenase", "nitrate reductase", "nitrite reductase", "nitric oxide reductase", "nitrous oxide reductase", "SOX", "sulfate adenylyltransferase", "sulfide quinone reductase", "sulfite reductase"))

# marker_sums <- aggregate(value ~ lake + annotations, table2keep, sum)
# marker_sums$value[which(marker_sums$lake == "ME_epi")] <- marker_sums$value[which(marker_sums$lake == "ME_epi")]/length(which(lakekey$site == "ME_epi"))
# marker_sums$value[which(marker_sums$lake == "TB_epi")] <- marker_sums$value[which(marker_sums$lake == "TB_epi")]/length(which(lakekey$site == "TB_epi"))
# marker_sums$value[which(marker_sums$lake == "TB_hypo")] <- marker_sums$value[which(marker_sums$lake == "TB_hypo")]/length(which(lakekey$site == "TB_hypo"))
# 
# # Remove categories that are too low in abundance to see on the plot
# marker_sums <- marker_sums[which(marker_sums$annotations != "soluble_methane_ammonia_monooxygenase"), ]
# marker_sums <- marker_sums[which(marker_sums$annotations != "methylamine_dehydrogenase"), ]
# 
# # Remove underscore separators and use factors to order properly
# marker_sums$annotations <- gsub("_", " ", marker_sums$annotations)
# marker_sums$annotations <- factor(marker_sums$annotations, levels = c("rubisco", "citrate lyase", "methane ammonia monooxygenase", "methanol dehydrogenase", "urease", "nitrogenase", "nitrate reductase", "nitrite reductase", "nitric oxide reductase", "nitrous oxide reductase", "SOX", "sulfate adenylyltransferase", "sulfide quinone reductase", "sulfite reductase"))

plot.boxes <- aggregate(value ~ variable + annotations + lake, table2keep, sum)
plot.boxes1 <- plot.boxes[which(plot.boxes$annotations == "sulfide quinone reductase" | plot.boxes$annotations == "sulfate adenylyltransferase" | plot.boxes$annotations == "SOX"), ]
plot.boxes2 <- plot.boxes[which(plot.boxes$annotations != "sulfide quinone reductase" & plot.boxes$annotations != "sulfate adenylyltransferase" & plot.boxes$annotations != "SOX"), ]

p1 <- ggplot(data = plot.boxes1, aes(x = annotations, y = value, fill = lake)) + geom_boxplot() + scale_fill_manual(values = c("#FF0000", "#8EEB00", "#00A287")) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = NULL, y = "Counts per Metagenome") + coord_flip()
p2 <- ggplot(data = plot.boxes2, aes(x = annotations, y = value, fill = lake)) + geom_boxplot() + scale_fill_manual(values = c("#FF0000", "#8EEB00", "#00A287")) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = NULL, y = "Counts per Metagenome") + coord_flip()

fig4 <- plot_grid(p1, p2, nrow = 2, axis = "l", rel_heights = c(1.75, 4), align = "v")

save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/marker_gene_hits.pdf", fig4, base_height = 8, base_aspect_ratio = 1)

# # Use this code to look at significance between lakes
# searchterm <- "rubisco"
# epilefse[grep(searchterm, epilefse$annotation),]
# TBlefse[grep(searchterm, TBlefse$annotation),]
# 
# # Output a wide version for testing averages in lefse
# sample_sums <- aggregate(value ~ variable + annotations, table2keep, sum)
# sample_sums <- sample_sums[which(sample_sums$annotations != "soluble_methane_ammonia_monooxygenase"), ]
# sample_sums <- sample_sums[which(sample_sums$annotations != "methylamine_dehydrogenase"), ]
# wide_sample_sums <- dcast(sample_sums, annotations ~ variable)
# lakerow <- as.character(lakekey$site[match(colnames(wide_sample_sums), as.character(lakekey$sample))])
# wide_sample_sums <- rbind(lakerow, wide_sample_sums)
# write.table(wide_sample_sums, "C:/Users/Alex/Desktop/MAGstravaganza/Data_files/funcational_lefse_input.txt")

# Use Kruskal-Wallis Test to look for significant differences in amounts between lakes
genes <- c("rubisco", "citrate lyase", "urease", "nitrogenase", "nitrate reductase", "nitrite reductase", "nitric oxide reductase", "nitrous oxide reductase", "SOX", "sulfate adenylyltransferase", "sulfide quinone reductase", "sulfite reductase")

for(i in 1:length(genes)){
  gene <- genes[i]
  genetable <- plot.boxes[which(plot.boxes$annotations == gene), ]
  x <- pairwise.wilcox.test(x = genetable$value, g = genetable$lake, p.adj = "bonf", paired = F)
  stat.table <- x$p.value
  item1 <- paste(gene, colnames(stat.table)[1], rownames(stat.table)[1], round(stat.table[1,1], 2))
  item2 <- paste(gene, colnames(stat.table)[2], rownames(stat.table)[2], round(stat.table[2,2], 2))
  print(item1)
  print(item2)
}

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
# Trace Cyanobacteria over time in Mendota

Mendota_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/Mendota_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
Mendota_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
Mendota_dates$sample <- as.character(Mendota_dates$sample)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

# Convert to RPKM (reads per kilobase per million reads) ish - not really RPKM because it's for metagenomes
# general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
size_vector <- metadata$Genome_Size[match(Mendota_raw_input$V1, metadata$IMG_OID)]
Mendota_RPKM <- Mendota_raw_input$V3/(size_vector/1000 * Mendota_raw_input$V4/1000000)
ME_results <- data.frame(Mendota_raw_input[,1:2], Mendota_RPKM)
colnames(ME_results) <- c("MAG", "metaG", "RPKM")
ME_results$Date <- as.Date(Mendota_dates$date[match(ME_results$metaG, Mendota_dates$sample)], format = "%m/%d/%y")
ME_results$Julian_Date <- as.numeric(format(ME_results$Date, "%j"))
ME_results <- ME_results[which(is.na(ME_results$Date) == F),]

# Get taxonomy
metadata$phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)

ME_results$phylum <- metadata$phylum[match(ME_results$MAG, metadata$IMG_OID)]
ME_results$year <- substr(ME_results$Date, start = 1, stop = 4)


seasons <- data.frame(x1 = c(75, 152, 244), x2 = c(151, 243, 334), y1 = c(0, 0, 0), y2 = c(Inf, Inf, Inf))
trace <- data.frame(input = ME_results$RPKM[which(ME_results$phylum == "Cyanobacteria")], dates = ME_results$Julian_Date[which(ME_results$phylum == "Cyanobacteria")], year = ME_results$year[which(ME_results$phylum == "Cyanobacteria")],OID = ME_results$MAG[which(ME_results$phylum == "Cyanobacteria")])
trace <- trace[order(trace$dates),]
trace2008 <- trace[which(trace$year == "2008"), ]
trace2009 <- trace[which(trace$year == "2009"), ]
trace2010 <- trace[which(trace$year == "2010"), ]
trace2011 <- trace[which(trace$year == "2011"), ]
trace2012 <- trace[which(trace$year == "2012"), ]

trace1 <- trace2008[which(trace2008$OID == "2582580518"), ]

trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], rollmean(trace1$input, 2))


p1 <- ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2008" & ME_results$MAG == "2582580518"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace1, aes(y = line, x = dates, color = OID), size = 2) + labs(title = "2008") + scale_color_manual(values = c("#a6611a")) + labs(x = NULL, y = NULL)

corrdata08 <- ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2008" & ME_results$MAG == "2582580518"), ]
corrdata08 <- corrdata08[order(corrdata08$Date), ]

trace4 <- trace2009[which(trace2009$OID == "2582580537"), ]

trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], rollmean(trace4$input, 2))

p2 <- ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2009" & ME_results$MAG == "2582580537"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace4, aes(y = line, x = dates, color = OID), size = 2) + labs(title = "2009") + scale_color_manual(values = c("#018571")) + labs(x = NULL, y = NULL)

corrdata09 <- ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2009" & ME_results$MAG == "2582580537"), ]
corrdata09 <- corrdata09[order(corrdata09$Date), ]

trace5 <- trace2010[which(trace2010$OID == "2582580551"), ]

trace5 <- trace5[order(trace5$dates),]
trace5$line <- c(trace5$input[1], rollmean(trace5$input, 2))


p3 <- ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2010" & ME_results$MAG == "2582580551"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace5, aes(y = line, x = dates, color = OID), size = 2) + labs(title = "2010") + scale_color_manual(values = c("#4dac26")) + labs(x = NULL, y = NULL)

corrdata10 <- ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2010" & ME_results$MAG == "2582580551"), ]
corrdata10 <- corrdata10[order(corrdata10$Date), ]

trace1 <- trace2011[which(trace2011$OID == "2582580518"), ]

trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], rollmean(trace1$input, 2))


p4 <- ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2011" & ME_results$MAG == "2582580518"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace1, aes(y = line, x = dates, color = OID), size = 2) + labs(title = "2011") + scale_color_manual(values = c("#a6611a")) + labs(x = NULL, y = NULL)

corrdata11 <- ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2011" & ME_results$MAG == "2582580518"), ]
corrdata11 <- corrdata11[order(corrdata11$Date), ]

trace10 <- trace2012[which(trace2012$OID == "2582580584"), ]

trace10 <- trace10[order(trace10$dates),]
trace10$line <- c(trace10$input[1], rollmean(trace10$input, 2))

p5 <- ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2012" & ME_results$MAG == "2582580584"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace10, aes(y = line, x = dates, color = OID), size = 2) + labs(title = "2012") + scale_color_manual(values = c("#d01c8b")) + labs(x = NULL, y = NULL)

corrdata12 <- ME_results[which(ME_results$phylum == "Cyanobacteria" & ME_results$year == "2012" & ME_results$MAG == "2582580584"), ]
corrdata12 <- corrdata12[order(corrdata12$Date), ]

all <- plot_grid(p1, p2, p3, p4, p5, nrow = 5, labels = c("A", "B", "C", "D", "E"))

# Add right panel of nitrogenase counts by year
table <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/marker_gene_table.txt", header = T, row.names = 1)
lakekey <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metagenome_metadata.csv", header = T)
lakerow <- as.character(lakekey$site[match(colnames(table), as.character(lakekey$sample))])
genedata <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metabolic_gene_info.txt", fill = TRUE)

genes2keep <- genedata$V1[grep("nitrogenase", as.character(genedata$V3))]
table2 <- table[match(genes2keep, rownames(table)), which(lakerow == "ME_epi")]
# Normalize by metagenome size
size_vector <- Mendota_raw_input$V4[match(colnames(table2), Mendota_raw_input$V2)]
table2 <- sweep(table2, 2, size_vector, "/")


table2$genes <- rownames(table2)
table2 <- melt(table2)
table2$date <- as.Date(Mendota_dates$date[match(table2$variable, Mendota_dates$sample)], format = "%m/%d/%y")
table2 <- table2[which(is.na(table2$date) == F),]
table2$Julian_Date <- as.numeric(format(table2$date, "%j"))
table2$year <- substr(table2$date, start = 1, stop = 4)

# Start plotting by year in the same format as above
trace <- table2
trace <- trace[order(trace$date),]
trace2008 <- trace[which(trace$year == "2008"), ]
trace2009 <- trace[which(trace$year == "2009"), ]
trace2010 <- trace[which(trace$year == "2010"), ]
trace2011 <- trace[which(trace$year == "2011"), ]
trace2012 <- trace[which(trace$year == "2012"), ]

table2 <- table2[which(table2$genes == "TIGR01282-sample1" | table2$genes == "TIGR01287-sample1" | table2$genes == "TIGR01286-sample1"),]
trace1 <- trace2008[which(trace2008$genes == "TIGR01286-sample1"), ]
trace2 <- trace2008[which(trace2008$genes == "TIGR01282-sample1"), ]
trace3 <- trace2008[which(trace2008$genes == "TIGR01287-sample1"), ]

trace1 <- trace1[order(trace1$date),]
trace1$line <- c(trace1$value[1], rollmean(trace1$value, 2))
trace2 <- trace2[order(trace2$date),]
trace2$line <- c(trace2$value[1], rollmean(trace2$value, 2))
trace3 <- trace3[order(trace3$date),]
trace3$line <- c(trace3$value[1], rollmean(trace3$value, 2))
trace <- rbind(trace1, trace2, trace3)


p6 <- ggplot(data = table2[which(table2$year == "2008"), ], aes(x = Julian_Date, y = value, color = genes)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + scale_x_continuous(breaks = pretty(table2$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = Julian_Date, color = genes), size = 2) + labs(title = "2008") + scale_color_manual(values = c("#b2182b", "#2166ac", "#542788")) + labs(x = NULL, y = NULL)

corrdata08$Gene1 <- trace1$value
corrdata08$Gene2 <- trace2$value
corrdata08$Gene3 <- trace3$value
cor(corrdata08$RPKM, corrdata08$Gene1)
cor(corrdata08$RPKM, corrdata08$Gene2)
cor(corrdata08$RPKM, corrdata08$Gene3)

trace1 <- trace2009[which(trace2009$genes == "TIGR01286-sample1"), ]
trace2 <- trace2009[which(trace2009$genes == "TIGR01282-sample1"), ]
trace3 <- trace2009[which(trace2009$genes == "TIGR01287-sample1"), ]

trace1 <- trace1[order(trace1$date),]
trace1$line <- c(trace1$value[1], rollmean(trace1$value, 2))
trace2 <- trace2[order(trace2$date),]
trace2$line <- c(trace2$value[1], rollmean(trace2$value, 2))
trace3 <- trace3[order(trace3$date),]
trace3$line <- c(trace3$value[1], rollmean(trace3$value, 2))
trace <- rbind(trace1, trace2, trace3)


p7 <- ggplot(data = table2[which(table2$year == "2009"), ], aes(x = Julian_Date, y = value, color = genes)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + scale_x_continuous(breaks = pretty(table2$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = Julian_Date, color = genes), size = 2) + labs(title = "2009") + scale_color_manual(values = c("#b2182b", "#2166ac", "#542788")) + labs(x = NULL, y = NULL)

corrdata09$Gene1 <- trace1$value
corrdata09$Gene2 <- trace2$value
corrdata09$Gene3 <- trace3$value
cor(corrdata09$RPKM, corrdata09$Gene1)
cor(corrdata09$RPKM, corrdata09$Gene2)
cor(corrdata09$RPKM, corrdata09$Gene3)

trace1 <- trace2010[which(trace2010$genes == "TIGR01286-sample1"), ]
trace2 <- trace2010[which(trace2010$genes == "TIGR01282-sample1"), ]
trace3 <- trace2010[which(trace2010$genes == "TIGR01287-sample1"), ]

trace1 <- trace1[order(trace1$date),]
trace1$line <- c(trace1$value[1], rollmean(trace1$value, 2))
trace2 <- trace2[order(trace2$date),]
trace2$line <- c(trace2$value[1], rollmean(trace2$value, 2))
trace3 <- trace3[order(trace3$date),]
trace3$line <- c(trace3$value[1], rollmean(trace3$value, 2))
trace <- rbind(trace1, trace2, trace3)


p8 <- ggplot(data = table2[which(table2$year == "2010"), ], aes(x = Julian_Date, y = value, color = genes)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + scale_x_continuous(breaks = pretty(table2$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = Julian_Date, color = genes), size = 2) + labs(title = "2010") + scale_color_manual(values = c("#b2182b", "#2166ac", "#542788")) + labs(x = NULL, y = NULL)

corrdata10$Gene1 <- trace1$value
corrdata10$Gene2 <- trace2$value
corrdata10$Gene3 <- trace3$value
cor(corrdata10$RPKM, corrdata10$Gene1)
cor(corrdata10$RPKM, corrdata10$Gene2)
cor(corrdata10$RPKM, corrdata10$Gene3)

trace1 <- trace2011[which(trace2011$genes == "TIGR01286-sample1"), ]
trace2 <- trace2011[which(trace2011$genes == "TIGR01282-sample1"), ]
trace3 <- trace2011[which(trace2011$genes == "TIGR01287-sample1"), ]

trace1 <- trace1[order(trace1$date),]
trace1$line <- c(trace1$value[1], rollmean(trace1$value, 2))
trace2 <- trace2[order(trace2$date),]
trace2$line <- c(trace2$value[1], rollmean(trace2$value, 2))
trace3 <- trace3[order(trace3$date),]
trace3$line <- c(trace3$value[1], rollmean(trace3$value, 2))
trace <- rbind(trace1, trace2, trace3)


p9 <- ggplot(data = table2[which(table2$year == "2011"), ], aes(x = Julian_Date, y = value, color = genes)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + scale_x_continuous(breaks = pretty(table2$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = Julian_Date, color = genes), size = 2) + labs(title = "2011") + scale_color_manual(values = c("#b2182b", "#2166ac", "#542788")) + labs(x = NULL, y = NULL)

corrdata11$Gene1 <- trace1$value
corrdata11$Gene2 <- trace2$value
corrdata11$Gene3 <- trace3$value
cor(corrdata11$RPKM, corrdata11$Gene1)
cor(corrdata11$RPKM, corrdata11$Gene2)
cor(corrdata11$RPKM, corrdata11$Gene3)

trace1 <- trace2012[which(trace2012$genes == "TIGR01286-sample1"), ]
trace2 <- trace2012[which(trace2012$genes == "TIGR01282-sample1"), ]
trace3 <- trace2012[which(trace2012$genes == "TIGR01287-sample1"), ]

trace1 <- trace1[order(trace1$date),]
trace1$line <- c(trace1$value[1], rollmean(trace1$value, 2))
trace2 <- trace2[order(trace2$date),]
trace2$line <- c(trace2$value[1], rollmean(trace2$value, 2))
trace3 <- trace3[order(trace3$date),]
trace3$line <- c(trace3$value[1], rollmean(trace3$value, 2))
trace <- rbind(trace1, trace2, trace3)


p10 <- ggplot(data = table2[which(table2$year == "2012"), ], aes(x = Julian_Date, y = value, color = genes)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + scale_x_continuous(breaks = pretty(table2$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = Julian_Date, color = genes), size = 2) + labs(title = "2012") + scale_color_manual(values = c("#b2182b", "#2166ac", "#542788")) + labs(x = NULL, y = NULL)

corrdata12$Gene1 <- trace1$value
corrdata12$Gene2 <- trace2$value
corrdata12$Gene3 <- trace3$value
cor(corrdata12$RPKM, corrdata12$Gene1)
cor(corrdata12$RPKM, corrdata12$Gene2)
cor(corrdata12$RPKM, corrdata12$Gene3)

all <- plot_grid(p1, p6, p2, p7, p3, p8, p4, p9, p5, p10, nrow = 5, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))


save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Manuscript_plots/Fig5.pdf", all, base_height = 10, base_aspect_ratio = 1.25)
