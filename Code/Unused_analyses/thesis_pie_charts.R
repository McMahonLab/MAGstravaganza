# MAGstravaganza thesis plots

#Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
############

# pie charts


# Get # of reads assigned to each phylum (mapped MAG) by lake-layer
Mendota_raw_input <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/Mendota_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
Mendota_dates <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
Mendota_dates$sample <- as.character(Mendota_dates$sample)
metadata <- read.csv("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

TE_raw_input <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/TE_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TE_dates <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
TE_dates$sample <- as.character(TE_dates$sample)

TH_raw_input <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/TH_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TH_dates <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")
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
metadata$Taxonomy <- gsub("Proteobacteria;", "", metadata$Taxonomy)
metadata$phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)

ME_results$phylum <- metadata$phylum[match(ME_results$MAG, metadata$IMG_OID)]
TE_results$phylum <- metadata$phylum[match(TE_results$MAG, metadata$IMG_OID)]
TH_results$phylum <- metadata$phylum[match(TH_results$MAG, metadata$IMG_OID)]

ME_results$Taxonomy <- metadata$Taxonomy[match(ME_results$MAG, metadata$IMG_OID)]
TE_results$Taxonomy <- metadata$Taxonomy[match(TE_results$MAG, metadata$IMG_OID)]
TH_results$Taxonomy <- metadata$Taxonomy[match(TH_results$MAG, metadata$IMG_OID)]

# Aggregate normalized reads by phylum
ME_phylum_results <- aggregate(RPKM ~ phylum, ME_results, sum)
TE_phylum_results <- aggregate(RPKM ~ phylum, TE_results, sum)
TH_phylum_results <- aggregate(RPKM ~ phylum, TH_results, sum)

# For each lake, make a pie chart of number of genes 

phylum_counts <- table(metadata$phylum[which(metadata$Lake == "Mendota")])
ME_pie <- ME_phylum_results
ME_pie$genomes <- as.numeric(phylum_counts)

phylum_counts <- table(metadata$phylum[which(metadata$Lake == "Trout Bog Epilimnion")])
TE_pie <- TE_phylum_results
TE_pie$genomes <- as.numeric(phylum_counts)

phylum_counts <- table(metadata$phylum[which(metadata$Lake == "Trout Bog Hypolimnion")])
TH_pie <- TH_phylum_results
TH_pie$genomes <- as.numeric(phylum_counts)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=18, face="bold")
  )

colors1 <- c("darkgrey", "#a50026", "#d6604d", "#fdae61",  "#e08214", "#fee08b", "#7fbc41", "#1a9850", "#abd9e9", "#4575b4", "#8073ac", "#b2abd2", "#542788")
colors2 <- c("#a50026", "#d6604d", "#fdae61",  "#e08214", "#fee08b", "#7fbc41", "#1a9850", "#abd9e9", "#8073ac", "#542788")

panelA <- ggplot(data = ME_pie, aes(x = "", y = genomes, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors1, name = NULL)
panelB <- ggplot(data = ME_pie, aes(x = "", y = RPKM, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors1) + theme(legend.position="none")

panelC <- ggplot(data = TE_pie, aes(x = "", y = genomes, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors2, name = NULL)
panelD <- ggplot(data = TE_pie, aes(x = "", y = RPKM, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors2) + theme(legend.position="none")

panelE <- ggplot(data = TH_pie, aes(x = "", y = genomes, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors1, name = NULL)
panelF <- ggplot(data = TH_pie, aes(x = "", y = RPKM, fill = phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = colors1) + theme(legend.position="none")

legend1 <- get_legend(panelA + theme_cowplot(font_size = 10))
legend2 <- get_legend(panelC + theme_cowplot(font_size = 10))
legend3 <- get_legend(panelE + theme_cowplot(font_size = 10))

panelA <- panelA + theme(legend.position="none")
panelC <- panelC + theme(legend.position="none")
panelE <- panelE + theme(legend.position="none")

pies <- plot_grid(panelA, panelB, legend1, panelC, panelD, legend2, panelE, panelF, legend3, labels = c("A. Mendota MAGs (99)", "B. Mendota Reads", NA, "C. Trout Bog Epi. MAGs (31)", "D. Trout Bog Epi. Reads", NA, "E. Trout Bog Hypo. MAGs (63)", "F. Trout Bog Hypo. Reads", NA), nrow = 3, label_size = 12, hjust = 0)
save_plot("C:/Users/Goose and Gander/Documents/magpies.png", pies, base_height = 8, base_aspect_ratio = 1)
