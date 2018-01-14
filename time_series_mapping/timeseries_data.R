# Plot abundance and variability of MAGs
library(raster)
library(ggplot2)
library(cowplot)
# Read data in
Mendota_raw_input <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/Mendota_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
Mendota_dates <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
Mendota_dates$sample <- as.character(Mendota_dates$sample)
metadata <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

TE_raw_input <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/TE_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TE_dates <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
TE_dates$sample <- as.character(TE_dates$sample)

TH_raw_input <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/TH_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TH_dates <- read.table("C:/Users/amlin/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")
TH_dates$sample <- as.character(TH_dates$sample)

# Convert to RPKM (reads per kilobase per million reads)
# general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
size_vector <- metadata$Genome_Size[match(Mendota_raw_input$V1, metadata$IMG_OID)]
Mendota_RPKM <- Mendota_raw_input$V3/(size_vector/1000 * Mendota_raw_input$V4/1000000)
ME_results <- data.frame(Mendota_raw_input[,1:2], Mendota_RPKM)
colnames(ME_results) <- c("MAG", "metaG", "RPKM")
ME_results$Date <- as.Date(Mendota_dates$date[match(ME_results$metaG, Mendota_dates$sample)], format = "%m/%d/%y")

size_vector <- metadata$Genome_Size[match(TE_raw_input$V1, metadata$IMG_OID)]
TE_RPKM <- TE_raw_input$V3/(size_vector/1000 * TE_raw_input$V4/1000000)
TE_results <- data.frame(TE_raw_input[,1:2], TE_RPKM)
colnames(TE_results) <- c("MAG", "metaG", "RPKM")
TE_results$Date <- as.Date(TE_dates$date[match(TE_results$metaG, TE_dates$sample)], format = "%m/%d/%y")

size_vector <- metadata$Genome_Size[match(TH_raw_input$V1, metadata$IMG_OID)]
TH_RPKM <- TH_raw_input$V3/(size_vector/1000 * TH_raw_input$V4/1000000)
TH_results <- data.frame(TH_raw_input[,1:2], TH_RPKM)
colnames(TH_results) <- c("MAG", "metaG", "RPKM")
TH_results$Date <- as.Date(TH_dates$date[match(TH_results$metaG, TH_dates$sample)], format = "%m/%d/%y")





# 
# # Calculate mean abundance and coefficient of variation
# 
# ME_MAGs <- unique(ME_results$MAG)
# abundance <- c()
# variation <- c()
# 
# for(i in 1:length(ME_MAGs)){
#   genome <- ME_MAGs[i]
#   data <- ME_results$RPKM[which(ME_results$MAG == genome)]
#   abundance[i] <- mean(data)
#   variation[i] <- cv(data)
# }
# 
# ME_traits <- data.frame(ME_MAGs, abundance, variation)
# 
# # Overlay other traits
# MAG_phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)
# 
# ME_rows <- match(ME_traits$ME_MAGs, metadata$IMG_OID)
# ME_traits$Genome_Size <- metadata$Genome_Size[ME_rows]/(metadata$Est_Completeness[ME_rows] / 100)
# ME_traits$Codon_Usage <- metadata$Codon_Usage[ME_rows]
# ME_traits$Amino_Acid_Bias <- metadata$Amino_Acid_Bias[ME_rows]
# ME_traits$Taxonomy <- metadata$Taxonomy[ME_rows]
# ME_traits$Phylum <- MAG_phylum[ME_rows]
# #Keep only phyla with > 5 genomes
# ME_traits <- ME_traits[which(ME_traits$Phylum == "Actinobacteria" | ME_traits$Phylum == "Bacteroidetes" | ME_traits$Phylum == "Cyanobacteria" | ME_traits$Phylum == "Planctomycetes" | ME_traits$Phylum == "Proteobacteria" | ME_traits$Phylum == "Verrucomicrobia"),]
# 
# # Plot
# #tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
# #rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
# rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
# #tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
# #dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
# #set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")
# 
# 
# 
# ggplot(data = ME_traits, aes(y = abundance, x = variation, color = Phylum)) + geom_point(size = 2)  + stat_density2d(geom = "polygon", aes(fill = Phylum), alpha = 0.1, bins = 4) + scale_color_manual(values = rich6equal)  + scale_fill_manual(values = rich6equal) + scale_y_continuous(limits = c(-0.01, 2)) + scale_x_continuous(limits = c(0, 700))
# 
# + scale_color_brewer(palette = "Pastel2") + scale_fill_brewer(palette = "Pastel2")
# 
# + scale_color_manual(values = c("darkgoldenrod4", "magenta3", "darkolivegreen3", "firebrick2", "orange1", "slateblue1")) + scale_fill_manual(values = c("darkgoldenrod4", "magenta3", "darkolivegreen3", "firebrick2", "orange1", "slateblue1"))
# 
# rich3equal = c("#000043", "#0033FF", "#FFCC00")
# ggplot(data = ME_traits[which(ME_traits$Phylum == "Actinobacteria" | ME_traits$Phylum == "Bacteroidetes" | ME_traits$Phylum == "Proteobacteria"),], aes(y = abundance, x = variation, color = Phylum)) + geom_point(size = 2) + scale_y_continuous(limits = c(-0.01, 0.15)) + scale_x_continuous(limits = c(0, 400)) + scale_color_manual(values = rich3equal)  + scale_fill_manual(values = rich3equal)+ stat_density2d(geom = "polygon", aes(fill = Phylum), alpha = 0.1, bins = 3)
# 
# #fix colors