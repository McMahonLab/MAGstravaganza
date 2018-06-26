# Make figure of genome size vs number of genomes, filled by phylum
library(ggplot2)
library(cowplot)
genome_data <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Data_files/Genome_info/MAG_phylogeny.csv", header = T)

taxonomy <- strsplit(as.character(genome_data$Taxonomy.README.), ";")
phylum <- c()
for(i in 1:length(taxonomy)){
  phylum[i] <- taxonomy[[i]][1]
}

genome_size <- genome_data$Genome_Size / (genome_data$Est_Completeness/100)
size_category <- c()
size_category[which(genome_size <= quantile(genome_size)[2])] <- "1st"
size_category[which(genome_size <= quantile(genome_size)[3] & genome_size > quantile(genome_size)[2])] <- "2nd"
size_category[which(genome_size <= quantile(genome_size)[4] & genome_size > quantile(genome_size)[3])] <- "3rd"
size_category[which(genome_size <= quantile(genome_size)[5]  & genome_size > quantile(genome_size)[4])] <- "4th"


figure_data <- data.frame(genome_size, size_category, phylum)
figure_data <- figure_data[which(figure_data$phylum == "[Blank]" | figure_data$phylum == "Acidobacteria" | figure_data$phylum == "Actinobacteria" | figure_data$phylum == "Bacteroidetes" | figure_data$phylum == "Cyanobacteria" | figure_data$phylum == "Planctomycetes" | figure_data$phylum == "Proteobacteria" | figure_data$phylum == "Verrucomicrobia"), ]
figure_data$phylum <- as.character(figure_data$phylum)
figure_data$phylum[which(figure_data$phylum == "[Blank]")] <- "Unclassified"

ggplot(data = figure_data) + geom_bar(aes(x = size_category, fill = phylum)) + scale_fill_brewer(palette = "Set1") + labs(x = "Quartile by Estimated Genome Size", y = "Number of Genomes")

ggplot(data = figure_data) + geom_histogram(aes(x = genome_size, fill = phylum), binwidth = 2000000) + scale_fill_brewer(palette = "Set1") + labs(x = "Estimated Genome Size", y = "Number of Genomes")
