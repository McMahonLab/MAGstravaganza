# Question: are number of pathways related to abundance trends in the MAGs?
library(raster)
library(dplyr)
# Get number of pathways
path2repo <- "C:/Users/amlin/Desktop/MAGstravaganza/"

metapathways <- read.csv(paste(path2repo, "Data_files/Metapathways_output/master_table.csv", sep = ""), header = T)
path_info <- metapathways[2:203, 4:434]
rownames(path_info) <- substr(metapathways[2:203, 1], start = 2, stop = 11)

pathway_count <- rowSums(path_info  == 1)
names(pathway_count) <- rownames(path_info)

# Adjust pathway counts by the completeness of each mag
metadata <- read.csv(paste(path2repo, "Data_files/Genome_info/revised_MAG_metadata.csv", sep = ""), header = T)

metadata$Est_Completeness <- metadata$Est_Completeness/100
pathway_count2 <- pathway_count/(metadata$Est_Completeness[match(names(pathway_count), metadata$IMG_Genome_ID)])
hist(pathway_count2)

# Now I need info on the persistance, abundance, and CV of each MAG

TBH_coverage <- read.csv(paste(path2repo, "Previous_MAG_analyses/TBH_averaged_readcoverage.csv", sep = ""), header = T)
rownames(TBH_coverage) <- TBH_coverage[,1]
TBH_coverage <- TBH_coverage[, 2:58]

#I'm going to consider anything less than 1x coverage to be an absence
abundance <- c()
persistence <- c()
variance <- c()

for(j in 1:dim(TBH_coverage)[1]){
  row <- TBH_coverage[j, ]
  abundance[j] <- sum(row[which(row > 1)])/length(which(row > 1))
  persistence[j] <- length(which(row > 1))/length(row)
  variance[j] <- cv(as.numeric(row))
}

TBH_metrics <- data.frame(abundance, persistence, variance)
TBH_metrics$Genome <- substr(rownames(TBH_coverage), start = 1, stop = 7)

TBE_coverage <- read.csv(paste(path2repo, "Previous_MAG_analyses/TBE_averaged_readcoverage.csv", sep = ""), header = T)
rownames(TBE_coverage) <- TBE_coverage[,1]
TBE_coverage <- TBE_coverage[, 2:45]

abundance <- c()
persistence <- c()
variance <- c()

for(j in 1:dim(TBE_coverage)[1]){
  row <- TBE_coverage[j, ]
  abundance[j] <- sum(row[which(row > 1)])/length(which(row > 1))
  persistence[j] <- length(which(row > 1))/length(row)
  variance[j] <- cv(as.numeric(row))
}

TBE_metrics <- data.frame(abundance, persistence, variance)
TBE_metrics$Genome <- substr(rownames(TBE_coverage), start = 1, stop = 7)

TB_metrics <- rbind(TBE_metrics, TBH_metrics)

#I want to average the metrics of multiple contigs
TB_metrics <- group_by(TB_metrics, Genome)
TB_metrics <- summarize(TB_metrics, mean(abundance), mean(persistence), mean(variance))

# I have all the metrics I need - now I just need to match up the genome names in TB_metrics with the IMG ID in pathway_count2
TB_metrics$IMG_ID <- metadata$IMG_Genome_ID[match(TB_metrics$Genome, metadata$Genome_Name)]
TB_metrics <- TB_metrics[which(is.na(TB_metrics$IMG_ID) == F), ]
TB_metrics$pathway_count <- pathway_count2[match(TB_metrics$IMG_ID, names(pathway_count2))]

# Now start looking for correlations with number of pathways
plot(TB_metrics$`mean(abundance)`, TB_metrics$pathway_count)
plot(TB_metrics$`mean(variance)`, TB_metrics$pathway_count)
plot(TB_metrics$`mean(persistence)`, TB_metrics$pathway_count)

#Nothing obvious. Drop this line of work.