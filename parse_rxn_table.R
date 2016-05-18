# Path - change as needed
path2repo <- "C:/Users/Alex/Desktop/MAGstravaganza/"

# Files
coverage <- read.table(paste(path2repo, "combined_MAGS_pathwaycoverage.txt", sep=""), sep = "\t", header = T, colClasses = c(rep("character", 3), rep("numeric", 3)))
rxns <- read.table(paste(path2repo, "combinedMAGS_pws_rxns.orfs.txt", sep=""), sep = "\t", header = T, colClasses = c("character"))
metadata <- read.csv(paste(path2repo, "revised_MAG_metadata.csv", sep=""), colClasses = c("character", "numeric", "numeric", "character"))
pathways <- read.csv(paste(path2repo, "redox_pathways.csv", sep=""), colClasses = c("character"))

# Remove genomes that are not in the metadata file, and remove pathways with less than 50% coverage

rxns2 <- rxns[which(substr(rxns$SAMPLE, start = 2, stop = 11) %in% metadata$IMG_Genome_ID),]

coverage$percent <- coverage$NUM_COVERED_REACTIONS / coverage$NUM_REACTIONS * 100
hits <- match(paste(rxns2$SAMPLE, rxns2$PWY_NAME, sep = ""), paste(coverage$SAMPLE, coverage$PWY_NAME, sep = ""))
rxns2$coverage <- coverage$percent[hits]
rxns3 <- rxns2[which(rxns2$coverage > 50),]

write.csv(rxns3, file = paste(path2repo, "curated_rxns.csv", sep = ""), row.names = F)

rxnkey <- unique(rxns3$PWY_COMMON_NAME)

write.csv(rxnkey, file = paste(path2repo, "rxn_key.csv", sep = ""))
