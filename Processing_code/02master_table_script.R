# Goal: make a script to generate a master table that I can use for descriptive analysis

# For example, a question I would want to answer with this table would be "what phylogenetic groups contribute to photosynthesis in different lakes?"

# The general format would be genomes in rows, with pathways in columns.

# I want 2 additional columns with the lake the genome came from and its taxonomic assignment
# I want at least one additional row with a higher-level classification for each pathway.

path2repo <- "C:/Users/amlin/Desktop/MAGstravaganza/"

# Files - using curated sets
rxns <- read.csv(paste(path2repo, "curated_rxns.csv", sep=""), header = T, colClasses = c("character"))
metadata <- read.csv(paste(path2repo, "revised_MAG_metadata.csv", sep=""), colClasses = c("character", "numeric", "numeric", "character"))
pathways <- read.csv(paste(path2repo, "metacyc_pathway_tree.csv", sep=""), colClasses = c("character"))
taxonomy <- read.csv(paste(path2repo, "Readme.csv", sep = ""), colClasses = c("character"))

# Make the base of the table - a presence/absence table

genomes <- unique(rxns$SAMPLE)
pwys <- unique(rxns$PWY_COMMON_NAME)

pa_table <- matrix(0, length(genomes), length(pwys))
rownames(pa_table) <- genomes
colnames(pa_table) <- pwys

for(i in 1:length(genomes)){
  genome_pathways <- rxns[grep(genomes[i], rxns$SAMPLE),]
  for(j in 1:length(pwys)){
    x <- which(genome_pathways$PWY_COMMON_NAME == pwys[j])
    if(length(x > 0)){
      pa_table[i, j] <- 1
    }
  }
}

# Make a lake key column

lakekey <- substr(metadata$Genome_Name, start = 1, stop = 2)[match(substr(rownames(pa_table), start = 2, stop = 11), metadata$IMG_Genome_ID)]

# Make a taxonomy column

taxonomy_strings <- do.call(paste, c(taxonomy[, 3:8], sep = ";"))
hits <- match(substr(rownames(pa_table), start = 2, stop = 11), taxonomy$IMG.OID)
taxass <- taxonomy_strings[hits]

master_table <- data.frame(lakekey, taxass, pa_table)

# Make a top row that is level 3 in the metacyc tree file

pasted_pathways <- do.call(paste, c(pathways, sep = "_"))

# escape special regex characters in the common name

regex_pwys <- gsub("\\(", "\\\\(", colnames(pa_table))
regex_pwys <- gsub("\\[", "\\\\[", regex_pwys)
regex_pwys <- gsub("\\\\", "\\\\\\", regex_pwys)
regex_pwys <- gsub("<sup>", "", regex_pwys)
regex_pwys <- gsub("</sup>", "", regex_pwys)

level3 <- c()
missing <- c()
for(i in 1:dim(pa_table)[2]){
  x <- grep(regex_pwys[i], pasted_pathways)
  if(length(x) == 1){
    level3[i] <- pathways[x, 3]
  }else if (length(x) > 1){
    level3[i] <- pathways[x[1], 3]
  }else{
    missing <- append(missing, regex_pwys[i], length(missing))
  }
  if(i %% 100 == 0){
    print(i)
  }
}

level3 <- c("NA", "NA", level3)
master_table2 <- rbind(level3, master_table)

write.csv(master_table2, file = paste(path2repo, "master_table.csv"), row.names = T)
