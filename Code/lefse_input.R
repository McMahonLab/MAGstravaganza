# Process time series cog blast results
library(reshape2)
# Load in each file and tabulate the COG IDs

mendota <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
tbe <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
tbh <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")
lakekey <- read.csv("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/metagenome_metadata.csv", header = T)

metaG <- c(paste(mendota$sample, ".len150", sep = ""), as.character(tbe$sample), as.character(tbh$sample))
metaG <- metaG[which(metaG != "PTZH.len150" & metaG != "IHSC" & metaG != "IHTG" & metaG != "IHWH" & metaG != "IHWX" & metaG != "IHWP")]

ID <- c()
COGs <- c()
NUM <- c()
for(i in 1:length(metaG)){
  sample <- metaG[i]
  results_file <- read.table(paste("C:/Users/Goose and Gander/Documents/lefse_blast/", sample, ".blast", sep = ""))
  counts <- aggregate(V3 ~ V2, results_file, sum)
  ID <- append(ID, rep(sample, dim(counts)[1]), length(ID))
  COGs <- append(COGs, as.character(counts$V2), length(COGs))
  NUM <- append(NUM, counts$V3, length(NUM))
}
long_results <- data.frame(ID, COGs, NUM)

# DO NOT NORMALIZE. LEFSE DOES THIS
# # Divide counts by size of metagenome
# long_results$ID <- gsub(".len150", "", long_results$ID)
# all_lengths <- rbind(mendota[,c(1,2)], tbe[,c(1,2)], tbh[,c(1,2)])
# meta_size <- all_lengths$reads[match(long_results$ID, all_lengths$sample)]
# long_results$NUM <- long_results$NUM/meta_size

# Convert to wide format for input into LefSe

lefse_input <- dcast(long_results, COGs ~ ID)

# Replace NAs with 0
for (i in 1:dim(lefse_input)[1]){
  row <- lefse_input[i,]
  row[which(is.na(row) == T)] <- 0
  lefse_input[i,] <- row
}

# Add a 1st colum of the lake for each metagenome
# Make products a row. since output is tabular, replace spaces with underscores
rownames(lefse_input) <- lefse_input[,1]
lefse_input <- lefse_input[,2:dim(lefse_input)[2]]
colnames(lefse_input) <- gsub(".len150", "", colnames(lefse_input))
lakerow <- as.character(lakekey$site[match(colnames(lefse_input), as.character(lakekey$sample))])
lefse_input <- rbind(lakerow, lefse_input)
write.table(lefse_input, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/lefse_input.txt", quote = F, row.names = T, col.names = T, sep = "\t")

### karthik's version
#load lake datasets

table <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/marker_gene_table.txt", header = T, row.names = 1)
lakerow <- as.character(lakekey$site[match(colnames(table), as.character(lakekey$sample))])
lefse_input <- table[,2:dim(table)[2]]
write.table(lefse_input, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/alllefse_input.txt", quote = F, row.names = T, col.names = T, sep = "\t")
epi_input <- lefse_input[, which(lefse_input[1,] != "TB_hypo")]
write.table(lefse_input, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/epilefse_input.txt", quote = F, row.names = T, col.names = T, sep = "\t")
