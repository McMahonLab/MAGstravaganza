#Kraken analysis results

#Make generalizable so that I can re-run this on the GEODES results when those come back.

path2results <- "D:/chtc_kraken_files/"
#metadata_file <- 

#From the metadata file, make a list of files to upload
#samplenames <- metadata_file$samplenames
samplenames <- c("IHWP", "IHWZ", "PTXB.len150")

for (i in 1:length(samplenames)){
  results <- read.table(paste(path2results, samplenames[i], ".fasta.report", sep = ""), header = F, fill = T, sep = "\t")
  colnames(results) <- c("Percent_reads", "Num_reads_root", "Num_reads_direct", "Rank_code", "NCBI_taxonomy", "Scientific_name")
  assign(samplenames[i], results)
}

IHWP_phyla <- IHWP[which(IHWP$Rank_code == "P"), ]
IHWP_class <- IHWP[which(IHWP$Rank_code == "C"), ]
IHWZ[which(IHWZ$Rank_code == "C"), ]
PTXB.len150[which(PTXB.len150$Rank_code == "C"), ]
