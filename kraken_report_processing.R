#Kraken analysis results
library(ggplot2)
library(cowplot)

#Make generalizable so that I can re-run this on the GEODES results when those come back.

metadata_file <- "C:/Users/Alex/Desktop/MAGstravaganza/Data_files/metagenome_metadata.csv"
path2results <- "D:/chtc_kraken_files/reports/"

#From the metadata file, make a list of files to upload
samplenames <- read.csv(metadata_file, header = T)
samples <- as.character(samplenames$sample)
samples[which(samplenames$site == "ME_epi")] <- paste(samples[which(samplenames$site == "ME_epi")], ".len150", sep = "")


for (i in 1:length(samples)){
  results <- read.table(paste(path2results, samples[i], ".fasta.report", sep = ""), header = F, sep = "\t")
  colnames(results) <- c("Percent_reads", "Num_reads_root", "Num_reads_direct", "Rank_code", "NCBI_taxonomy", "Scientific_name")
  assign(samples[i], results)
}
class_data <- CSNA.len150[which(CSNA.len150$Rank_code == "C" & CSNA.len150$Percent_reads > 0), c(1, 6)]
class_data$Sample <- rep(samples[1], dim(class_data)[1])

for(j in 2:length(samples)){
  datapart <- get(samples[i])
  sub_class_data <- datapart[which(datapart$Rank_code == "C" & datapart$Percent_reads > 0), c(1, 6)]
  sub_class_data$Sample <- rep(samples[j], dim(sub_class_data)[1])
  class_data <- rbind(class_data, sub_class_data)
}

class_data$Lake <- samplenames$site[match(class_data$Sample, samples)]
agg_class <- aggregate(Percent_reads ~ Lake + Scientific_name, data = class_data, FUN = sum)

num_ME <- length(which(samplenames$site == "ME_epi"))
num_TE <- length(which(samplenames$site == "TB_epi"))
num_TH <- length(which(samplenames$site == "TB_hypo"))

agg_class$Percent_reads[which(agg_class$Lake == "ME_epi")] <- agg_class$Percent_reads[which(agg_class$Lake == "ME_epi")]/num_ME
agg_class$Percent_reads[which(agg_class$Lake == "TB_epi")] <- agg_class$Percent_reads[which(agg_class$Lake == "TB_epi")]/num_TE
agg_class$Percent_reads[which(agg_class$Lake == "TB_hypo")] <- agg_class$Percent_reads[which(agg_class$Lake == "TB_hypo")]/num_TH

agg_class_reduced <- agg_class[which(agg_class$Percent_reads > 0.05), ]

ggplot(data = agg_class_reduced, aes(x = Lake, y = Percent_reads, fill = Scientific_name)) + geom_bar(stat = "identity")

# Repeat by order
order_data <- CSNA.len150[which(CSNA.len150$Rank_code == "O" & CSNA.len150$Percent_reads > 0), c(1, 6)]
order_data$Sample <- rep(samples[1], dim(order_data)[1])

for(j in 2:length(samples)){
  datapart <- get(samples[i])
  sub_order_data <- datapart[which(datapart$Rank_code == "O" & datapart$Percent_reads > 0), c(1, 6)]
  sub_order_data$Sample <- rep(samples[j], dim(sub_order_data)[1])
  order_data <- rbind(order_data, sub_order_data)
}

order_data$Lake <- samplenames$site[match(order_data$Sample, samples)]
agg_order <- aggregate(Percent_reads ~ Lake + Scientific_name, data = order_data, FUN = sum)
agg_order$Percent_reads[which(agg_order$Lake == "ME_epi")] <- agg_order$Percent_reads[which(agg_order$Lake == "ME_epi")]/num_ME
agg_order$Percent_reads[which(agg_order$Lake == "TB_epi")] <- agg_order$Percent_reads[which(agg_order$Lake == "TB_epi")]/num_TE
agg_order$Percent_reads[which(agg_order$Lake == "TB_hypo")] <- agg_order$Percent_reads[which(agg_order$Lake == "TB_hypo")]/num_TH

agg_order_reduced <- agg_order[which(agg_order$Percent_reads > 0.025), ]

myplot <- ggplot(data = agg_order_reduced, aes(x = Lake, y = Percent_reads, fill = Scientific_name)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Mean Percent Reads", title = "kraken classification") + theme(legend.title = element_blank())
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/plots/kraken_results.pdf", myplot, base_height = 3, base_aspect_ratio = 2)


