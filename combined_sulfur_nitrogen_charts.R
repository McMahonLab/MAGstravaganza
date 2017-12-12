# Make panels of nitrogen and sulfur pathways by lake
library(ggplot2)
library(cowplot)
# Use data from Excel pathway analysis

table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/consolidated_pathway_data.csv", header = T, row.names = 2)

colnames(table) <- gsub("X", "", colnames(table))
table <- table[1:61,]
for(i in 2:dim(table)[2]){
  column <- table[,i]
  column[which(is.na(column) == T)] <- 0
  table[,i] <- column
}

# Read in metadata info to find out which MAGs came from which lake

metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

# Mendota Nitrogen
MEnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Nitrogen Degradation"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
MEnitrogen_totals <- rowSums(MEnitrogen > 0)
MEnitrogen_df <- data.frame(names(MEnitrogen_totals), MEnitrogen_totals)
colnames(MEnitrogen_df) <- c("Pathway", "MAGs")
MEN <- ggplot(data = MEnitrogen_df, aes(x = Pathway, y = MAGs/102*100)) + geom_bar(stat = "identity", fill = "#FF0000", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(95, 0))

#Trout Bog totals
TEnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Nitrogen Degradation"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TEnitrogen_totals <- rowSums(TEnitrogen > 0)
TEnitrogen_df <- data.frame(names(TEnitrogen_totals), TEnitrogen_totals)
colnames(TEnitrogen_df) <- c("Pathway", "MAGs")

THnitrogen <- table[which(table$Category == "Nitrogen Metabolism" | table$Category == "Nitrogen Degradation"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
THnitrogen_totals <- rowSums(THnitrogen > 0)
THnitrogen_df <- data.frame(names(THnitrogen_totals), THnitrogen_totals)
colnames(THnitrogen_df) <- c("Pathway", "MAGs")

TB_nitrogen_df <- rbind(TEnitrogen_df, THnitrogen_df)
TB_nitrogen_df$Layer <- rep(c("Epilimnion", "Hypolimnion"), each = 9)

TBN <- ggplot(data = TB_nitrogen_df, aes(x = Pathway, y = MAGs, fill = Layer)) + geom_bar(stat = "identity", color = "black") + labs(x = "", y = "") + coord_flip() + scale_fill_manual(values = c("#8EEB00", "#00A287"))  + theme(legend.position = "none")  + scale_y_continuous(limits = c(0, 95))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/ME_nitrogen.pdf", MEN, base_height = 4, base_aspect_ratio = 2)   
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/TB_nitrogen.pdf", TBN, base_height = 4, base_aspect_ratio = 2) 

#Repeat with sulfur
MEsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Mendota")])) == F)]
MEsulfur_totals <- rowSums(MEsulfur > 0)
MEsulfur_df <- data.frame(names(MEsulfur_totals), MEsulfur_totals)
colnames(MEsulfur_df) <- c("Pathway", "MAGs")
MES <- ggplot(data = MEsulfur_df, aes(x = Pathway, y = MAGs/102*100)) + geom_bar(stat = "identity", fill = "#FF0000", color = "black") + labs(x = "", y = "") + coord_flip() + scale_y_reverse(limits = c(95, 0))

#Trout Bog totals
TEsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Epilimnion")])) == F)]
TEsulfur_totals <- rowSums(TEsulfur > 0)
TEsulfur_df <- data.frame(names(TEsulfur_totals), TEsulfur_totals)
colnames(TEsulfur_df) <- c("Pathway", "MAGs")

THsulfur <- table[which(table$Category == "Sulfur Metabolism"), which(is.na(match(colnames(table), metadata$IMG_OID[which(metadata$Lake == "Trout Bog Hypolimnion")])) == F)]
THsulfur_totals <- rowSums(THsulfur > 0)
THsulfur_df <- data.frame(names(THsulfur_totals), THsulfur_totals)
colnames(THsulfur_df) <- c("Pathway", "MAGs")

TB_sulfur_df <- rbind(TEsulfur_df, THsulfur_df)
TB_sulfur_df$Layer <- rep(c("Epilimnion", "Hypolimnion"), each = 3)

TBS <- ggplot(data = TB_sulfur_df, aes(x = Pathway, y = MAGs, fill = Layer)) + geom_bar(stat = "identity", color = "black") + labs(x = "", y = "") + coord_flip() + scale_fill_manual(values = c("#8EEB00", "#00A287"))  + theme(legend.position = "none") + scale_y_continuous(limits = c(0, 95))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/ME_sulfur.pdf", MES, base_height = 2.5, base_aspect_ratio = 3.2)   
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/TB_sulfur.pdf", TBS, base_height = 2.5, base_aspect_ratio = 3.2) 

