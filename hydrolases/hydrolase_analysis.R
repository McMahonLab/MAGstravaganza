# Analyze carbohydrate active enzymes in the MAGs

# Read in list of MAGs

mag_data <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)

gh_df <- mag_data[,c(1,3,4,10)]

gh_density <- c()
gh_diversity <- c()

cazy_MAGs <- rep(mag_data$IMG_OID, each = 7)
cazy_lakes <- rep(mag_data$Lake, each = 7)
cazy_names <- c("GH", "GT", "PL", "CE", "CB", "AA", "CS")
cazy_names_all <- rep(cazy_names, dim(mag_data)[1])
cazy_counts <- c()

for(i in 1:dim(mag_data)[1]){
  datafile <- read.table(paste("C:/Users/Alex/Desktop/MAGstravaganza/hydrolases/", mag_data$IMG_OID[i], ".txt", sep = ""))
  hits <- sub(".hmm", "", datafile$V2)
  types <- substr(hits, start = 1, stop = 2)
  
  add2vector <- c()
  add2vector[1] <- length(which(types == cazy_names[1]))/mag_data$Gene_Count[i]*100
  add2vector[2] <- length(which(types == cazy_names[2]))/mag_data$Gene_Count[i]*100
  add2vector[3] <- length(which(types == cazy_names[3]))/mag_data$Gene_Count[i]*100
  add2vector[4] <- length(which(types == cazy_names[4]))/mag_data$Gene_Count[i]*100
  add2vector[5] <- length(which(types == cazy_names[5]))/mag_data$Gene_Count[i]*100
  add2vector[6] <- length(which(types == cazy_names[6]))/mag_data$Gene_Count[i]*100
  add2vector[7] <- length(which(types == "do" | types == "co" | types == "SL"))/mag_data$Gene_Count[i]*100
  cazy_counts <- append(cazy_counts, add2vector, length(cazy_counts))
  
  gh_hits <- hits[which(types == "GH")]
  gh_density[i] <- length(gh_hits)/mag_data$Gene_Count[i]*100
  gh_diversity[i] <- length(unique(gh_hits))
}

gh_df$Diversity <- gh_diversity
gh_df$Density <- gh_density

cazy_df <- data.frame(cazy_MAGs, cazy_lakes, cazy_names_all, cazy_counts)
colnames(cazy_df) <- c("MAG", "Lake", "Enzyme", "Density")

gh_df$Order <- sapply(strsplit(mag_data$Taxonomy, ";"), "[", 3)
gh_df$Order <- factor(gh_df$Order, levels = rev(c("Actinomycetales", "Chlamydiales", "Cytophagales", "Mycoplasmatales", "Planctomycetales", "Puniceicoccales", "Rhodocyclales", "Sphingomonadales", "Xanthomonadales", "Solibacterales", "Bdellovibrionales", "Chlorobiales", "Holophagales", "Rickettsiales", "Bacteroidales", "Campylobacterales", "Desulfobacterales", "Desulfuromonadales", "Gallionellales", "Ignavibacteriales", "Legionellales", "Nitrosomonadales", "Pseudomonadales", "Rhizobiales", "Solirubrobacterales", "Flavobacteriales", "Acidimicrobiales", "Burkholderiales", "Methylococcales", "Methylophilales", "Sphingobacteriales", "Verrucomicrobiales")))
# Plot the cazy by lake

library(ggplot2)
library(cowplot)

ggplot(data = cazy_df, aes(x = Enzyme, y = Density, fill = Lake)) + geom_boxplot() + scale_fill_manual(values=c("#35978f", "#74add1", "#bf812d")) + background_grid(major = "xy", minor = "none")

ggplot(data = gh_df, aes(x = Density, y = Diversity, shape = Lake)) + geom_point(size = 3, alpha = 0.5) + scale_color_manual(values=c("#35978f", "#74add1", "#bf812d")) + background_grid(major = "xy", minor = "none")

ggplot(data = gh_df, aes(x = Order, y = Density, fill = Lake)) + geom_bar(stat = "summary", fun.y = "mean", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

[which(gh_df$Order == "Sphingobacteriales" | gh_df$Order == "Verrucomicrobiales" | gh_df$Order == "Bacteroidales" | gh_df$Order == "Flavobacteriales" | gh_df$Order == "Planctomycetales"), ]

ggplot(data = gh_df[which(gh_df$Order != "[Blank]" & is.na(gh_df$Order) == F), ], aes(x = Lake, y = Order, fill = Density)) + geom_tile() + scale_fill_gradient2(low = "#efedf5", mid = "#bcbddc", high = "#756bb1", midpoint = 4) + labs(title = "Mean Glucoside Hydrolase Coding Density", x = "", y = "")

ggplot(data = gh_df[which(gh_df$Order != "[Blank]" & is.na(gh_df$Order) == F), ], aes(x = Lake, y = Order, fill = Diversity)) + geom_tile() + scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#4575b4", midpoint = 40) + labs(title = "Mean Glucoside Hydrolase Coding Diversity", x = "", y = "")
