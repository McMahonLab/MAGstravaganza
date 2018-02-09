# Analyze carbohydrate active enzymes in the MAGs

# Read in list of MAGs

mag_data <- read.csv("C:/Users/amlin/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)
mag_data$Lake <- as.character(mag_data$Lake)

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

##This one top contender for figure
gh.agg <- aggregate(Diversity ~ Lake + Order, data = gh_df, mean)
ggplot(data = gh.agg[which(gh.agg$Order != "[Blank]" & is.na(gh.agg$Order) == F), ], aes(x = Lake, y = Order, fill = Diversity)) + geom_tile() + scale_fill_gradient2(low = "white", mid = "#8c96c6", high = "#810f7c", midpoint = 40) + labs(title = "Mean Glucoside Hydrolase Coding Diversity", x = "", y = "") + background_grid(major = "xy")

gh.agg <- aggregate(Density ~ Lake + Order, data = gh_df, mean)
gh_density_plot <- ggplot(data = gh.agg[which(gh.agg$Order != "[Blank]" & is.na(gh.agg$Order) == F), ], aes(x = Lake, y = Order, fill = Density)) + geom_tile(color = "black") + scale_fill_gradient2(low = "white", mid = "#8c96c6", high = "#810f7c", midpoint = 4) + labs(x = "", y = "") + background_grid(major = "xy") + theme(axis.text.y = element_text(size = 10))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/Plots/GH.density.heatmap.pdf", gh_density_plot, base_height = 5, base_aspect_ratio = 1)
###

cazy_df$MAG <- factor(cazy_df$MAG)
ggplot(data = cazy_df[which(cazy_df$Enzyme == "GT"),], aes(x = Lake, y = MAG, fill = Density)) + geom_tile() + scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#4575b4", midpoint = 20) + labs(title = "GT density", x = "", y = "")


# Abundance of GHs
gh_names <- c()
gh_counts <- c()
magIDs <- c()
lakekey <- c()
taxonomy <- c()

for(i in 1:dim(mag_data)[1]){
  datafile <- read.table(paste("C:/Users/amlin/Desktop/MAGstravaganza/hydrolases/", mag_data$IMG_OID[i], ".txt", sep = ""))
  hits <- sub(".hmm", "", datafile$V2)
  types <- substr(hits, start = 1, stop = 2)
  gh_hits <- hits[which(types == "GH")]
  gh_table <- table(gh_hits)
  ghs <- names(gh_table)
  genome <- rep(mag_data$IMG_OID[i], length(ghs))
  lake <- rep(mag_data$Lake[i], length(ghs))
  tax <- rep(mag_data$Taxonomy[i], length(ghs))
  gh_names <- append(gh_names, ghs, length(gh_names))
  gh_counts <- append(gh_counts, gh_table, length(gh_counts))
  magIDs <- append(magIDs, genome, length(magIDs))
  lakekey <- append(lakekey, lake, length(lakekey))
  taxonomy <- append(taxonomy, tax, length(taxonomy))
}

gh_types <- data.frame(gh_names, gh_counts, magIDs, lakekey, taxonomy)
colnames(gh_types) <- c("GH", "Counts", "OID", "Lake", "Taxonomy")

#What are the 10 most aboundant GHs by lake?
gh_lake <- aggregate(Counts ~ GH + Lake, data = gh_types, sum)
ME <- gh_lake[which(gh_lake$Lake == "Mendota"),]
ME <- ME[order(ME$Counts, decreasing = T), ]
part1 <- ME[1:10,]

TE <- gh_lake[which(gh_lake$Lake == "Trout Bog Epilimnion"),]
TE <- TE[order(TE$Counts, decreasing = T), ]
part2 <- TE[1:10,]

TH <- gh_lake[which(gh_lake$Lake == "Trout Bog Hypolimnion"),]
TH <- TH[order(TH$Counts, decreasing = T), ]
part3 <- TH[1:10,]

#Make wide format to compare profiles
wide_gh_lake <- reshape(gh_lake, idvar = "Lake", timevar = "GH", direction = "wide")
for(i in 1:3){
  row <- wide_gh_lake[i, ]
  row[which(is.na(row) == T)] <- 0
  wide_gh_lake[i, ] <- row
}
rownames(wide_gh_lake) <- wide_gh_lake$Lake
wide_gh_lake <- wide_gh_lake[,2:177]
wide_gh_lake <- wide_gh_lake[,which(colSums(wide_gh_lake) > 50)]
heatmap(as.matrix(wide_gh_lake))

#Which hydrolases are unique to each region?

for(i in 1:dim(ME)[1]){
  inTE <- which(TE$GH == ME$GH[i])
  inTH <- which(TH$GH == ME$GH[i])
  if(length(inTE) == 0 && length(inTH) == 0){
    print(ME[i,])
  }
}
# unique (greater than 1) is largely GH13 families - 21, 5, 2, 2

for(i in 1:dim(TE)[1]){
  inME <- which(ME$GH == TE$GH[i])
  inTH <- which(TH$GH == TE$GH[i])
  if(length(inME) == 0 && length(inTH) == 0){
    print(TE[i,])
  }
}

# Only one is GH62 with 2 counts

for(i in 1:dim(TH)[1]){
  inTE <- which(TE$GH == TH$GH[i])
  inME <- which(ME$GH == TH$GH[i])
  if(length(inTE) == 0 && length(inME) == 0){
    print(TH[i,])
  }
}

#More than in the other lakes
#GH129 (11), GH43_12, GH89, GH44, GH66, GH67 (for ones greater 4 counts but there are many more)

wilcox.test(x = mag_data$Amino_Acid_Bias[which(mag_data$Lake == "Mendota")], y = mag_data$Amino_Acid_Bias[which(mag_data$Lake != "Mendota")], paired = F)

wilcox.test(x = mag_data$GC_Content[which(mag_data$Lake == "Mendota")], y = mag_data$GC_Content[which(mag_data$Lake != "Mendota")], paired = F)

mag_data$Est_Genome_Size <- mag_data$Genome_Size/mag_data$Est_Completeness
wilcox.test(x = mag_data$Est_Genome_Size[which(mag_data$Lake == "Mendota")], y = mag_data$Est_Genome_Size[which(mag_data$Lake != "Mendota")], paired = F)

# Make barplot of the top 10 GHs in each lake across lakes
# from line 112ish

top_gh <- c(as.character(part1$GH), as.character(part2$GH), as.character(part3$GH))
uniq_top_gh <- unique(top_gh)

#Make Mendota barplot
Mendota_data <- ME[match(uniq_top_gh, as.character(ME$GH)), ]
Mendota_data <- Mendota_data[order(Mendota_data$Counts, decreasing = T), ]
Mendota_data$GH <- factor(Mendota_data$GH, levels = Mendota_data$GH)

p1 <- ggplot(data = Mendota_data, aes(x = GH, y = Counts)) + geom_bar(stat = "identity", fill = "#FF0000", color = "black") + labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_y_continuous(limits = c(0, 670))
save_plot("C:/Users/amlin/Desktop/MAGstravaganza/Plots/ME_GH.pdf", p1, base_height = 4, base_aspect_ratio = 2.5) 

# Repeat with TB stacked bars

TH_data <- TH[match(uniq_top_gh, as.character(TH$GH)), ]
TE_data <- TE[match(uniq_top_gh, as.character(TE$GH)), ]
TB_data <- rbind(TH_data, TE_data)
TB_data <- TB_data[which(is.na(TB_data$GH) == F), ]
TB_data$GH <- factor(TB_data$GH, levels = Mendota_data$GH)

p2 <- ggplot(data = TB_data, aes(x = GH, y = Counts, fill = Lake)) + geom_bar(stat = "identity", color = "black") + labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(legend.position = "none") + scale_fill_manual(values = c("#8EEB00", "#00A287")) + scale_y_continuous(limits = c(0, 670))
save_plot("C:/Users/amlin/Desktop/MAGstravaganza/Plots/TH_GH.pdf", p2, base_height = 4, base_aspect_ratio = 2.5) 
