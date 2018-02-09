# Plot abundance and variability of MAGs
#library(raster)
library(ggplot2)
library(cowplot)
library(zoo)
# Read data in
Mendota_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/Mendota_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
Mendota_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/ME.sampledata.filtered.tsv", sep = "\t", header = T)
Mendota_dates$sample <- as.character(Mendota_dates$sample)
metadata <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)

TE_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/TE_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TE_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_epi.sample_data.filtered.txt", header = T, sep = "\t")
TE_dates$sample <- as.character(TE_dates$sample)

TH_raw_input <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/TH_results.txt", colClasses = c("character", "character", "numeric", "numeric"))
TH_dates <- read.table("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/tb_hyp.sample_data.filtered.txt", header = T, sep = "\t")
TH_dates$sample <- as.character(TH_dates$sample)

# Convert to RPKM (reads per kilobase per million reads)
# general form: RPKM =   numReads / ( genomeLength/1000 * totalNumReads/1,000,000 )
size_vector <- metadata$Genome_Size[match(Mendota_raw_input$V1, metadata$IMG_OID)]
Mendota_RPKM <- Mendota_raw_input$V3/(size_vector/1000 * Mendota_raw_input$V4/1000000)
ME_results <- data.frame(Mendota_raw_input[,1:2], Mendota_RPKM)
colnames(ME_results) <- c("MAG", "metaG", "RPKM")
ME_results$Date <- as.Date(Mendota_dates$date[match(ME_results$metaG, Mendota_dates$sample)], format = "%m/%d/%y")
ME_results$Julian_Date <- as.numeric(format(ME_results$Date, "%j"))
ME_results <- ME_results[which(is.na(ME_results$Date) == F),]

size_vector <- metadata$Genome_Size[match(TE_raw_input$V1, metadata$IMG_OID)]
TE_RPKM <- TE_raw_input$V3/(size_vector/1000 * TE_raw_input$V4/1000000)
TE_results <- data.frame(TE_raw_input[,1:2], TE_RPKM)
colnames(TE_results) <- c("MAG", "metaG", "RPKM")
TE_results$Date <- as.Date(TE_dates$date[match(TE_results$metaG, TE_dates$sample)], format = "%m/%d/%y")
TE_results$Julian_Date <- as.numeric(format(TE_results$Date, "%j"))
TE_results <- TE_results[which(is.na(TE_results$Date) == F),]

size_vector <- metadata$Genome_Size[match(TH_raw_input$V1, metadata$IMG_OID)]
TH_RPKM <- TH_raw_input$V3/(size_vector/1000 * TH_raw_input$V4/1000000)
TH_results <- data.frame(TH_raw_input[,1:2], TH_RPKM)
colnames(TH_results) <- c("MAG", "metaG", "RPKM")
TH_results$Date <- as.Date(TH_dates$date[match(TH_results$metaG, TH_dates$sample)], format = "%m/%d/%y")
TH_results$Julian_Date <- as.numeric(format(TH_results$Date, "%j"))
TH_results <- TH_results[which(is.na(TH_results$Date) == F),]

metadata$phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)
metadata$class<- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 2)
metadata$order <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 3)

ME_results$phylum <- metadata$phylum[match(ME_results$MAG, metadata$IMG_OID)]
ME_results$class <- metadata$class[match(ME_results$MAG, metadata$IMG_OID)]
ME_results$order <- metadata$order[match(ME_results$MAG, metadata$IMG_OID)]

TE_results$phylum <- metadata$phylum[match(TE_results$MAG, metadata$IMG_OID)]
TE_results$class <- metadata$class[match(TE_results$MAG, metadata$IMG_OID)]
TE_results$order <- metadata$order[match(TE_results$MAG, metadata$IMG_OID)]

TH_results$phylum <- metadata$phylum[match(TH_results$MAG, metadata$IMG_OID)]
TH_results$class <- metadata$class[match(TH_results$MAG, metadata$IMG_OID)]
TH_results$order <- metadata$order[match(TH_results$MAG, metadata$IMG_OID)]
# Make traces of select groups by lake
# Primary producers - Cyano in Mendota, Chlorobi in Trout Bog

# ggplot(data = ME_results[which(ME_results$phylum == "Cyanobacteria"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point() + scale_y_continuous(limits = c(0, 30))
# ggplot(data = TE_results[which(TE_results$phylum == "Chlorobi"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# ggplot(data = TH_results[which(TH_results$phylum == "Chlorobi"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# # Conclusion: lots of cyano strains all over the place. One dominant Chlorobi that increases towards fall
# 
# # Glycoside hydrolasers
# ggplot(data = ME_results[which(ME_results$phylum == "Bacteroidetes"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point() #there are 40ish MAGs!
# ggplot(data = ME_results[which(ME_results$phylum == "Verrucomicrobia"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# ggplot(data = ME_results[which(ME_results$phylum == "Planctomycetes"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# 
# ggplot(data = TE_results[which(TE_results$phylum == "Bacteroidetes"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# ggplot(data = TE_results[which(TE_results$phylum == "Verrucomicrobia"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point() + geom_line()
# 
# ggplot(data = TH_results[which(TH_results$phylum == "Bacteroidetes"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()
# ggplot(data = TH_results[which(TH_results$phylum == "Verrucomicrobia"), ], aes(x = Date, y = RPKM, color = MAG)) + geom_point()

#These are all pretty terrible

# Methylotrophs
seasons <- data.frame(x1 = c(75, 152, 244), x2 = c(151, 243, 334), y1 = c(0, 0, 0), y2 = c(Inf, Inf, Inf))
trace <- data.frame(input = ME_results$RPKM[which(ME_results$order == "Methylococcales")], dates = ME_results$Julian_Date[which(ME_results$order == "Methylococcales")])
trace <- trace[order(trace$dates),]
trace$line <- c(trace$input[1], mean(trace$input[1:2]), mean(trace$input[1:3]), mean(trace$input[1:4]), mean(trace$input[1:5]), mean(trace$input[1:6]), mean(trace$input[1:7]), mean(trace$input[1:8]), rollmean(trace$input, 9))
trace$MAG <- rep("2582580566", length(trace$line))

p1 <- ggplot(data = ME_results[which(ME_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = ME_results$RPKM[which(ME_results$order == "Methylophilales")], dates = ME_results$Julian_Date[which(ME_results$order == "Methylophilales")], MAG =  ME_results$MAG[which(ME_results$order == "Methylophilales")])
trace1 <- trace[which(trace$MAG == "2582580561"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580598"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace <- rbind(trace1, trace2)
            
p2 <- ggplot(data = ME_results[which(ME_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TE_results$RPKM[which(TE_results$order == "Methylococcales")], dates = TE_results$Julian_Date[which(TE_results$order == "Methylococcales")], MAG =  TE_results$MAG[which(TE_results$order == "Methylococcales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580615"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580639"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2582580646"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace <- rbind(trace1, trace2, trace3)

p3 <- ggplot(data = TE_results[which(TE_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TE_results$RPKM[which(TE_results$order == "Methylophilales")], dates = TE_results$Julian_Date[which(TE_results$order == "Methylophilales")], MAG =  TE_results$MAG[which(TE_results$order == "Methylophilales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580623"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580649"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace <- rbind(trace1, trace2)

p4 <- ggplot(data = TE_results[which(TE_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TH_results$RPKM[which(TH_results$order == "Methylococcales")], dates = TH_results$Julian_Date[which(TH_results$order == "Methylococcales")], MAG =  TH_results$MAG[which(TH_results$order == "Methylococcales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580677"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2593339178"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2593339191"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace <- rbind(trace1, trace2, trace3)

p5 <- ggplot(data = TH_results[which(TH_results$order == "Methylococcales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

trace <- data.frame(input = TH_results$RPKM[which(TH_results$order == "Methylophilales")], dates = TH_results$Julian_Date[which(TH_results$order == "Methylophilales")], MAG =  TH_results$MAG[which(TH_results$order == "Methylophilales")])
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580712"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2593339176"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2593339182"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace4 <- trace[which(trace$MAG == "2593339192"), ]
trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], mean(trace4$input[1:2]), mean(trace4$input[1:3]), mean(trace4$input[1:4]), mean(trace4$input[1:5]), mean(trace4$input[1:6]), mean(trace4$input[1:7]), mean(trace4$input[1:8]), rollmean(trace4$input, 9))
trace <- rbind(trace1, trace2, trace3, trace4)
               
p6 <- ggplot(data = TH_results[which(TH_results$order == "Methylophilales"), ], aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334)) + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")

meth <- plot_grid(p1, p3, p5, p2, p4, p6, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/methylotroph_trace.pdf", meth, base_height = 5, base_aspect_ratio = 3)

# Much better
# one last thing to look at - freshwater favorites
# acI is the best bet - no Pnecs in ME or Limnohabitans anywhere
acI <- c("ME00283", "ME03864", "ME04252", "ME00885", "TE02057v2", "TE02754", "TE03207", "TE03475", "TE04208", "TH02152", "TH03463", "TH03838v2", "TH00680")
acI_MAGs <- metadata$IMG_OID[match(acI, metadata$Genome_Name)]
acI_ME <- ME_results[which(is.na(match(ME_results$MAG, acI_MAGs[1:4])) == F),]
acI_TE <- TE_results[which(is.na(match(TE_results$MAG, acI_MAGs[5:9])) == F),]
acI_TH <- TH_results[which(is.na(match(TH_results$MAG, acI_MAGs[10:13])) == F),]

trace <- data.frame(input = acI_ME$RPKM, dates = acI_ME$Julian_Date, MAG =  acI_ME$MAG)
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580555"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580572"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2582580575"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace4 <- trace[which(trace$MAG == "2582580608"), ]
trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], mean(trace4$input[1:2]), mean(trace4$input[1:3]), mean(trace4$input[1:4]), mean(trace4$input[1:5]), mean(trace4$input[1:6]), mean(trace4$input[1:7]), mean(trace4$input[1:8]), rollmean(trace4$input, 9))
trace <- rbind(trace1, trace2, trace3, trace4)

p1 <- ggplot(data = acI_ME, aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2")  + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334))

trace <- data.frame(input = acI_TE$RPKM, dates = acI_TE$Julian_Date, MAG =  acI_TE$MAG)
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580627"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580630"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2582580632"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace4 <- trace[which(trace$MAG == "2582580635"), ]
trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], mean(trace4$input[1:2]), mean(trace4$input[1:3]), mean(trace4$input[1:4]), mean(trace4$input[1:5]), mean(trace4$input[1:6]), mean(trace4$input[1:7]), mean(trace4$input[1:8]), rollmean(trace4$input, 9))
trace <- rbind(trace1, trace2, trace3, trace4)

p2 <- ggplot(data = acI_TE, aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334))

trace <- data.frame(input = acI_TH$RPKM, dates = acI_TH$Julian_Date, MAG =  acI_TH$MAG)
trace <- trace[which(is.na(trace$input) == F), ]
trace1 <- trace[which(trace$MAG == "2582580657"), ]
trace1 <- trace1[order(trace1$dates),]
trace1$line <- c(trace1$input[1], mean(trace1$input[1:2]), mean(trace1$input[1:3]), mean(trace1$input[1:4]), mean(trace1$input[1:5]), mean(trace1$input[1:6]), mean(trace1$input[1:7]), mean(trace1$input[1:8]), rollmean(trace1$input, 9))
trace2 <- trace[which(trace$MAG == "2582580675"), ]
trace2 <- trace2[order(trace2$dates),]
trace2$line <- c(trace2$input[1], mean(trace2$input[1:2]), mean(trace2$input[1:3]), mean(trace2$input[1:4]), mean(trace2$input[1:5]), mean(trace2$input[1:6]), mean(trace2$input[1:7]), mean(trace2$input[1:8]), rollmean(trace2$input, 9))
trace3 <- trace[which(trace$MAG == "2582580705"), ]
trace3 <- trace3[order(trace3$dates),]
trace3$line <- c(trace3$input[1], mean(trace3$input[1:2]), mean(trace3$input[1:3]), mean(trace3$input[1:4]), mean(trace3$input[1:5]), mean(trace3$input[1:6]), mean(trace3$input[1:7]), mean(trace3$input[1:8]), rollmean(trace3$input, 9))
trace4 <- trace[which(trace$MAG == "2593339186"), ]
trace4 <- trace4[order(trace4$dates),]
trace4$line <- c(trace4$input[1], mean(trace4$input[1:2]), mean(trace4$input[1:3]), mean(trace4$input[1:4]), mean(trace4$input[1:5]), mean(trace4$input[1:6]), mean(trace4$input[1:7]), mean(trace4$input[1:8]), rollmean(trace4$input, 9))
trace <- rbind(trace1, trace2, trace3, trace4)

p3 <- ggplot(data = acI_TH, aes(x = Julian_Date, y = RPKM, color = MAG)) + geom_point(size = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + labs(x = "Julian Date") + geom_rect(data = seasons, inherit.aes = FALSE, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), alpha = 0.1, fill = c("skyblue", "springgreen", "goldenrod"), color = "black") + geom_line(data = trace, aes(y = line, x = dates), size = 2) + scale_color_brewer(palette = "Set2") + scale_x_continuous(breaks = pretty(ME_results$Julian_Date, 30), limits = c(75, 334))

acI_plot <- plot_grid(p1, p2, p3, nrow = 3, labels = c("A", "B", "C"))
save_plot("C:/Users/Alex/Desktop/MAGstravaganza/time_series_mapping/acI_trace.pdf", acI_plot, base_height = 8, base_aspect_ratio = 1)

# 
# # Calculate mean abundance and coefficient of variation
# 
# ME_MAGs <- unique(ME_results$MAG)
# abundance <- c()
# variation <- c()
# 
# for(i in 1:length(ME_MAGs)){
#   genome <- ME_MAGs[i]
#   data <- ME_results$RPKM[which(ME_results$MAG == genome)]
#   abundance[i] <- mean(data)
#   variation[i] <- cv(data)
# }
# 
# ME_traits <- data.frame(ME_MAGs, abundance, variation)
# 
# # Overlay other traits
# MAG_phylum <- sapply(strsplit(as.character(metadata$Taxonomy),";"), `[`, 1)
# 
# ME_rows <- match(ME_traits$ME_MAGs, metadata$IMG_OID)
# ME_traits$Genome_Size <- metadata$Genome_Size[ME_rows]/(metadata$Est_Completeness[ME_rows] / 100)
# ME_traits$Codon_Usage <- metadata$Codon_Usage[ME_rows]
# ME_traits$Amino_Acid_Bias <- metadata$Amino_Acid_Bias[ME_rows]
# ME_traits$Taxonomy <- metadata$Taxonomy[ME_rows]
# ME_traits$Phylum <- MAG_phylum[ME_rows]
# #Keep only phyla with > 5 genomes
# ME_traits <- ME_traits[which(ME_traits$Phylum == "Actinobacteria" | ME_traits$Phylum == "Bacteroidetes" | ME_traits$Phylum == "Cyanobacteria" | ME_traits$Phylum == "Planctomycetes" | ME_traits$Phylum == "Proteobacteria" | ME_traits$Phylum == "Verrucomicrobia"),]
# 
# # Plot
# #tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
# #rainbow6equal = c("#BF4D4D", "#BFBF4D", "#4DBF4D", "#4DBFBF", "#4D4DBF", "#BF4DBF")
# rich6equal = c("#000043", "#0033FF", "#01CCA4", "#BAFF12", "#FFCC00", "#FF3300")
# #tim6equal = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")
# #dark6equal = c("#1B9E77", "#66A61E", "#7570B3", "#D95F02", "#E6AB02", "#E7298A")
# #set6equal = c("#66C2A5", "#8DA0CB", "#A6D854", "#E78AC3", "#FC8D62", "#FFD92F")
# 
# 
# 
# ggplot(data = ME_traits, aes(y = abundance, x = variation, color = Phylum)) + geom_point(size = 2)  + stat_density2d(geom = "polygon", aes(fill = Phylum), alpha = 0.1, bins = 4) + scale_color_manual(values = rich6equal)  + scale_fill_manual(values = rich6equal) + scale_y_continuous(limits = c(-0.01, 2)) + scale_x_continuous(limits = c(0, 700))
# 
# + scale_color_brewer(palette = "Pastel2") + scale_fill_brewer(palette = "Pastel2")
# 
# + scale_color_manual(values = c("darkgoldenrod4", "magenta3", "darkolivegreen3", "firebrick2", "orange1", "slateblue1")) + scale_fill_manual(values = c("darkgoldenrod4", "magenta3", "darkolivegreen3", "firebrick2", "orange1", "slateblue1"))
# 
# rich3equal = c("#000043", "#0033FF", "#FFCC00")
# ggplot(data = ME_traits[which(ME_traits$Phylum == "Actinobacteria" | ME_traits$Phylum == "Bacteroidetes" | ME_traits$Phylum == "Proteobacteria"),], aes(y = abundance, x = variation, color = Phylum)) + geom_point(size = 2) + scale_y_continuous(limits = c(-0.01, 0.15)) + scale_x_continuous(limits = c(0, 400)) + scale_color_manual(values = rich3equal)  + scale_fill_manual(values = rich3equal)+ stat_density2d(geom = "polygon", aes(fill = Phylum), alpha = 0.1, bins = 3)
# 
# #fix colors