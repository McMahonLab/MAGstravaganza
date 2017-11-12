#Categories by lake 
library(reshape2)
library(dplyr)
library(ggplot2)
table <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/rxn_categories_by_lake.csv", header = T, colClasses = "character")
categores <- unlist(table[1,], use.names = F)
table2 <- table[2:dim(table)[1], ]
#table2$X <- as.character(table2$X)
table2[, 2:dim(table2)[2]] <- as.numeric(as.matrix(table2[, 2:dim(table2)[2]]))
table3 <- melt(table2)
table3$cat <- categores[match(table3$variable, colnames(table))]
table3 <- group_by(table3, X, cat)
sum_table <- summarize(table3, sums = sum(value))

ME_total <- sum(sum_table$sums[which(sum_table$X == "ME")])
TE_total <- sum(sum_table$sums[which(sum_table$X == "TE")])
TH_total <- sum(sum_table$sums[which(sum_table$X == "TH")])

sum_table$sums[which(sum_table$X == "ME")] <- sum_table$sums[which(sum_table$X == "ME")]/ME_total
sum_table$sums[which(sum_table$X == "TE")] <- sum_table$sums[which(sum_table$X == "TE")]/TE_total
sum_table$sums[which(sum_table$X == "TH")] <- sum_table$sums[which(sum_table$X == "TH")]/TH_total

sum_table$cat <- toupper(sum_table$cat)

ggplot(sum_table, aes(x = X, y = cat, fill = sums)) + geom_tile() + theme_bw() + labs(x = NULL, y = NULL) + scale_fill_gradient2(low = "white", mid = "goldenrod", high = "darkred", midpoint = 0.10)




