minreads <- 0
minsamples <- .01
lag <- 0
minlsa <- 0

library(OTUtable)
table <- read.csv(file.choose(), row.names=1)
table <- table[which(rowSums(table)> 10),]
dates <- as.Date(substr(colnames(table), start=2, stop=12), format = "%Y.%m.%d")
table <- table[,order(dates)]

temp <- zscore(table)

out <- "C:/Users/amlinz16/Documents/hypo_gene_network.txt"

write.table(temp, file="temp.txt", row.names=F, col.names=F, sep = "\t")
system(paste("C:/Users/amlinz16/fastlsa_win/fastlsa.exe -i C:/Users/amlinz16/Documents/temp.txt", " -o ", out, " -d ", lag, " -m ", minlsa, sep=""))

network <- read.table(file=out, header=T)
network$index1 <- rownames(table)[network$index1 + 1]
network$index2 <- rownames(table)[network$index2 + 1]

write.table(network, file=out, row.names=F, col.names=T, sep = "\t", quote=F)

network2 <- network[which(network$LSA >= 0.95),]
write.table(network2, file="C:/Users/amlinz16/Documents/hypo_gene_network.90.txt", row.names=F, col.names=T, sep = "\t", quote=F)
