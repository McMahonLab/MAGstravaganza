# Identify hubs of connectivity  by performing network analysis on the read coverage data for Trout bog
library(GGally)
library(network)
library(sna)
library(ggplot2)
path2repo <- "C:/Users/Alex/Desktop/MAGstravaganza/"
setwd(path2repo)

TBH_coverage <- read.csv(paste(path2repo, "Previous_MAG_analyses/TBH_averaged_readcoverage.csv", sep = ""), header = T)
rownames(TBH_coverage) <- TBH_coverage[,1]
TBH_coverage <- TBH_coverage[, 2:58]

TBE_coverage <- read.csv(paste(path2repo, "Previous_MAG_analyses/TBE_averaged_readcoverage.csv", sep = ""), header = T)
rownames(TBE_coverage) <- TBE_coverage[,1]
TBE_coverage <- TBE_coverage[, 2:45]

minreads <- 1
minlsa <- .75

lsa_prep <- function(table){
  ztable <- matrix(NA, ncol=length(table[1,]), nrow=length(table[,1]))
  for (i in 1:length(table[,1])){
    otu <- as.numeric(table[i,])
    zotu <- (otu- mean(otu))/sd(otu)
    ztable[i,] <- zotu
  }
  return(ztable)
}

temp <- lsa_prep(TBE_coverage)
write.table(temp, file="temp.txt", row.names=F, col.names=F, sep = "\t")
out <- "TBE_MAGS_network_25Oct16.txt"

#My desktop path
#system(paste("D:/fastlsa_win/fastlsa.exe -i temp.txt", " -o ", out, " -d 0", " -m ", minlsa, sep=""))
#My laptop path
system(paste("C:/Users/Alex/fastlsa_win/fastlsa.exe -i temp.txt", " -o ", out, " -d 0", " -m ", minlsa, sep=""))

network <- read.table(file=out, header=T)
network$index1 <- rownames(TBE_coverage)[network$index1 + 1]
network$index2 <- rownames(TBE_coverage)[network$index2 + 1]
write.table(network, file=out, row.names=F, col.names=T, sep = "\t", quote=F)

#Number edges:
dim(network)[1]
#Number nodes:
length(unique(c(network$index1, network$index2)))
network$index1[which(table(network$index1) > 10)]
#Greater than 10 connections in TBE - 0031, 0439, 0037, 0567
#IMG numbers: the only one in my metadata file is 567 2582580642
#It's a 70% complete Burkholderiales



#x <- network(network)
#ggnet(x, node.size = 0.5)

temp <- lsa_prep(TBH_coverage)
write.table(temp, file="temp.txt", row.names=F, col.names=F, sep = "\t")
out <- "TBH_MAGS_network_25Oct16.txt"

system(paste("C:/Users/Alex/fastlsa_win/fastlsa.exe -i temp.txt", " -o ", out, " -d 0", " -m ", minlsa, sep=""))

network <- read.table(file=out, header=T)
network$index1 <- rownames(TBH_coverage)[network$index1 + 1]
network$index2 <- rownames(TBH_coverage)[network$index2 + 1]
write.table(network, file=out, row.names=F, col.names=T, sep = "\t", quote=F)
network$index1[which(table(network$index1) > 40)]

