# convert Robin's R objects to csv format for input into OTUtable

TBtags <- readRDS("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/bogPyroTags_noFilter_tax_list.rds")
MEtags <- readRDS("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/metaGiTags_nofilter_list.rds")

# Get the otu table portions

TBtable <- TBtags$seqID
MEtable <- MEtags$seqID

# Remove Chloroplasts

TBtable <- TBtable[which(TBtable$class != "c__Chloroplast"), ]
MEtable <- MEtable[which(MEtable$class != "c__Chloroplast"), ]

# Remove singletons

persistence <- c()
for(i in 1:dim(TBtable)[1]){
  row <- TBtable[i, 9:dim(TBtable)[2]]
  persistence[i] <- length(which(row != 0))
}

TBtable <- TBtable[which(persistence > 1), ]

persistence <- c()
for(i in 1:dim(MEtable)[1]){
  row <- MEtable[i, 9:dim(MEtable)[2]]
  persistence[i] <- length(which(row != 0))
}
MEtable <- MEtable[which(persistence > 1), ]

# Output csv files

write.csv(TBtable, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Supplemental/TroutBog_OTUtable.csv", quote = F, row.names = F)
write.csv(MEtable, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Supplemental/Mendota_OTUtable.csv", quote = F, row.names = F)
