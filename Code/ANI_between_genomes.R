# Which MAGs are 95% ANI?

#Load up the ANI matrix

ANImatrix <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/ANI_matrix.csv", header = T, row.names = 1)

# remove x and .fna from row and column names

rownames(ANImatrix) <- gsub(".fna", "", rownames(ANImatrix))
colnames(ANImatrix) <- gsub("X", "", colnames(ANImatrix))
colnames(ANImatrix) <- gsub(".fna", "", colnames(ANImatrix))

# For every MAG, find the other MAGs with greater than 95% ANI - excluding the match to itself

hits <- rep("none", dim(ANImatrix)[2])
for(i in 1:dim(ANImatrix)[2]){
  MAG <- ANImatrix[, i]
  sim <- which(MAG > 95)
  # Remove the match to itself
  sim <- sim[which(sim != i)]
  if(length(sim) > 0){
    names <- rownames(ANImatrix)[sim]
    hits[i] <- paste(names, collapse = ", ")
  }
}

#Output a csv file so I can copy/paste into my combined genome data file

write.csv(data.frame(colnames(ANImatrix), hits), "C:/Users/Alex/Desktop/MAGstravaganza/Data_files/ANI_95.csv")
