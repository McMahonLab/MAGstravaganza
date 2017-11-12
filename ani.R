ani <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/ANI_matrix.csv", header = T, row.names = 1)

planctos <- c("2582580522", "2582580611", "2582580560", "2582580578", "2582580606", "2582580536", "2582580547", "2582580529", "2582580558", "2582580512", "2582580542", "2582580569", "2582580604")

mps <- c("2582580712", "2593339176", "2582580623", "2582580649", "2582580561", "2593339182", "2582580598", "2593339192")

mcs <- c("2582580677", "2593339178", "2582580639", "2582580566", "2582580615", "2582580646", "2593339191")

location <- c()
for(i in 1:length(mcs)){
  location[i] <- grep(mcs[i], rownames(ani))
}

mcs_ani <- ani[location, location]
heatmap(as.matrix(mcs_ani), Rowv = T, Colv = NA)

library(ape)
plancto_ani <- as.matrix(plancto_ani)
plancto_ani[which(plancto_ani == 0)] <- NA
tree <- njs(plancto_ani)
plot(tree)

library(ggtree)
tree <- NJ(plancto_ani)
print(tree)

tree <- hclust(as.dist(plancto_ani))
plot(tree)
