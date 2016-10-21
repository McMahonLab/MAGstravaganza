#Do patterns in genes truly repeat year after year?

#2008 is the weird year so leaving that out.

#Test correlation of genes in 2007 vs 2009
hypo.data <- read.csv("C:/Users/amlinz16/Dropbox/hypo_metagenome_genesurvey.csv", row.names=1)
hypo.dates <- as.Date(substr(colnames(hypo.data), start=2, stop=12), format = "%Y.%m.%d")

groups <- factor(substr(colnames(hypo.data), start = 4, stop = 5), levels = c("07", "08", "09"))

year07 <- hypo.data[,which(groups == "07")]
year08 <- hypo.data[,which(groups == "08")]
year09 <- hypo.data[,which(groups == "09")]

dates07 <- as.Date(substr(colnames(year07), start=2, stop=12), format = "%Y.%m.%d")
dates09 <- as.Date(substr(colnames(year09), start=2, stop=12), format = "%Y.%m.%d")
dates08 <- as.Date(substr(colnames(year08), start=2, stop=12), format = "%Y.%m.%d")

plot(dates07, year07[1,], type = "l")

epi.data <- read.csv("C:/Users/amlinz16/Dropbox/epi_metagenome_genesurvey.csv", row.names=1)
epi.dates <- as.Date(substr(colnames(epi.data), start=2, stop=12), format = "%Y.%m.%d")

find <- grep("cytochrome c oxidase", rownames(hypo.data))
cytc <- hypo.data[find,]
hypo.cytc <- colSums(cytc)/colSums(hypo.data)

find <- grep("cytochrome c oxidase", rownames(epi.data))
cytc <- epi.data[find,]
epi.cytc <- colSums(cytc)/colSums(hypo.data)

plot(hypo.dates, hypo.cytc, type = "l")
plot(epi.dates, epi.cytc, type= "l")
