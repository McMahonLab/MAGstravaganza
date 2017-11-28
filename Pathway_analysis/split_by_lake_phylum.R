# Split MAGs IMG pathway info by lake and phylum
# Make sure there are no more than 30 genomes per group
# Rename with lowest classification available and IMG OID

# Read in MAG info file
mag_data <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Supplemental/MAG_information.csv", header = T)
mag_data$Taxonomy <- as.character(mag_data$Taxonomy)

# Make vector of last known classification
ID <- c()
for(i in 1:dim(mag_data)[1]){
  taxonomy <- strsplit(mag_data$Taxonomy[i], ";")
  ID[i] <- taxonomy[[1]][max(which(taxonomy[[1]] != "" & taxonomy[[1]] != "[Blank]"))]
}
ID[which(is.na(ID) == T)] <- "Unclassified"

mag_data$ID <- paste(ID, mag_data$IMG_OID, sep = "_")

phylum <- sapply(strsplit(mag_data$Taxonomy, ";"), '[', 1)
mag_data$Phylum <- phylum

# Read in IMG function data

enzymes <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/IMG_function_data.csv")
enzymes[,2] <- as.character(enzymes[,2])
colnames(enzymes) <- enzymes[1,]
rownames(enzymes) <- enzymes[,1]
enzymes <- enzymes[2:17895, 2:206]

# What phyla are in Lake Mendota and how many MAGs are in each phylum?

mendota <- mag_data[which(mag_data$Lake == "Mendota"), ]
table(mendota$Phylum)

# Start outputting data for input into the pathway analysis sheets

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "[Blank]")]
newnames <- mendota$ID[which(mendota$Phylum == "[Blank]")]
ME_unclassified <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_unclassified) <- newnames
write.csv(ME_unclassified, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_unclassified.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Actinobacteria")]
newnames <- mendota$ID[which(mendota$Phylum == "Actinobacteria")]
ME_actino <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_actino) <- newnames
write.csv(ME_actino, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_actino.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Bacteroidetes")]
newnames <- mendota$ID[which(mendota$Phylum == "Bacteroidetes")]
ME_bacteroid <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_bacteroid) <- newnames
write.csv(ME_bacteroid, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_bacteroidetes.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Chlamydiae" | mendota$Phylum == "Tenericutes" | mendota$Phylum == "Chloroflexi")]
newnames <- mendota$ID[which(mendota$Phylum == "Chlamydiae" | mendota$Phylum == "Tenericutes" | mendota$Phylum == "Chloroflexi")]
ME_other <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_other) <- newnames
write.csv(ME_other, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_other.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Cyanobacteria")]
newnames <- mendota$ID[which(mendota$Phylum == "Cyanobacteria")]
ME_cyano <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_cyano) <- newnames
write.csv(ME_cyano, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_cyano.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Planctomycetes")]
newnames <- mendota$ID[which(mendota$Phylum == "Planctomycetes")]
ME_plancto <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_plancto) <- newnames
write.csv(ME_plancto, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_plancto.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Proteobacteria")]
newnames <- mendota$ID[which(mendota$Phylum == "Proteobacteria")]
ME_proteo <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_proteo) <- newnames
write.csv(ME_proteo, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_proteo.csv")

IMGOID <- mendota$IMG_OID[which(mendota$Phylum == "Verrucomicrobia")]
newnames <- mendota$ID[which(mendota$Phylum == "Verrucomicrobia")]
ME_verruco <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(ME_verruco) <- newnames
write.csv(ME_verruco, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/ME_verruco.csv")

#TBE next
tbe <- mag_data[which(mag_data$Lake == "Trout Bog Epilimnion"), ]
table(tbe$Phylum)

IMGOID <- tbe$IMG_OID[which(tbe$Phylum == "Actinobacteria")]
newnames <- tbe$ID[which(tbe$Phylum == "Actinobacteria")]
TE_actino <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TE_actino) <- newnames
write.csv(TE_actino, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TE_actino.csv")

IMGOID <- tbe$IMG_OID[which(tbe$Phylum == "Proteobacteria")]
newnames <- tbe$ID[which(tbe$Phylum == "Proteobacteria")]
TE_proteo <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TE_proteo) <- newnames
write.csv(TE_proteo, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TE_proteo.csv")

IMGOID <- tbe$IMG_OID[which(tbe$Phylum == "Acidobacteria" | tbe$Phylum == "Bacteroidetes" | tbe$Phylum == "Chlorobi" | tbe$Phylum == "Verrucomicrobia")]
newnames <- tbe$ID[which(tbe$Phylum == "Acidobacteria" | tbe$Phylum == "Bacteroidetes" | tbe$Phylum == "Chlorobi" | tbe$Phylum == "Verrucomicrobia")]
TE_other <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TE_other) <- newnames
write.csv(TE_other, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TE_other.csv")

#TBH
tbh <- mag_data[which(mag_data$Lake == "Trout Bog Hypolimnion"), ]
table(tbh$Phylum)

IMGOID <- tbh$IMG_OID[which(tbh$Phylum == "Actinobacteria")]
newnames <- tbh$ID[which(tbh$Phylum == "Actinobacteria")]
TH_actino <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TH_actino) <- newnames
write.csv(TH_actino, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TH_actino.csv")

IMGOID <- tbh$IMG_OID[which(tbh$Phylum == "Bacteroidetes")]
newnames <- tbh$ID[which(tbh$Phylum == "Bacteroidetes")]
TH_bacteroid <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TH_bacteroid) <- newnames
write.csv(TH_bacteroid, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TH_bacteroid.csv")

IMGOID <- tbh$IMG_OID[which(tbh$Phylum == "Proteobacteria")]
newnames <- tbh$ID[which(tbh$Phylum == "Proteobacteria")]
TH_proteo <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TH_proteo) <- newnames
write.csv(TH_proteo, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TH_proteo.csv")

IMGOID <- tbh$IMG_OID[which(tbh$Phylum == "Verrucomicrobia")]
newnames <- tbh$ID[which(tbh$Phylum == "Verrucomicrobia")]
TH_verruco <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TH_verruco) <- newnames
write.csv(TH_verruco, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TH_verruco.csv")

IMGOID <- tbh$IMG_OID[which(tbh$Phylum == "Acidobacteria" | tbh$Phylum == "Chlorobi" | tbh$Phylum == "Elusimicrobia" | tbh$Phylum == "Ignavibacteria")]
newnames <- tbh$ID[which(tbh$Phylum == "Acidobacteria" | tbh$Phylum == "Chlorobi" | tbh$Phylum == "Elusimicrobia" | tbh$Phylum == "Ignavibacteria")]
TH_other <- enzymes[,match(IMGOID, colnames(enzymes))]
colnames(TH_other) <- newnames
write.csv(TH_other, "C:/Users/Alex/Desktop/MAGstravaganza/Pathway_analysis/TH_other.csv")


