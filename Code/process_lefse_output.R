# Process LEfSe results

lefse <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/Galaxy36-[B)_LDA_Effect_Size_(LEfSe)_on_data_35].lefse_internal_res", fill = TRUE)
genedata <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/metabolic_gene_info.txt", fill = TRUE)
lefse_input <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/alllefse_input.txt", row.names = 1)

# Make supplemental table of gene, function, and if/where it was significantly more abundant

supp_markers <- lefse[,c(1,3)]
genedata$V1 <- gsub("-", "_", genedata$V1)
supp_markers$annotation <- genedata$V3[match(as.character(supp_markers$V1), as.character(genedata$V1))]
#Fill in ones that got screwed up in all the program transfers
supp_markers$annotation[which(supp_markers$V1 == "gi.118567952.gb.ABL02757_1.sqr")] <- "sulfide_quinone_reductase"
supp_markers$annotation[which(supp_markers$V1 == "gi.500800570.ref.WP_011998629_1.4_hydroxybutyryl_CoA_dehydratase")] <- "alt_carbon_fixation"
supp_markers$annotation[which(supp_markers$V1 == "sor_PF07682_full_stockholm_sample1")] <- "sulfur_oxygenase_reductase"
supp_markers$annotation[which(supp_markers$V1 == "gi.118567952.gb")] <- "sulfide_quinone_reductase"
supp_markers$annotation[which(supp_markers$V1 == " f_4_hydroxybutyryl_CoA_synthetase")] <- "alt_carbon_fixation"
supp_markers$annotation[which(supp_markers$V1 == "NrfA_PF02335_full__stockholm_sample1")] <- "cytochrome_c_nitrite_reductase"
supp_markers$annotation[which(supp_markers$V1 == "NrfA_PF02335_full__stockholm_sample1")] <- "cytochrome_c_nitrite_reductase"
supp_markers$annotation[which(supp_markers$V1 == "MmoB_PF02406_full_stockholm_sample1")] <- "methane_monooxygenase"
supp_markers$annotation[which(supp_markers$V1 == "NapB_PF03892_full_stockholm_sample1")] <- "nitrate_reductase_cytochrome_c-type"

supp_markers <- supp_markers[which(is.na(supp_markers$annotation) == F), ]

#Sort by significance
supp_markers <- supp_markers[order(supp_markers$V3), ]

# Output as supplemental doc
write.csv(supp_markers, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Supplemental/marker_genes_LDA_significance.csv", row.names = F)


# Do the same with epi/hypo TB file

lefse <- read.table("C:/Users/Goose and Gander/Desktop/MAGstravaganza/Data_files/Galaxy44-[B)_LDA_Effect_Size_(LEfSe)_on_data_43].csv", fill = TRUE)


# Make supplemental table of gene, function, and if/where it was significantly more abundant

supp_markers <- lefse[,c(1,3)]
genedata$V1 <- gsub("-", "_", genedata$V1)
supp_markers$annotation <- genedata$V3[match(as.character(supp_markers$V1), as.character(genedata$V1))]
#Fill in ones that got screwed up in all the program transfers
supp_markers$annotation[which(supp_markers$V1 == "gi.118567952.gb.ABL02757_1.sqr")] <- "sulfide_quinone_reductase"
supp_markers$annotation[which(supp_markers$V1 == "gi.500800570.ref.WP_011998629_1.4_hydroxybutyryl_CoA_dehydratase")] <- "alt_carbon_fixation"
supp_markers$annotation[which(supp_markers$V1 == "sor_PF07682_full_stockholm_sample1")] <- "sulfur_oxygenase_reductase"
supp_markers$annotation[which(supp_markers$V1 == "gi.118567952.gb")] <- "sulfide_quinone_reductase"
supp_markers$annotation[which(supp_markers$V1 == " f_4_hydroxybutyryl_CoA_synthetase")] <- "alt_carbon_fixation"
supp_markers$annotation[which(supp_markers$V1 == "NrfA_PF02335_full__stockholm_sample1")] <- "cytochrome_c_nitrite_reductase"
supp_markers$annotation[which(supp_markers$V1 == "NrfA_PF02335_full__stockholm_sample1")] <- "cytochrome_c_nitrite_reductase"
supp_markers$annotation[which(supp_markers$V1 == "MmoB_PF02406_full_stockholm_sample1")] <- "methane_monooxygenase"
supp_markers$annotation[which(supp_markers$V1 == "NapB_PF03892_full_stockholm_sample1")] <- "nitrate_reductase_cytochrome_c-type"

supp_markers <- supp_markers[which(is.na(supp_markers$annotation) == F), ]

#Sort by significance
supp_markers <- supp_markers[order(supp_markers$V3), ]

# Output as supplemental doc
write.csv(supp_markers, file = "C:/Users/Goose and Gander/Desktop/MAGstravaganza/Supplemental/TBmarker_genes_LDA_significance.csv", row.names = F)
