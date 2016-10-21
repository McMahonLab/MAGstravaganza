library(OTUtable)
#The metagenomes show that Trout Bog in 2008 behaved very differently than 05, 07, and 09.
#As far as I can tell, the phenology of hypolimnion development holds in 08
#By PCoA, community composition in 08 hypolimion is different than the other 3 years. However, when I ran significance tests, every year was significantly different than the other years.
#In fact, the 08 epilimion looks more similar to the hypolimnion than in other years.

#Check TB08 water column - I checked hypolimnion for anoxic conditions, but what if the epilimnion was unusually oxygen poor?

data(metadata)

TB08 <- make_do_matrix("TBE.....08", metadata)
plot_column(TB08, "Trout Bog 2008 Oxygen")

#Oxygen looks pretty stable all summer. Nothing remarkable here.
#Check temp

TB08 <- make_temp_matrix("TBE.....08", metadata)
plot_column(TB08, "Trout Bog 2008 Temperature")

#Also normal. Looks like we caught the spring mixing event, too.

#run indicator taxa analysis by year to see if there is anything notable absent in 2008 or only present in 2008

library(indicspecies)
data(otu_table)
data(taxonomy)

TBE <- bog_subset("TBE", otu_table)
TBH <- bog_subset("TBH", otu_table)

#Make tables at each phylogenetic level
Etribe_table <- combine_otus("Tribe", TBE, taxonomy)
Eclade_table <- combine_otus("Clade", TBE, taxonomy)
Elineage_table <- combine_otus("Lineage", TBE, taxonomy)
Eorder_table <- combine_otus("Order", TBE, taxonomy)
Eclass_table <- combine_otus("Class", TBE, taxonomy)
Ephylum_table <- combine_otus("Phylum", TBE, taxonomy)

#Change OTU number designations to full taxonomic assignment
named_otu_table <- TBE
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

#Combine tables at each level into one giant table
Efull_table <- rbind(named_otu_table, Etribe_table, Eclade_table, Elineage_table, Eorder_table, Eclass_table, Ephylum_table)

#Remove groups that are unclassified at any level
classified <- grep("unclassified", rownames(Efull_table))
classified1 <- grep("__$", rownames(Efull_table))
Einput_table <- Efull_table[-c(classified, classified1),]

#Keep only the top quantile in abundance
threshold <- quantile(rowSums(Einput_table))[4]
Einput_table <- Einput_table[which(rowSums(Einput_table) >= threshold),]

#Format for input into mulitplatt() function
Einput_table <- t(Einput_table)
Einput_table <- as.data.frame(Einput_table)

yearid <- c("05", "07", "08", "09")
years <- substr(rownames(Einput_table), start=9, stop=10)

Egroups <- c()
Egroups[which(years == "05")] <- 1
Egroups[which(years == "07")] <- 2
Egroups[which(years == "08")] <- 3
Egroups[which(years == "09")] <- 4

epi_years <- multipatt(x = Einput_table, cluster = Egroups, func = "r.g", control = how(nperm = 9999))

write.csv(epi_years$sign, file = "C:/Users/amlinz16/Desktop/TBEyears_indicators.csv")

#Make tables at each phylogenetic level
Htribe_table <- combine_otus("Tribe", TBH, taxonomy)
Hclade_table <- combine_otus("Clade", TBH, taxonomy)
Hlineage_table <- combine_otus("Lineage", TBH, taxonomy)
Horder_table <- combine_otus("Order", TBH, taxonomy)
Hclass_table <- combine_otus("Class", TBH, taxonomy)
Hphylum_table <- combine_otus("Phylum", TBH, taxonomy)

#Change OTU number designations to full taxonomic assignment
named_otu_table <- TBH
fullnames <- c()
for(i in 1:dim(taxonomy)[1]){
  fullnames[i] <- paste(taxonomy[i,], collapse = ";")
}
fullnames <- make.unique(fullnames)
rownames(named_otu_table) <- fullnames

#Combine tables at each level into one giant table
Hfull_table <- rbind(named_otu_table, Htribe_table, Hclade_table, Hlineage_table, Horder_table, Hclass_table, Hphylum_table)

#Remove groups that are unclassified at any level
classified <- grep("unclassified", rownames(Hfull_table))
classified1 <- grep("__$", rownames(Hfull_table))
Hinput_table <- Hfull_table[-c(classified, classified1),]

#Keep only the top quantile in abundance
threshold <- quantile(rowSums(Hinput_table))[4]
Hinput_table <- Hinput_table[which(rowSums(Hinput_table) >= threshold),]

#Format for input into mulitplatt() function
Hinput_table <- t(Hinput_table)
Hinput_table <- as.data.frame(Hinput_table)

yearid <- c("05", "07", "08", "09")
years <- substr(rownames(Hinput_table), start=9, stop=10)

Hgroups <- c()
Hgroups[which(years == "05")] <- 1
Hgroups[which(years == "07")] <- 2
Hgroups[which(years == "08")] <- 3
Hgroups[which(years == "09")] <- 4

hypo_years <- multipatt(x = Hinput_table, cluster = Hgroups, func = "r.g", control = how(nperm = 9999))

write.csv(hypo_years$sign, file = "C:/Users/amlinz16/Desktop/TBHyears_indicators.csv")

#Results
#Epi
#2008 does have odd anaerobic indicators - Desulfobulbus, Chlorobi, and Chloroflexi. All anaerobic sulfur reducers that are not normally found in the epilimnion.
#Methylophilaceae and Methylococcaceae are indicators of other years, esp 2007
#Things NOT present in 2008:
#acV, Geobacter, and bacI

#Hypo
#Lots of indicators here. Bunch of Pnec, OP3. Syntrophorhabdaceae, Helicobacter sulfurimonas, LD28, gamIV, Caulobacter, Luna1-A, Mycobacterium, Deinococcus, Sphingomonas, Chlorobi, gamIII, Clostridia, Pedosphaera, Flavobacteriia, Rhodocyclaceae
#NOT present
#Alphaproteobacteria, betI, Plancotmyces;phycisphaera, Chlamydiae

#The presence of sulfur reducers as abundant indicators of the epilimnion is odd. Check their abundance in Trout Bog manually.

TBE.dates <- extract_date(colnames(TBE))
chlorobi <- grab_group("Chlorobi", "Phylum", TBE, taxonomy)
plot(TBE.dates, colSums(chlorobi), pch=16)
#Chlorobi holds true. One bloom in 2007, some abundance in 2009, but predominantly in 2008.

chloroflexi <- grab_group("Chloroflexi", "Phylum", TBE, taxonomy)
plot(TBE.dates, colSums(chloroflexi), pch=16)
#Not enough abundance to say anything.

desulfo <- grab_group("Desulfobulbus", "Clade", TBE, taxonomy)
plot(TBE.dates, colSums(desulfo), pch=16)
#Wow, really striking here. Similar trend of bloom at the end of 07, abundant in 08, low abundnce in 09.

#The big question - WHY are there sulfur reducers in the epilimnion in 2008?
#May have even started fall of 2007 and carried into 2009.
#If sulfur reducers are abundant, then SO4 must be abundant, too.
#What would make SO4 abundant? Oxidation of H2S (leading to a cycle where SO4 builds up until reducers take over? Would this carry over through mixing events?)
#Do the genes tell the same story?

epi.data <- read.csv("C:/Users/Alex/Dropbox/Methylotroph_analysis/epi_metagenome_genesurvey.csv", row.names=1)
library(ggplot2)
genes <- c("sulfate adenylyltransferase", "sulfite reductase", "methane monooxygenase", "formaldehyde activating", "dissimilatory adenylylsulfate reductase", "sulfide dehydrogenase")

find <- grep(genes[1], rownames(epi.data))
sulfate <- epi.data[find,]
total.sulfate <- colSums(sulfate)/colSums(epi.data)

find <- grep(genes[2], rownames(epi.data))
sulfite <- epi.data[find,]
total.sulfite <- colSums(sulfite)/colSums(epi.data)

find <- grep(genes[3], rownames(epi.data))
meth <- epi.data[find,]
total.meth <- colSums(meth)/colSums(epi.data)

find <- grep(genes[4], rownames(epi.data))
form <- epi.data[find,]
total.form <- colSums(form)/colSums(epi.data)

find <- grep(genes[5], rownames(epi.data))
diss <- epi.data[find,]
total.diss <- colSums(diss)/colSums(epi.data)

find <- grep(genes[6], rownames(epi.data))
sulfide <- epi.data[find,]
total.sulfide <- colSums(sulfide)/colSums(epi.data)

totals <- c(total.sulfate, total.sulfite, total.meth, total.form, total.diss, total.sulfide)
dates <- rep(as.Date(substr(names(total.sulfate), start=2, stop=12), format = "%Y.%m.%d"), 6)
key <- rep(genes, each=45)

data <- data.frame(key, dates, totals)

ggplot(data, aes(x=dates, y=totals, color=key)) + geom_line(size=1.5)

#No change in sulfur related gene abundance in the epilimion that I can see. I'll try hypolimnion next.

hypo.data <- read.csv("C:/Users/Alex/Dropbox/Methylotroph_analysis/hypo_metagenome_genesurvey.csv", row.names=1)

find <- grep(genes[1], rownames(hypo.data))
sulfate <- hypo.data[find,]
total.sulfate <- colSums(sulfate)/colSums(hypo.data)

find <- grep(genes[2], rownames(hypo.data))
sulfite <- hypo.data[find,]
total.sulfite <- colSums(sulfite)/colSums(hypo.data)

find <- grep(genes[3], rownames(hypo.data))
meth <- hypo.data[find,]
total.meth <- colSums(meth)/colSums(hypo.data)

find <- grep(genes[4], rownames(hypo.data))
form <- hypo.data[find,]
total.form <- colSums(form)/colSums(hypo.data)

find <- grep(genes[5], rownames(hypo.data))
diss <- hypo.data[find,]
total.diss <- colSums(diss)/colSums(hypo.data)

find <- grep(genes[6], rownames(hypo.data))
sulfide <- hypo.data[find,]
total.sulfide <- colSums(sulfide)/colSums(hypo.data)

totals <- c(total.sulfate, total.sulfite, total.meth, total.form, total.diss, total.sulfide)
dates <- rep(as.Date(substr(names(total.sulfate), start=2, stop=12), format = "%Y.%m.%d"), 6)
key <- rep(genes, each=45)

data <- data.frame(key, dates, totals)

ggplot(data, aes(x=dates, y=totals, color=key)) + geom_line(size=1.5)

#Not necessarily a higher abundance of any one sulfur genes. There is that same odd matching in trends in 2008 between sulfate adenylyltransferase, sulfide dehydrogenase, and sulfite reductase.

#Maybe sulfate was at a high enough concentration to no longer be limiting, allowing sulfur reducers to outcompete in other ways?

#Did anyone measure sulfur levels in Trout Bog during that time period?

#Alternative hypothesis: the chlorobi sweep occurred in 2008. We know chlorobi can have massive layers right at the interface of the epilimnion and hypolimnion. this could be consuming oxygen and producing C1 compounds such as formate, formaldehyde, methanol, and acetate, leading to the consistency of gene and contig correlations in 2008. Can our chlorobi MAGs produced this?
#No production genes (but I'm not sure what these would look like.) Not a single consumption gene though. It's entirely possible that C1 compounds are being produced in Chlorobi that it cannot use, so it spews them out for Methyloversatilis to use.
#Methylococcales, whose superpower is methane to C1, therefore does not compete as well.

#this is a beautiful story. What do I need to prove it?
#Get readcoverage of Chlorobi - is there truly a freakishly high abundance in the epilimnion in the MAGs?
#Does it make sense for Chlorobi to be hanging out at the interface?
#Can Chlorobi isolates produce C1 compounds? It is known to hang out in consortia, possibly with sulfate-reducing Deltaproteobacteria (Desulfobulbus?)
#Are other taxa correlated to Chlorobi in 2008 but not in other years?
#Bonus: what is it about this clone that allows it to outcompete other groups and grow unchecked (acquistion of CBB pathway?)
#Likely not CBB, the C. clathratiformes type strain has this. However H111 (the sweep strain) is the only one that has it.
#If I align 111 and 211, I can see if there are any major differences.
#The Switzerland C. clathratiformes clonal populations uses reductive TCA to fix carbon
#No RuBisCos anywhere in IMG chlorobis or in the literature that I can see
#One other MAG uses reductive TCA as well.
