library(ggplot2)
library(cowplot)

trait_info <- read.csv("C:/Users/Alex/Desktop/MAGstravaganza/Data_files/Genome_info/MAG_phylogeny.csv", header = T)
trait_info$Phylum <- sapply(strsplit(as.character(trait_info$Taxonomy.README.),";"), `[`, 1)
trait_info$Lake <- substr(trait_info$GenomeName, start = 1, stop = 2)


abun_trait_info <- trait_info[which(trait_info$Phylum == "Actinobacteria" | trait_info$Phylum == "Bacteroidetes" | trait_info$Phylum == "Planctomycetes" | trait_info$Phylum == "Verrucomicrobia" | trait_info$Phylum == "Acidobacteria" | trait_info$Phylum == "Proteobacteria" | trait_info$Phylum == "Cyanobacteria" ), ]

# I've already looked at genome size vs phylum and it has trends.
# Look at codon usage by lake and by phylum

ggplot(data = abun_trait_info, aes(x = CodonUsage, fill = Lake)) + geom_density(alpha = 0.5)
#Not a huge difference. By phylum:

ggplot(data = abun_trait_info, aes(x = CodonUsage, fill = Phylum)) + geom_density(alpha = 0.5)
# Now we're talking! Split into two plots to make it easier to see

ggplot(data = abun_trait_info[which(abun_trait_info$Phylum == "Actinobacteria" | abun_trait_info$Phylum == "Bacteroidetes" | abun_trait_info$Phylum == "Planctomycetes" | abun_trait_info$Phylum == "Verrucomicrobia"), ], aes(x = CodonUsage, fill = Phylum)) + geom_density(alpha = 0.5)
ggplot(data = abun_trait_info[which(abun_trait_info$Phylum == "Acidobacteria" | abun_trait_info$Phylum == "Proteobacteria" | abun_trait_info$Phylum == "Cyanobacteria" ), ], aes(x = CodonUsage, fill = Phylum)) + geom_density(alpha = 0.5)

# Most are pretty spread out, but Actinos are predominantly low codon bias (slow growing) and Acidos are fast

#Scatterplot of genome size vs codon usage, colored by various things

ggplot(data = abun_trait_info, aes(x = Est_Genome_Size, y = CodonUsage, color = Phylum)) + geom_point(size = 3, alpha = 0.5) + theme_bw() + scale_color_brewer(palette = "Set2")
ggplot(data = abun_trait_info, aes(x = Est_Genome_Size, y = CodonUsage, color = Motility)) + geom_point(size = 3, alpha = 0.5) + theme_bw()
ggplot(data = abun_trait_info, aes(x = Est_Genome_Size, y = CodonUsage, color = Carbon_Source)) + geom_point(size = 3, alpha = 0.5) + theme_bw()
ggplot(data = abun_trait_info, aes(x = Est_Genome_Size, y = CodonUsage, color = Num_C_MetaP)) + geom_point(size = 3, alpha = 0.5) + theme_bw() + scale_color_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 20)
hist(abun_trait_info$Num_C_MetaP)


#Things that I can work with:
#large genome size -> more likely to be motile
#Genome size and codon usage vary by phylum
#Carbon source and number of pathways don't really help much
#Codon usage bias and genome size not correlated
