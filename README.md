# MAGstravaganza
Investigation of how freshwater nutrient cycling is linked at the organism level.

The goal of this analysis is to determine pathway presence/absence in genomes assembled from metagenomes (MAGs) in two lakes with very different environmental conditions.
We hypothesize that the pathways present in genomes will reflect the functioning of the ecosystem where those genomes belong.

Copyright (c) 2016, Katherine McMahon, Alex Linz, and friends

Repo Structure
------------------------------
/MAGstravagansa						A party of MAGs
- README.md						This file
- paper_outline.docx					The current state of the manuscript on this project. Up for review at lab meeting on Oct 6.
- .RData, .Rhistory					Rstudio log files
- /Plots
	- MEabundance.jpeg				This file and those below are barplots of read phylogenetic assignments in the pooled samples for each lake.
	- MEabundance.pdf
	- TBabundance.pdf
	- TBEabundance.jpeg
	- TBHabundnace.jpeg
	- /ISME_poster
		- Linz_ISME16.ai			Pre-lab review
		- Linz_ISME16.pdf
		- Linz_ISME16_v2.ai			Version presented at ISME
		- Linz_ISME16.v2.pdf
		- poster_citations.docx			Citations for final poster
		- poster_feedback.docx			Changes from v1 to v2

- /Data_files
	- metacyc_pathway_tree.csv			Shows all broader categories of metacyc functions. I made this by hand. You're welcome.
	- /Genome_info
		- MAG_metadata.csv			Info about genome size and completeness
		- revised_MAG_metadata.csv		Old genomes removed
	- /Metapathways_output
		- combined_MAGs_pathwaycoverage.txt	Coverage information about each pathway in each genome
		- combined_MAGs_pws_rxns.orfs		Location and reactions of each pathway in each genome
		- curated_rxns.csv			Only pathways with 75% coverage remain
		- master_table.csv			Wide format of curated_rxns.csv
		- master_table.xls			Additional sheets of exploratory analysis
		- rxn_key				Reaction names of PWY numbers in output 	
	-/Phylogeny_of_metagenomes			IMG results of phylogeny BLAST on samples pooled by lake
		- ME_IMGblast_phylogeny.csv
		- TBepi_IMGblast_phylogeny.csv
		- TBhypo_IMGblast_phylogeny
- /Processing_code
	- MAGstravaganza_workflow_and_notes.txt		Commands for Metapathways/Pathways Tools processing
	- 01parse_rxn_table.R				Makes curated_rxns.csv (pathways < 75% coverage removed)
	- 02master_table_script.R			Converts to wide format in master_table.csv
- /TB_vs_ME_analysis
	- combined_assembly_read_abundances.R		Plots IMG phylogeny BLAST data
	- /Tree_building
		- 00MoveFastaFiles.sh			Script to move fasta files out of zipped IMG folders
		- 01FastaSequenceMerger.pl		From J.Hamilton at McMahonLab/Scripts repo. Combines contigs into a single sequence
		- 02RunPhylosift.pl			From J.Hamilton at McMahonLab/Scripts repo. Runs search and alignment steps of Phylosift on multiple sequences
		- 03CreateAlignmentFile.pl		From J.Hamilton at McMahonLab/Scripts repo. Creates a single alignment file of all genomes for input into FastTree
		- ME_MAGS.nwk				Output of FastTree on Mendota MAGS
		- TBE_MAGS.nwk				Output of FastTree on Trout Bog Epi MAGS
		- TBH_MAGS.nwk				Output of FastTree on Trout Bog Hypo MAGS
		- test_tree.R				First pass at code for plotting the trees with pathway presence/absence indicators
- /OLDFunction_correlation_analysis
	My original plan was to correlated presence/absence of function on genomes.
	However, we don't have enough datapoints for this.
	I'm mostly disregarding this folder, but saving it in case I want to recycle something.