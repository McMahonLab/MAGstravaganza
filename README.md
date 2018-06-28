# Connections between freshwater carbon and nutrient cycles revealed through reconstructed community genomes

We propose that the name for a group of MAGs (metagenome-assembled genomes) is a MAGstravaganza, like a gaggle of geese or a congress of owls. Using time-series metagenomes and time-series resolved MAGs from three different freshwater environments, we analyze predicted functions and how these function co-occur in genomes. This reveals how carbon and nutrient cycles are connected by microbes in freshwater.

The publication associated with this GitHub repo is:
"Connections between freshwater carbon and nutrient cycles revealed through time series metagenomics." Linz, A.M., He, S., Stevens, S.L.R., Anantharaman, K. Rohwer, R.R., Malmstrom, R.R., Bertilsson, S., McMahon, K.D. Submitted. 2018.

Copyright (c) 2018, Katherine McMahon, Alex Linz, and friends

#### Are you looking to use data or code from our paper?
You're in the right place. Raw data is located in IMG (https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) or on our lab fileshare.
- Use "IMG_Genome_ID" from the file "MAGs_Manuscript_2018-06-05/supplemental_documents/Supplemental_Data_S1.csv" to download the metagenomic time series.
- Use "IMG_Genome_ID" from the file "MAGs_Manuscript_2018-06-05/supplemental_documents/Supplemental_Data_S1.csv" to download MAGs.
- Use "IMG_Genome_ID" from the file "MAGs_Manuscript_2018-06-05/supplemental_documents/Supplemental_Table_S2.csv" to download pooled metagenome assemblies.
- If you're looking for intermediate files or code, check the list of files in this repo below.

Not seeing what you're looking for? Feel free to email us at amlinz@wisc.edu or trina.mcmahon@wisc.edu.


Repo Structure
------------------------------
					
- README.md: This file
- Code/
  - Tree_building/: Scripts to build a rough tree from MAGs. Full instructions available at https://github.com/McMahonLab/Scripts/tree/master/Phylogeny
  - Functional_marker_genes/
    - lefse_input.R: Code to take BLAST output table and format for input into LEfSe
    - process_lefse_output.R: Code to turn LEfSe output into a dataframe suitable for R plotting
  - Unused_analyses/
    - chtc-kraken.sh: Code to run kraken on metagenomes in high throughput (not included in manuscript)
    - kraken_report_processing.R: Code to take kraken output and make it into dataframes or plots
    - thesis_pie_charts.R: Alternative view of MAG diversity made for A. Linz's thesis
    - genome_size_figure.R: Initial testing of genome size vs. number of proposed substrates
  - ANI_between_genomes.R: Make supplemental table of ANI from calculation output
  - MAGstravaganza_manuscript_plots.R: Code to generate all plots in the manuscript
  - amino_acid_bias.sh: Code to calculate amino acid bias from open reading frames
  - process_tag_data.R: process .rData files that were output by TaxAss 16S classification
- Manuscript_drafts/
  - See the evolution of the manuscript!
- time_series_mapping/
  - Data files with metagenome information (ie dates, sizes, and library codes)
  - Lake_results.txt files have the number of reads mapping to each MAG from each metagenome, and the size of the metagenome
- dbCAN_results/
  - Output of annotating CAZy enzymes in the MAGs, one file per MAG
- Mansucript_plots/
  - Intermediate files of plots in the mansucript
- Pathway_analysis
  - IMG_function_data.csv: annotations downloaded from IMG
  - consolidated_pathway_data.csv/xls: aggregated pathway prediction data for input into Figure 3
  - *.xls files with lake, phylum, and "pathways in the name are the pathwya prediction files
  - *.csv files with lake and phylum have pathway data aggregated by phylum
  - split_by_lake_phylum.R: split the IMG_function_data by lake and phylum
  - pathways_template.xls: Template for running pathway prediction. Paste IMG function data into the 1st worksheet to use.
- Supplemental/
  - Assorted files used for supplemental tables, figures, and datasets. Mostly redundant with supplemental info accompanying the manuscript.
