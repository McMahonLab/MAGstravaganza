# MAGstravaganza
We propose that the name for a group of MAGs (metagenome-assembled genomes) is a MAGstravaganza, like a gaggle of geese or a congress of owls. Using time-series metagenomes and time-series resolved MAGs from three different freshwater environments, we analyze predicted functions and how these function co-occur in genomes. This reveals how carbon and nutrient cycles are connected by microbes in freshwater.

The publication associated with this GitHub repo is:
"Connections between freshwater carbon and nutrient cycles revealed through time series metagenomics." Linz, A.M., He, S., Stevens, S.L.R., Anatharaman, K. Rohwer, R.R., Malmstrom, R.R., Bertilsson, S., McMahon, K.D. Submitted. 2018.

Copyright (c) 2018, Katherine McMahon, Alex Linz, and friends

Repo Structure
------------------------------
					
- README.md: This file
- Code/
  - Tree_building/
  - Functional_marker_genes/
    - lefse_input.R: Code to take BLAST output table and format for input into LEfSe
    - process_lefse_output.R: Code to turn LEfSe output into a dataframe suitable for R plotting
  - Unused_analyses/
    - chtc-kraken.sh: Code to run kraken on metagenomes in high throughput (not included in manuscript)
    - kraken_report_processing.R: Code to take kraken output and make it into dataframes or plots
    - thesis_pie_charts.R: Alternative view of MAG diversity made for A. Linz's thesis
  - ANI_between_genomes.R: Make supplemental table of ANI from calculation output
  - MAGstravaganza_manuscript_plots.R: Code to generate all plots in the manuscript
  - amino_acid_bias.sh: Code to calculate amino acid bias from open reading frames


- Manuscript_drafts/					See the evolution of the manuscript!
- time_series_mapping/					Results and preliminary analysis of mapping the time series metagenomes to the MAGs
- dbCAN_results/					Output of annotating CAZy enzymes in the MAGs, one file per MAG
- Mansucript_plots/					Intermediate and file versions of plots in the mansucript
- Pathway_analysis					Data sheets and results used for calculating pathway presence/absence
- Supplemental/						Supplemental files accompanying the manuscript
