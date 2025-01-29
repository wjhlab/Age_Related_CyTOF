# Overview

This repository includes all relevant custom R scripts and associated configuration files for the suspension mass cytometry analysis of “Age-related divergence of circulating immune responses in patients with solid tumors treated with immune checkpoint inhibitors”

# Input Data

All debarcoded fcs files and a fully annotated data frame ("backup_output.rds") are available at 10.5281/zenodo.14755936. 

The fully annotated data frame can be loaded onto the R script to generate the published figures in the manuscript. In Config, there are metadata, panel, and merged (annotation) files necessary to generate the heatmap and other plots for the entire dataset. Fcs files will also need to be in a folder named "Data" in working directory as well. A UMAP dataframe is also available at Zenodo ("backup_UMAP.rds").

For PCA and correlation plots, excel files containing raw data necessary to generate these plots are provided.

# R Scripts

Scripts stored in Rscript format.
