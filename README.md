# A single cell atlas of human teeth
Healthy periodontium and pulp single-cell data

# Citation

The code in this repository pertains to publication:

Pierfrancesco Pagella, Laura de Vargas Roditi, Bernd Stadlinger, Andreas E. Moor, Thimios A. Mitsiadis. A single cell atlas of human teeth.

doi: https://doi.org/10.1101/2021.02.19.431962 

# Data

Data pertaining the above publication can be found at GEO. Access number: GSE161267

# Code description: 
R:

20191129_healthy_pulp.R       : Pre-processing of single-cell data from pulp samples only

20191204_perio.R              : Pre-processing of single-cell data from periodontium samples only

20191204_merged_perio_pulp.R  : Pre-processing of single-cell data from pulp and periodontium samples together

Dis_graph.R                   : Extended Jaccard similarity between pulp and perio sampels plotted as a force directed graph layout 

Dis_graph_review.R            : Same as above but pairs of smallest distances were ranked from smallest to largest for visualization purposes

DoMultiBarHeatmap.R           : Code to allow multiple identities to be visualized as a bar on for of an ordered heatmap

Mean_median_table_pulp_perio.R: Statistics on pulp and periodontium cell population sizes and genes per cell 

boxplot_byMCT.R               : Boxplot of major cell types in each pulp and periodontium dataset

diff_exp.R                    : Differential gene expression

umap_heatmaps.R               : UMAPs and heatmaps of pulp and periodontium datasets

velocity_perio_step1_R.R      : Data prep for velocity estimate of periodontium samples with scvelo in python

velocity_pulp_step1_R.R       : Data prep for velocity estimate of pulp samples with scvelo in python

umaps_supp.R                  : Supplementary umap plots

perio_review.R and pulp_review.R were run to adjust the original labels in 20191129_healthy_pulp.R and 20191204_perio.R after re-naming of original file names to Perio 1-5 and Pulp 1-5 according to GEO submission

Python:

velocity_perio_step2_scvelo.py              : Velocity estimate for periodontium samples (part 1 done in R)

velocity_pulp_step2_scvelo.py               : Velocity estimate for pulp samples (part 1 done in R)
