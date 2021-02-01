# organoid-paper
This repository contains code for scRNAseq and Nanostring analysis for the following manuscript : "Transplanted human organoids empower PK/PD assessment of drug candidate for the clinic".

All analyses were run within docker containers.

Description of the folders/files:

**Nanostring/heatmap.ipynb**
Hierarchical clustering of and visualization of Nanostring gene expression on heatmaps (Sup Figure 1F and Sup Figure 5)

**Nanostring/data/**
selected podocyte genes, Nanostring gene list

**scRNAseq/R/**
helper functions for plot generation, saving etc...

**scRNAseq/data/**
sample IDs list, cluster annotations

**scRNAseq/analysis/heatmap_extract_clusteraverages.R**
extract averages for Nanostring genes from all single cell in vitro clusters

**scRNAseq/analysis/heatmap_plotting.R**
plotting cluster averages of Nanostring genes in heatmaps (Sup Figure 1F and Sup Figure 5)

**scRNAseq/analysis/dirichlet.Rmd**
Dirichelet regression for comparing cell proportions between days (Figure 1C).

**scRNAseq/analysis/figS7.Rmd**
Plots for Sup Figure 7

**scRNAseq/analysis/figS8.Rmd**
Plots for Sup Figure 8

**scRNAseq/analysis/human_qc.Rmd**
QC and filtering for publically sources human fetal and adult kidney scRNAseq data (Figure 1D and Sup Figure 7 and 8)

**scRNAseq/analysis/in_vitro_analysis**
QC and filtering of raw matrices for all in vitro experiments

**scRNAseq/analysis/in_vitro_cell_types.Rmd**
in vitro organoid stacked bargraphs (Sup Figure 1D) and Jenson-Shannon Divergences (Sup Figure 3)

**scRNAseq/analysis/in_vitro_UMAPs.Rmd**
UMAP plots for in vitro organoids (Sup Figure 2)

**scRNAseq/analysis/in_vivo_cell_types.Rmd**
in vivo organoid stacked bargraphs (Sup Figure 1B) and Jenson-Shannon Divergences (Sup Figure 6)

**scRNAseq/analysis/in_vivo_UMAPs.Rmd**
UMAP plots for in vivo organoids (Sup Figure 2)

**scRNAseq/analysis/pseudotime_analysis**
Pseudotime analysis and plot generation for in vivo podocytes (Figure 1D and Sup Figure 7) and collecting duct like cells (Sup Figure 8).

**scRNAseq/analysis/transplant_integration**
QC and filtering of raw matrices for all in vivo experiments

