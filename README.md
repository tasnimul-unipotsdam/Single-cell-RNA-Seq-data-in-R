# Single-cell RNA-Seq Analysis & Visualization in R
This repository contains a comprehensive pipeline for processing single-cell RNA sequencing (scRNA-seq) data using Seurat and visualizing gene expression profiles across diverse conditions using Tidyverse.

## 1. Single-Cell Pipeline
The singleCell_standard_workflow.R script performs:

Quality Control: Filtering cells based on mitochondrial percentage (percent.mt) and unique feature counts.

Normalization: Global-scaling normalization (LogNormalize).

Dimensionality Reduction: PCA and UMAP visualization.

Clustering: Graph-based clustering at multiple resolutions (0.1 to 1.0).

Doublet Detection & Quality Control: It ensures that clusters represent single cells rather than "multiplets" or "doublets" (two or more cells captured in one droplet).

## 2. Visualization & Statistics
The visualize-gene-expression.R script provides publication-ready plots:

Violin & Boxplots: Includes automated p-value brackets and significance stars using ggpubr.

Density Plots: Tissue-specific expression distributions.

Heatmaps: Multi-gene expression patterns across samples.

Scatter Plots: Gene-gene correlation analysis with linear regression.

### Prerequisites

Ensure you have the following R packages installed:

```R
install.packages(c("tidyverse", "Seurat", "ggpubr", "ggplot2", "hdf5r"))
