# Doublet-finder.R
# Purpose: filter out doublets using DoubletFinder

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

# 1. Load Data -------------------------------------------------------------
# Provide paths to your matrix, features, and barcodes files
cts <- ReadMtx(mtx = 'data/matrix.mtx.gz',
               features = 'data/features.tsv.gz',
               cells = 'data/barcodes.tsv.gz')

# Initialize Seurat object
sc_obj <- CreateSeuratObject(counts = cts)

# 2. QC and Filtering ------------------------------------------------------
# Calculate mitochondrial percentage
sc_obj$mitoPercent <- PercentageFeatureSet(sc_obj, pattern = '^MT-')

# Filter cells based on count, feature, and mitochondrial thresholds
sc_obj_filtered <- subset(sc_obj, subset = nCount_RNA > 800 &
                            nFeature_RNA > 500 &
                            mitoPercent < 10)

# 3. Standard Pre-processing -----------------------------------------------
# These steps are required before running DoubletFinder
sc_obj_filtered <- NormalizeData(object = sc_obj_filtered)
sc_obj_filtered <- FindVariableFeatures(object = sc_obj_filtered)
sc_obj_filtered <- ScaleData(object = sc_obj_filtered)
sc_obj_filtered <- RunPCA(object = sc_obj_filtered)
sc_obj_filtered <- FindNeighbors(object = sc_obj_filtered, dims = 1:20)
sc_obj_filtered <- FindClusters(object = sc_obj_filtered)
sc_obj_filtered <- RunUMAP(object = sc_obj_filtered, dims = 1:20)

# 4. pK Identification -----------------------------------------------------
# Find the optimal pK neighborhood size
sweep.res.list <- paramSweep_v3(sc_obj_filtered, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Select pK with the maximum BCmetric
pK_value <- bcmvn %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK_value <- as.numeric(as.character(pK_value[[1]]))

# 5. Doublet Proportion Estimation -----------------------------------------
annotations <- sc_obj_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# Estimate doublet rate (e.g., 7.5% for ~10k cells) and adjust for homotypic doublets
nExp_poi <- round(0.076 * nrow(sc_obj_filtered@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# 6. Run DoubletFinder -----------------------------------------------------
sc_obj_filtered <- doubletFinder_v3(sc_obj_filtered, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK_value, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)

# 7. Visualization and Cleanup ---------------------------------------------
# Identify the dynamic column name created by DoubletFinder
df_col <- grep("DF.classifications", colnames(sc_obj_filtered@meta.data), value = TRUE)

# Visualize results
DimPlot(sc_obj_filtered, reduction = 'umap', group.by = df_col)

# Print summary table
table(sc_obj_filtered@meta.data[[df_col]])

# Optional: Filter for singlets only
# sc_obj_singlets <- subset(sc_obj_filtered, cells = WhichCells(sc_obj_filtered, expression = !!sym(df_col) == "Singlet"))
