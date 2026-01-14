# Standard workflow for single-cell RNA-Seq analysis

# Load libraries
library(Seurat)
library(tidyverse)

# 1. Load Data ---------------------------------------
# Replace the path below with your actual file path
sparse_matrix <- Read10X_h5(filename = '../data/sample_feature_bc_matrix.h5')
raw_counts <- sparse_matrix$`Gene Expression`

# Initialize Seurat object
sc_obj <- CreateSeuratObject(counts = raw_counts, 
                             project = "SC_Project", 
                             min.cells = 3, 
                             min.features = 200)

# 2. Quality Control ---------------------------------
# Calculate mitochondrial percentage
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')

# 3. Filtering ---------------------------------------
# Adjust these thresholds based on your QC plots above
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & 
                 nFeature_RNA < 2500 & 
                 percent.mt < 5)

# 4. Normalization & Feature Selection ---------------
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

# Identify and plot top 10 variable genes
top10_genes <- head(VariableFeatures(sc_obj), 10)
plot1 <- VariableFeaturePlot(sc_obj)
LabelPoints(plot = plot1, points = top10_genes, repel = TRUE)

# 5. Scaling & PCA -----------------------------------
all_genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all_genes)
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))

# Determine dimensionality
ElbowPlot(sc_obj)

# 6. Clustering --------------------------------------
# Using 15 dimensions based on original script; adjust after checking ElbowPlot
sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# Set default identity to a specific resolution
Idents(sc_obj) <- "RNA_snn_res.0.5"

# 7. UMAP Visualization ------------------------------
sc_obj <- RunUMAP(sc_obj, dims = 1:15)
DimPlot(sc_obj, reduction = "umap", label = TRUE)
