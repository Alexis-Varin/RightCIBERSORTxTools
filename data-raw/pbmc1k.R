library(Seurat)
library(cowplot)
ctrl.data <- read.table(file = "../data/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "../data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

# Integrate objects
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# UMAP and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

# Downsampling to about 1000 cells
pbmc1k = WhichCells(immune.combined, downsample = 75)
pbmc1k = immune.combined[,pbmc1k]

# Rename clusters of seurat_clusters
clusters.names = paste(rep("Cluster", length(levels(pbmc1k))),levels(pbmc1k),sep=".")
names(clusters.names) = levels(pbmc1k)
pbmc1k = RenameIdents(pbmc1k, clusters.names)
pbmc1k@meta.data$seurat_clusters = pbmc1k@active.ident
DimPlot(pbmc1k, reduction = "umap", label = TRUE)
Idents(pbmc1k) = "stim"
DimPlot(pbmc1k, reduction = "umap", label = TRUE)
Idents(pbmc1k) = "seurat_clusters"
DimPlot(pbmc1k, reduction = "umap", label = TRUE)

# Seurat object ready to be used for CIBERSORTx_Reference_Matrix_Builder examples
Idents(pbmc1k) = "stim"
usethis::use_data(pbmc1k)
