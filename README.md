# cellstruct

cellstruct provides three distinct metric scores to quantify the preservation of biological information between any two embeddings.
* global single-cell (GS) - Pearson's correlation of two cell-cell distance matrices from reference (e.g. PCA) and reduced (e.g. UMAP) embeddings using 10,000 waypoint cells.
* local single-cell (LS) - ratio of average squared distances of reference and reduced neighbors from the target cell, in reference embedding.
* global cluster (GC) - Pearson's correlation of two centroid-centroid distance matrices from reference and reduced embeddings (centroid refers to cluster centroid).

cellstruct is available as an R package. It is developed in Ouyang Lab by Jui Wan Loh and John Ouyang.

## Installation
The following dependencies are required to install and run cellstruct:
* reticulate
* Matrix
* Seurat
* ComplexHeatmap 
```R
install.packages("reticulate") #do it for Matrix and Seurat too
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("the-ouyang-lab/cellstruct", build_vignettes=TRUE)
```

## Usage
cellstruct requires at least two embeddings for comparison. GC scores will be biologically meaningful, if cell type annotation is provided. 
By default, GC and GS metrics are calculated. LS scores will be computed if `plot_LS` is set to true. An example of running script:
```R
library(cellstruct)
ref.proj <- "pca"
target.proj <- c("umap","tsne","fdl")
clusterVar <- "cell.type.fine"
seuName <- paste0("~/data/wilk_2020")
seu <- readRDS(paste0(seuName,".rds"))
seu <- run_cellstruct(seu, seuName, ref.proj, target.proj, clusterVar, nCores = 4) #GC and GS are saved as metadata of returned object, and a figure with panels of GC, GS, and cells population distribution is generated

#for visualization purpose
coi <- colnames(seu)[which(seu$cell.type.fine=="CD14 Monocyte")] #a vector of cell names belonging to CD14 Monocytes
heatmapGS(seu, ref.proj, target.proj, coi = coi) #generating heatmaps illustrating normalized cell-cell distances between provided cells and randomly selected 1,000 cells, in respective embeddings
dimredGS(seu, "pca", "umap", coi[1]) #plotting a UMAP dimension reduction showing the PCA distances between the marked cell (i.e. provided cell name) and all the other cells
heatmapGC(seu, ref.proj, target.proj) #generating heatmaps illustrating centroid-centroid distances between any two clusters, in respective embeddings
```
