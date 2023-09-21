#' featurePlot with Centroid marked
#'
#' Generates featurePlot for GC and LC score with centroids marked
#'
#' @param seu seurat object
#' @param embed low-dimensional embedding of reference, a matrix with N cells
#' and d1 dimensions.
#' @param df a matrix with one column storing 'seurat_clusters'
#' @param proj dimension reduction that is stored in the seurat object.
#' Choices are "umap", "fdl", and "tsne".
#' @param score metadata stored in the seurat object that we want to show on
#' featurePlot. Eg. "UMAP_gc"
#'
#' @return featurePlot
#'
#' @author Jui Wan Loh
#'
#' @import data.table Seurat
#'
#' @examples
#' plot_withCentroid(pbmc, inpUMAP, clusID, "umap", "UMAP_gc")
#'
#' @export

plot_withCentroid <- function(seu, embed, df, proj, score){
  centroid_embed <- as.data.frame(findCentroid(embed, df))
  colnames(centroid_embed) <- c("dim1","dim2")
  p1 <- FeaturePlot(seu, reduction=proj, features = score, order = T) +
    scale_color_distiller(palette = "RdYlBu", limits = c(0,1)) +
    geom_point(data = centroid_embed, aes(x=dim1, y=dim2), shape=4, size=1, stroke=1, colour="black")
  return(p1)
}

