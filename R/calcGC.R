#' Global Cluster score
#'
#' Calculates the global cluster score by correlating the pairwise distances
#' between any two cluster centroids
#'
#' @param embed1 low-dimensional embedding of reference, a matrix with N cells
#'   and d1 dimensions.
#' @param embed2 low-dimensional embedding of target, a matrix with N cells
#'   and d2 dimensions. Note that d1 and d2 need not be the same.
#' @param df a matrix with one column storing 'seurat_clusters'
#' @param cor_method correlation method in calculating GC score
#'
#' @return metric score
#'
#' @author Jui Wan Loh
#'
#' @import data.table
#'
#' @examples
#' gcScore <- calcGC(inpPCA, inpUMAP, clusID)
#'
#' @export


calcGC <- function(embed1, embed2, df, cor_method = "pearson"){
  # Do checks
  if (nrow(embed1) != nrow(embed2)){
    stop("embed1 and embed2 do not have same number of rows!")
  }
  if (nrow(embed1) != nrow(df)){
    stop("embed1 and cluster matrix do not have same number of rows!")
  }

  #Start calculation
  centroid_embed1 <- findCentroid(embed1, df)
  centroid_embed2 <- findCentroid(embed2, df)
  dist.centroid_embed1 <- as.matrix(dist(centroid_embed1, method = "euclidean"))
  dist.centroid_embed2 <- as.matrix(dist(centroid_embed2, method = "euclidean"))
  corr.mat <- cor(dist.centroid_embed1, dist.centroid_embed2, method = cor_method)
  return (corr.mat)
}
