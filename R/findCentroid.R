#' Find Centroid of a cluster
#'
#' Identifies the centroid of a cluster
#'
#' @param embed low-dimensional embedding of reference, a matrix with N cells
#'   and d1 dimensions.
#' @param df a matrix with one column storing 'seurat_clusters'
#'
#' @return a matrix with M clusters and d1 dimensions.
#'
#' @author Jui Wan Loh
#'
#' @import data.table
#'
#' @examples
#'centroid_embed1 <- findCentroid(inpPCA, clusID)
#'
#' @export

findCentroid <- function(embed, df){
  # Do checks
  if (nrow(embed) != nrow(df)){
    stop("embed and cluster matrix do not have same number of rows!")
  }

  addClus_embed <- merge(embed, df, by = 'row.names', all = T) #the order changes from embed
  centroid_embed <- c()
  clusNames <- unique(sort(df[,1]))
  for (i in 1:length(clusNames)){
    centroid_embed <- rbind(centroid_embed, colMeans(addClus_embed[which(addClus_embed[,ncol(addClus_embed)]==clusNames[i]),2:(ncol(addClus_embed)-1)]))
  }
  rownames(centroid_embed) <- clusNames
  return (centroid_embed)
}
