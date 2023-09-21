#' dimension reduction GS
#'
#' Plotting the distance between a selected cell and all other cells in reference embedding,
#' on a target projection
#'
#' @param seu seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a target projection, eg. UMAP
#' @param coi cell_ID of a selected cell
#' @param dist_percentile the percentile used to normalize the distance between
#' a selected cell and all other cells. Defaults to .99
#'
#' @return NULL
#'
#' @author Jui Wan Loh
#'
#' @import data.table pdist Seurat ggplot2
#'
#' @examples
#' dimredGS(pbmc, "pca", "umap","cellID_1")
#'
#' @export

dimredGS <- function(seu, ref.proj, target.proj, coi, dist_percentile = .99){
  embed1 <- seu@reductions[[ref.proj]]@cell.embeddings
  dist.embed1 <- as.matrix(pdist(embed1[match(coi,rownames(embed1)),],embed1))
  prop.embed1 <- dist.embed1/apply(dist.embed1, 1, function (x) quantile(x, dist_percentile))
  prop.embed1[which(prop.embed1>1)] <- 1
  col_lim <- c(0,1)
  seu$PCAdist <- t(prop.embed1)
  if (target.proj == "spatial"){
    p1 <- SpatialDimPlot(seu, cells.highlight = coi) + NoLegend()
    p2 <- SpatialFeaturePlot(seu, features = "PCAdist") + scale_fill_distiller(palette = "RdYlBu", limits = col_lim)
    plot_grid(plotlist = list(p1,p2), nrow = 1)
  }else{
    FeaturePlot(seu, reduction=target.proj, features = "PCAdist", order = T, raster = T) +
      scale_color_distiller(palette = "RdYlBu", limits = col_lim) +
      geom_point(aes(x=seu@reductions[[target.proj]]@cell.embeddings[match(coi,rownames(embed1)),1],
                     y=seu@reductions[[target.proj]]@cell.embeddings[match(coi,rownames(embed1)),2]),
                 shape=4, size=1, stroke=1, colour="black")
  }
}
