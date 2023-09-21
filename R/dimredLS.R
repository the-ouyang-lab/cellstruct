#' dimension reduction LS
#'
#' plotting the corresponding LS contribution of each neighbor on respective target projection
#'
#' @param seu seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a target projection, eg. UMAP
#' @param coi cell_ID of a (or a few) selected cell(s)
#' @param num_neighbor number of nearest neighbours. Defaults to 30
#' @param dist_power exponent to raise distance to. Either 1 or 2
#'
#' @return metric score
#'
#' @author Jui Wan Loh
#'
#' @import data.table RANN Seurat pdist ggrepel
#'
#' @examples
#' dimredLS(pbmc, "pca", "umap", c(cellID_1, cellID_2))
#'
#' @export

dimredLS <- function(seu, ref.proj, target.proj, coi, num_neighbor = 30, dist_power = 2){
  embed1 <- seu@reductions[[ref.proj]]@cell.embeddings
  embed2 <- seu@reductions[[target.proj]]@cell.embeddings[,1:2]

  # Do checks
  if (nrow(embed1) != nrow(embed2)){
    stop("embed1 and embed2 does not have same number of rows!")
  }
  if (num_neighbor < 15){stop("num_neighbor must be >= 15!")}
  if (num_neighbor > 500){stop("num_neighbor must be <= 500!")}
  if (!dist_power %in% c(1,2)){stop("dist_power should be either 1 or 2!")}
  if (num_neighbor >= nrow(embed1)){num_neighbor <- nrow(embed1)-1}

  nearest.embed1 <- nn2(embed1, k = num_neighbor+1)
  nearest.embed2 <- nn2(embed2, k = num_neighbor+1)
  coi_idx <- match(coi,rownames(embed2))
  embed1.neighbor_idx <- matrix(nearest.embed1$nn.idx[coi_idx,2:ncol(nearest.embed1$nn.idx)], nrow = length(coi))
  embed2.neighbor_idx <- matrix(nearest.embed2$nn.idx[coi_idx,2:ncol(nearest.embed2$nn.idx)], nrow = length(coi))
  seu$neighbors <- rep(NA,length(colnames(seu)))
  for (i in 1:length(coi_idx)){
    dist.embed1 <- pdist(embed1[coi_idx[i],],embed1[embed1.neighbor_idx[i,],])@dist**dist_power
    dist.embed2 <- pdist(embed1[coi_idx[i],],embed1[embed2.neighbor_idx[i,],])@dist**dist_power
    prop.embed <- (1/dist.embed2)/mean(1/dist.embed1)
    prop.embed[which(prop.embed>1)] <- 1
    seu$neighbors[embed2.neighbor_idx[i,]] <- t(prop.embed)
  }
  col_lim <- c(0,1)
  if (target.proj == "spatial"){
    SpatialFeaturePlot(seu, features = "neighbors") + scale_fill_distiller(palette = "RdYlBu", limits = col_lim)
  }else{
    if (length(coi_idx) == 1){markedCells.loc <- as.data.frame(t(embed2[coi_idx,]))}
    else{markedCells.loc <- as.data.frame(embed2[coi_idx,])}
    colnames(markedCells.loc) <- c("dim1","dim2")
    FeaturePlot(seu, reduction=target.proj, features = "neighbors", order = T, pt.size = 0.3) +
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = col_lim, na.value = "grey90") +
      geom_point(data = markedCells.loc, aes(x=dim1, y=dim2), shape=4, size=0.1, stroke=0.5, colour="black") +
      geom_label_repel(data = markedCells.loc, aes(label=coi, x=dim1, y=dim2), box.padding = 3)
  }
}

