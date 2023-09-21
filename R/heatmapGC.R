#' heatmap GC
#'
#' plotting the centroid-centroid distance in respective embeddings, in a heatmap
#'
#' @param seu seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a vector of target projections, eg. UMAP and tSNE
#' @param dist_percentile the percentile used to normalize the centroid-centroid
#' distances. Defaults to 1
#'
#' @return NULL
#'
#' @author Jui Wan Loh
#'
#' @import data.table viris Seurat ComplexHeatmap circlize
#'
#' @examples
#' heatmapGC(pbmc, "pca", c("umap","tsne"))
#'
#' @export

heatmapGC <- function (seu, ref.proj, target.proj, dist_percentile = 1){
  embed1 <- seu@reductions[[ref.proj]]@cell.embeddings
  embed2.list <- list()
  for (p in 1:length(target.proj)){
    embed2.list[[p]] <- seu@reductions[[target.proj[p]]]@cell.embeddings[,1:2]
    # Do checks
    if (nrow(embed1) != nrow(embed2.list[[p]])){
      stop("ref embedding and one of the target embeddings do not have same number of rows!")
    }
  }

  exprs <- data.frame(FetchData(object = seu, vars = c(paste0(toupper(target.proj),rep("_gc",length(target.proj))),
                                                       "seurat_clusters")))
  #Start calculation
  centroid_embed1 <- findCentroid(embed1, exprs[ncol(exprs)])
  centroid_embed2.list <- list()
  for (p in 1:length(embed2.list)){
    centroid_embed2.list[[p]] <- findCentroid(embed2.list[[p]], exprs[ncol(exprs)])
  }
  exprs <- unique(exprs)
  exprs <- exprs[match(unique(sort(exprs$seurat_clusters)),exprs$seurat_clusters),]

  dist.centroid_embed1 <- as.matrix(dist(centroid_embed1, method = "euclidean"))
  prop.centroid_embed1 <- dist.centroid_embed1/quantile(dist.centroid_embed1, dist_percentile)
  prop.centroid_embed1[which(prop.centroid_embed1>1)] <- 1
  dist.centroid_embed2.list <- list()
  prop.centroid_embed2.list <- list()
  for (p in 1:length(embed2.list)){
    dist.centroid_embed2.list[[p]] <- as.matrix(dist(centroid_embed2.list[[p]], method = "euclidean"))
    prop.centroid_embed2.list[[p]] <- dist.centroid_embed2.list[[p]]/quantile(dist.centroid_embed2.list[[p]], dist_percentile)
    prop.centroid_embed2.list[[p]][which(prop.centroid_embed2.list[[p]]>1)] <- 1
  }

  #Start plotting
  col_ann = colorRamp2(c(seq(0, 1, by = 0.25)), inferno(5))
  col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  lgd_name <- "Norm_dist"

  plot.embed1 <- Heatmap(as.matrix(prop.centroid_embed1), name = lgd_name,
                         col = col_fun, column_title = toupper(ref.proj))
  cellOrder <- row_order(plot.embed1)
  plot.list <- NULL
  for (p in 1:length(embed2.list)){
    plot.list <- plot.list + Heatmap(as.matrix(prop.centroid_embed2.list[[p]]), col = col_fun,
                                     show_heatmap_legend = F, column_order = column_order(plot.embed1),
                                     row_order = cellOrder, column_title = toupper(target.proj[p]))
  }
  for (ann in 1:length(embed2.list)){
    plot.list <- plot.list + Heatmap(exprs[,ann], name = paste0(toupper(target.proj[ann]),"_GC"),
                                     col = col_ann, width = unit(5, "mm"))
  }
  plot.list <- plot.list + rowAnnotation(rn = anno_text(rownames(prop.centroid_embed1)))
  draw(plot.embed1 + plot.list, ht_gap = unit(1, "mm"))
}

