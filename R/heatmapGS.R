#' heatmap GS
#'
#' plotting the distance between selected and target cells in respective embeddings, in a heatmap
#'
#' @param seu seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a vector of target projections, eg. UMAP and tSNE
#' @param coi list of cell_ID of selected cells
#' @param dist_percentile the percentile used to normalize the distance between
#' a selected cell and target cells, i.e. row-normalization. Defaults to .99
#' @param target list of cell_ID of target cells. Defaults to NULL, i.e. list of random cells
#' @param num_waypoint number of waypoint cells used in calcGS. Defaults to 1000
#' @param param2 description of param2
#'
#' @return heatmap order of coi in its index format
#'
#' @author Jui Wan Loh
#'
#' @import data.table Seurat pdist ComplexHeatmap circlize viridis
#'
#' @examples
#' coiOrder <- heatmapGS(pbmc, "pca", c("umap","tsne"), c(cell_ID1, cell_ID2, cell_ID3))
#'
#' @export

heatmapGS <- function (seu, ref.proj, target.proj, coi, dist_percentile = .99,
                               target = NULL, num_waypoint = 1000, seed = 42){
  embed1 <- seu@reductions[[ref.proj]]@cell.embeddings
  embed2.list <- list()
  for (p in 1:length(target.proj)){
    embed2.list[[p]] <- seu@reductions[[target.proj[p]]]@cell.embeddings[,1:2]
    # Do checks
    if (nrow(embed1) != nrow(embed2.list[[p]])){
      stop("ref embedding and one of the target embeddings do not have same number of rows!")
    }
    if (nrow(embed2.list[[p]]) < num_waypoint){num_waypoint <- nrow(embed2.list[[p]])}
  }

  #Start calculation
  set.seed(seed)
  if (is.null(target)){randomCell <- sample(c(1:nrow(embed1)), size = num_waypoint)}
  else{randomCell <- match(target,rownames(embed1))}
  dist.embed1 <- as.matrix(pdist(embed1[match(coi,rownames(embed1)),],embed1[randomCell,]))
  prop.embed1 <- dist.embed1/apply(dist.embed1, 1, function (x) quantile(x, dist_percentile)) # distance percentile max per row
  prop.embed1[which(prop.embed1>1)] <- 1
  dist.embed2.list <- list()
  prop.embed2.list <- list()
  for (p in 1:length(embed2.list)){
    dist.embed2.list[[p]] <- as.matrix(pdist(embed2.list[[p]][match(coi,rownames(embed2.list[[p]])),],embed2.list[[p]][randomCell,]))
    prop.embed2.list[[p]] <- dist.embed2.list[[p]]/apply(dist.embed2.list[[p]], 1, function (x) quantile(x, dist_percentile))
    prop.embed2.list[[p]][which(prop.embed2.list[[p]]>1)] <- 1
  }

  #Start plotting
  col_ann = colorRamp2(c(seq(0, 1, by = 0.25)), inferno(5))
  col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  lgd_name <- "Norm_dist"

  exprs <- data.frame(FetchData(object = seu, vars = c(paste0(toupper(target.proj),rep("_gs",length(target.proj))),
                                                       "seurat_clusters")))
  randomCell.celltype <- exprs$seurat_clusters[randomCell]
  plot.embed1 <- Heatmap(as.matrix(prop.embed1), name = lgd_name, col = col_fun, column_title = toupper(ref.proj),
                         bottom_annotation = HeatmapAnnotation(Target_cluster=randomCell.celltype))
  cellOrder <- row_order(plot.embed1)
  plot.list <- NULL
  for (p in 1:length(prop.embed2.list)){
    plot.list <- plot.list + Heatmap(as.matrix(prop.embed2.list[[p]]), col = col_fun, show_heatmap_legend = F,
                                     column_order = column_order(plot.embed1), row_order = cellOrder,
                                     column_title = toupper(target.proj[p]))
  }

  exprs <- exprs[match(coi,rownames(exprs)),]
  for (ann in 1:length(prop.embed2.list)){
    plot.list <- plot.list + Heatmap(exprs[,ann], name = paste0(toupper(target.proj[ann]),"_GS"),
                                     col = col_ann, width = unit(5, "mm"))
  }
  plot.list <- plot.list + Heatmap(exprs$seurat_clusters, name = "Cluster",
                                   col = discreteColor(unique(exprs$seurat_clusters)), width = unit(5, "mm"))
  draw(plot.embed1 + plot.list, ht_gap = unit(1, "mm"))
  return(cellOrder)
}


