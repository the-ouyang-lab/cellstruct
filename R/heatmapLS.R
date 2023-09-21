#' heatmap LS
#'
#' plotting the corresponding LS contribution of each ordered neighbor for respective target projection,
#' in respective heatmaps
#'
#' @param seu seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a vector of target projections, eg. UMAP and tSNE
#' @param coi list of cell_ID of selected cells
#' @param num_neighbor number of nearest neighbours. Defaults to 30
#' @param dist_power exponent to raise distance to. Either 1 or 2
#'
#' @return heatmap order of coi in its index format
#'
#' @author Jui Wan Loh
#'
#' @import data.table RANN Seurat ComplexHeatmap circlize pdist viridis
#'
#' @examples
#' coiOrder <- heatmapLS(pbmc, "pca", c("umap","tsne"), c(cell_ID1, cell_ID2, cell_ID3))
#'
#' @export

heatmapLS <- function (seu, ref.proj, target.proj, coi, num_neighbor = 30, dist_power = 2){
  embed1 <- seu@reductions[[ref.proj]]@cell.embeddings
  embed2.list <- list()
  for (p in 1:length(target.proj)){
    embed2.list[[p]] <- seu@reductions[[target.proj[p]]]@cell.embeddings[,1:2]
    # Do checks
    if (nrow(embed1) != nrow(embed2.list[[p]])){
      stop("ref embedding and one of the target embeddings do not have same number of rows!")
    }
  }
  if (num_neighbor < 15){stop("num_neighbor must be >= 15!")}
  if (num_neighbor > 500){stop("num_neighbor must be <= 500!")}
  if (!dist_power %in% c(1,2)){stop("dist_power should be either 1 or 2!")}
  if (num_neighbor >= nrow(embed1)){num_neighbor <- nrow(embed1)-1}

  #Start calculation
  coi_idx <- match(coi,rownames(embed1))
  nearest.embed1 <- nn2(embed1, k = num_neighbor+1)
  nearest.embed2.list <- list()
  for (p in 1:length(embed2.list)){
    nearest.embed2.list[[p]] <- nn2(embed2.list[[p]], k = num_neighbor+1)
  }
  embed2.neighbor_idx.list <- vector("list",length(embed2.list))
  embed1.neighbor_idx.list <- list()
  for (n in 1:length(coi_idx)){
    embed1.neighbor_idx.list[[n]] <- nearest.embed1$nn.idx[coi_idx[n],2:ncol(nearest.embed1$nn.idx)]
    for (p in 1:length(embed2.list)){
      embed2.neighbor_idx.list[[p]][[n]] <- nearest.embed2.list[[p]]$nn.idx[coi_idx[n],2:ncol(nearest.embed2.list[[p]]$nn.idx)]
    }
  }

  prop.embed2.list <- vector("list", length(embed2.list))
  for (p in 1:length(embed2.list)){
    for (n in 1:length(coi_idx)){
      #distance of coi and UMAP neighbor in PCA space
      numerator <- 1/pdist(embed1[coi_idx[n],],embed1[embed2.neighbor_idx.list[[p]][[n]],])@dist**dist_power
      #distance of coi and PCA neighbor in PCA space
      denominator <- mean(1/pdist(embed1[coi_idx[n],],embed1[embed1.neighbor_idx.list[[n]],])@dist**dist_power)
      tmp <- numerator/denominator
      prop.embed2.list[[p]] <- rbind(prop.embed2.list[[p]],tmp)
    }
    prop.embed2.list[[p]][which(prop.embed2.list[[p]]>1)] <- 1
  }

  #Start plotting
  col_ann = colorRamp2(c(seq(0, 1, by = 0.25)), inferno(5))
  col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  lgd_name <- "Norm_ratio"

  ht1 <- Heatmap(as.matrix(prop.embed2.list[[1]]), col = col_fun, name = lgd_name,
                 cluster_columns = F, column_title = toupper(target.proj[1]))
  plot.list <- NULL
  if (length(prop.embed2.list) > 1){
    for (p in 2:(length(prop.embed2.list))){
      plot.list <- plot.list + Heatmap(as.matrix(prop.embed2.list[[p]]), col = col_fun,
                                       show_heatmap_legend = F, cluster_columns = F, column_title = toupper(target.proj[p]))
    }
  }
  exprs <- data.frame(FetchData(object = seu, vars = c(paste0(toupper(target.proj),rep("_ls",length(target.proj))),
                                                       "seurat_clusters")))
  exprs <- exprs[match(coi,rownames(exprs)),]
  for (ann in 1:length(target.proj)){
    plot.list <- plot.list + Heatmap(exprs[,ann], name = paste0(toupper(target.proj[ann]),"_LS"),
                                     col = col_ann, width = unit(5, "mm"))
  }
  plot.list <- plot.list + Heatmap(exprs$seurat_clusters, name = "cluster",
                                   col = discreteColor(unique(exprs$seurat_clusters)), width = unit(5, "mm"))
  draw(ht1 + plot.list, ht_gap = unit(1, "mm"))
  return (row_order(ht1))
}

