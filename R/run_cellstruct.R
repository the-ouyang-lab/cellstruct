#' run cellstruct
#'
#' Calculates the four metric scores
#'
#' @param seu seurat object
#' @param seuName name of seurat object, including directory
#' @param ref.proj the reference projection, eg. PCA
#' @param target.proj a vector of target projections, eg. UMAP and tSNE
#' @param clusterVar the metadata that is to be used as clusterID for GC
#' @param num_waypoint number of waypoint cells used in calcGS. Defaults to 10000
#' @param cor_method correlation method in calculating GS and GC score
#' @param num_neighbor number of nearest neighbours in calcLS. Defaults to 30
#' @param dist_power exponent to raise distance to in calcLS. Either 1 or 2
#' @param chunkSize size of partition for the dataset to be run in parallel
#' @param nCores number of cores used in parallel running
#' @param plot_LS boolean value deciding whether to calculate and plot LS.
#' Defaults to FALSE
#'
#' @return seurat object with scores saved as metadata (eg. UMAP_gs)
#'
#' @author Jui Wan Loh
#'
#' @import Seurat ggplot2 cowplot
#'
#' @examples
#' new_seu <- run_cellstruct(pbmc,"./data/seu","pca",c("umap","tsne"),"clusterID")
#'
#' @export


run_cellstruct <- function(seu, seuName, ref.proj, target.proj, clusterVar,
                           num_waypoint = 10000, cor_method = "pearson",
                           num_neighbor = 30, dist_power = 2, chunkSize = 1000,
                           nCores = 4, plot_LS = FALSE, seed = 42){
  ref <- seu@reductions[[ref.proj]]@cell.embeddings
  seu$seurat_clusters <- seu@meta.data[clusterVar]  #the cluster ID must be called "seurat_clusters"

  clusID <- data.frame(FetchData(object = seu, vars = c("seurat_clusters")))
  final.list <- list()
  for(proj in target.proj){
    print(proj)
    print(date())
    data <- seu@reductions[[proj]]@cell.embeddings

    seu <- AddMetaData(object = seu,
                      metadata = calcGS(ref, data, num_waypoint, cor_method, nCores = nCores, seed = seed),
                      col.name = paste0(toupper(proj),"_gs"))
    if(plot_LS){
      seu <- AddMetaData(object = seu,
                         metadata = calcLS(ref, data, num_neighbor, dist_power, nCores = nCores),
                         col.name = paste0(toupper(proj),"_ls"))
    }
    gc_pearson <- calcGC(ref,data,clusID, cor_method) #~1 sec
    diag(gc_pearson)[which(diag(gc_pearson)<0)] <- 0
    seu <- AddMetaData(object = seu, metadata = rep(NA,ncol(seu)),
                       col.name = paste0(toupper(proj),"_gc"))
    for (id in unique(clusID[,1])){
      seu@meta.data[paste0(toupper(proj),"_gc")][which(seu$seurat_clusters==id),] <- diag(gc_pearson)[id]
    }

    if(plot_LS){plot_var = paste0(rep(toupper(proj),3),c("_gc","_gs","_ls"))}
    else{plot_var = paste0(rep(toupper(proj),2),c("_gc","_gs"))}
    exprs <- data.frame(FetchData(seu, vars = plot_var))
    plot.list <- list()
    for (p in 1:length(plot_var)){
      if (p == 1){
        plot.list[[p]] <- plot_withCentroid(seu, data, clusID, proj, plot_var[p]) +
          ggtitle(paste0(plot_var[p],": ",round(mean(exprs[,p]),4))) + NoLegend()
      }
      else{
        plot.list[[p]] <- FeaturePlot(seu, reduction = proj, features = plot_var[p], order = T) +
          ggtitle(paste0(plot_var[p],": ",round(mean(exprs[,p]),4))) + NoLegend() +
          scale_color_distiller(palette = "RdYlBu", limits = c(0,1))
      }
    }
    plot.list[[p+1]] <- DimPlot(seu, reduction = proj, group.by = clusterVar) + NoLegend()
    final.list <- c(final.list, plot.list)
  }
  lgd <- get_legend(
    FeaturePlot(seu, reduction = proj, features = plot_var[p], order = T) +
      ggtitle(paste0(plot_var[p],": ",round(mean(exprs[,p]),4))) +
      scale_color_distiller(palette = "RdYlBu", limits = c(0,1)) +
      theme(legend.position = "bottom",legend.text=element_text(size=8))
  )
  pdf(paste0(seuName,"_metricScores.pdf"), w=(3*(length(plot_var)+1)), h=(3*length(target.proj)))
  pt <- plot_grid(plotlist = final.list, nrow = length(target.proj))
  print(plot_grid(pt, lgd, ncol = 1, rel_heights = c(1, .1)))
  dev.off()
  return (seu)
}

