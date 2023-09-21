#' tune runUMAP
#'
#' Changes the parameters of Seurat::runUMAP function to select the best
#' parameters of metrics, n.neighbors, and min.dist based on GS score
#'
#' @param seu seurat object
#' @param seuName name of seurat object
#' @param ref.proj the reference projection, eg. PCA
#' @param outFile full directory of the output file
#' @param clusterVar the metadata that is to be used as clusterID for GC
#' @param weight the assigned weight for each cluster in the metadata provided as
#' clusterID for GC. Defaults to NULL (i.e. equal weight for all clusters)
#' @param num_waypoint number of waypoint cells used in calcGS. Defaults to 10000
#' @param cor_method correlation method in calculating GS score
#' @param num_neighbor number of nearest neighbours in calcLS. Defaults to 30
#' @param dist_power exponent to raise distance to in calcLS. Either 1 or 2
#' @param chunkSize size of partition for the dataset to be run in parallel
#' @param nCores number of cores used in parallel running
#' @param plot_LS boolean value deciding whether to calculate and plot LS.
#' Defaults to FALSE
#'
#'
#' @return seurat object with runUMAP from best parameters, along with metric
#' scores saved as metadata
#'
#' @author Jui Wan Loh
#'
#' @import data.table Seurat
#'
#' @examples
#' new_seu <- tuneUMAP(pbmc,"pbmc","pca","./tuningUMAP_param.txt","clusterID")
#'
#' @export

tuneUMAP <- function(seu, seuName, ref.proj, outFile, clusterVar, weight = NULL, num_waypoint = 10000,
                     cor_method = "pearson",num_neighbor = 30, dist_power = 2,
                     chunkSize = 1000, nCores = 4, plot_LS = FALSE, seed = 42){
  outname <- paste0(outFile,".weighted_param.txt")
  if(is.null(weight)){
    df <- data.frame(FetchData(object = seu, vars = clusterVar))
    weight <- cbind.data.frame(clus=unique(sort(df[,1])),weight=rep(1,length(unique(df[,1]))))
    outname <- paste0(outFile,".unweighted_param.txt")
  }
  ref <- seu@reductions[[ref.proj]]@cell.embeddings

  metric <- c("euclidean","cosine")
  umap_neighbors <- c(15,30,50)
  min_dist <- c(0.05, 0.1, 0.3, 0.5)
  out.mat <- c()
  for (x in metric){
    for (y in umap_neighbors){
      for (z in min_dist){
        singleTuneUMAP(seu, x, y, z, ref.proj, outFile, clusterVar, num_waypoint,
                       cor_method, num_neighbor, dist_power, chunkSize, nCores, plot_LS, seed)
        tmp <- read.table(paste0(outFile,".",x,"_",y,"_",z,".txt"), sep = "\t", header = T)
        tmp$weight <- weight$weight[match(tmp[[clusterVar]],weight$clus)]
        if(plot_LS){
          out.mat <- rbind(out.mat, cbind.data.frame(x,y,z,mean(tmp$gs*tmp$weight),mean(tmp$ls*tmp$weight)))
        }else{
          out.mat <- rbind(out.mat, cbind.data.frame(x,y,z,mean(tmp$gs*tmp$weight)))
        }
      }
    }
  }
  if(plot_LS){colnames(out.mat) <- c("metric","n.neighbors","min.dist","mean_GS","mean_LS")}
  else{colnames(out.mat) <- c("metric","n.neighbors","min.dist","mean_GS")}
  write.table(out.mat, outname, sep = "\t", row.names = F)

  idx <- which(out.mat$mean_GS == max(out.mat$mean_GS))
  seu <- RunUMAP(seu, reduction = ref.proj, dims = 1:ncol(ref),
                 metric = out.mat[idx,1], n.neighbors = out.mat[idx,2], min.dist = out.mat[idx,3])
  final_seu <- run_cellstruct(seu, seuName, ref.proj, "umap", clusterVar,
                              num_waypoint, cor_method, num_neighbor, dist_power,
                              chunkSize, nCores, plot_LS, seed)
  return(final_seu)
}

