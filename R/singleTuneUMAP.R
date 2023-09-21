#' single run of runUMAP
#'
#' run Seurat::runUMAP function and write out the GS and LS for the parameters given
#'
#' @param seu seurat object
#' @param metric runUMAP parameters. Either "cosine" or "euclidean"
#' @param umap_nn runUMAP parameters for number of neighbors. Either 5, 15,30, or 50
#' @param min_dist runUMAP parameters for min distance. Either 0.05, 0.1, 0.3, or 0.5
#' @param ref.proj the reference projection, eg. PCA
#' @param outDir full directory of the output file
#' @param clusterVar the metadata that is to be used as clusterID for GC
#' @param num_waypoint number of waypoint cells used in calcGS. Defaults to 10000
#' @param cor_method correlation method in calculating GS score
#' @param num_neighbor number of nearest neighbours in calcLS. Defaults to 30
#' @param dist_power exponent to raise distance to in calcLS. Either 1 or 2
#' @param chunkSize size of partition for the dataset to be run in parallel
#' @param nCores number of cores used in parallel running
#' @param plot_LS boolean value deciding whether to calculate and plot LS.
#' Defaults to FALSE
#'
#' @return NULL
#'
#' @author Jui Wan Loh
#'
#' @import data.table Seurat
#'
#' @examples
#' singleTuneUMAP(PBMC,"euclidean",30,0.1,"pca","./tuneUMAP/","clusterID")
#'
#' @export

singleTuneUMAP <- function(seu, metric, umap_nn, min_dist, ref.proj, outDir, clusterVar, num_waypoint = 10000,
                     cor_method = "pearson",num_neighbor = 30, dist_power = 2,
                     chunkSize = 1000, nCores = 4, plot_LS = FALSE, seed = 42){
  ref <- seu@reductions[[ref.proj]]@cell.embeddings
  seu <- RunUMAP(seu, reduction = ref.proj, dims = 1:ncol(ref), metric = metric,
                 n.neighbors = umap_nn, min.dist = min_dist)
  data <- seu@reductions[["umap"]]@cell.embeddings[,1:2]
  seu$gs <- calcGS(ref, data, num_waypoint, cor_method, chunkSize, nCores, seed)
  outFile <- paste0(outDir,".",metric,"_",umap_nn,"_",min_dist,".txt")
  if (plot_LS){
    seu$ls <- calcLS(ref, data, num_neighbor, dist_power, chunkSize, nCores)
    write.table(data.frame(FetchData(object = seu, vars = c("gs","ls",clusterVar))), outFile, sep = "\t")
  }else{
    write.table(data.frame(FetchData(object = seu, vars = c("gs",clusterVar))), outFile, sep = "\t")
  }
}

