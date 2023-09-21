#' Global Single cell score
#'
#' Calculates the global single cell score by correlating the distance between a
#' single cell and other waypoint cells
#'
#' @param embed1 low-dimensional embedding of reference, a matrix with N cells
#'   and d1 dimensions.
#' @param embed2 low-dimensional embedding of target, a matrix with N cells
#'   and d2 dimensions. Note that d1 and d2 need not be the same.
#' @param num_waypoint number of waypoint cells. Defaults to 10000
#' @param cor_method correlation method in calculating GS score
#' @param chunkSize size of partition for the dataset to be run in parallel
#' @param nCores number of cores used in parallel running
#'
#'
#' @return vector of GS score
#'
#' @author Jui Wan Loh
#'
#' @import data.table pdist parallel
#'
#' @examples
#' gsScore <- calcGS(inpPCA, inpUMAP)
#'
#' @export

calcGS <- function(embed1, embed2, num_waypoint = 10000, cor_method = "pearson",
                   chunkSize = 1000, nCores = 4, seed = 42){
  # Do checks
  if (nrow(embed1) != nrow(embed2)){
    stop("embed1 and embed2 does not have same number of rows!")
  }
  if (nrow(embed2) < num_waypoint){num_waypoint <- nrow(embed2)}

  # Start calculation
  set.seed(seed)
  randomCell <- sample(c(1:nrow(embed1)), size = num_waypoint)
  nChunk = nrow(embed1) %/% chunkSize
  if(nrow(embed1) %% chunkSize > 5){nChunk = nChunk + 1}
  # Merge last chunk if too small
  tmpChunk1 = seq(1, nChunk * chunkSize, chunkSize)
  tmpChunk2 = tmpChunk1 + chunkSize - 1
  tmpChunk2[nChunk] = nrow(embed1)
  corr.vec <- do.call(
    c, mclapply(seq(nChunk), function(iChunk){
      dist.embed1 <- t(as.matrix(pdist(embed1[tmpChunk1[iChunk]:tmpChunk2[iChunk],],embed1[randomCell,])))
      dist.embed2 <- t(as.matrix(pdist(embed2[tmpChunk1[iChunk]:tmpChunk2[iChunk],],embed2[randomCell,])))
      vec <- numeric(length(tmpChunk1[iChunk]:tmpChunk2[iChunk]))
      for (i in 1:length(vec)){  # calculating correlation for diagonal of two matrices
        vec[i] <- cor(dist.embed1[,i], dist.embed2[,i], method = cor_method)
      }
      return(vec)
    }, mc.cores = nCores))

  corr.vec[corr.vec < 0] <- 0
  names(corr.vec) <- rownames(embed1)
  return (corr.vec)
}


