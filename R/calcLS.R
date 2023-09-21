#' Local Single cell score
#'
#' Calculates the local single cell score by comparing the distance of the k
#' nearest neighbours in two embeddings.
#'
#' @param embed1 low-dimensional embedding of reference, a matrix with N cells
#'   and d1 dimensions.
#' @param embed2 low-dimensional embedding of target, a matrix with N cells
#'   and d2 dimensions. Note that d1 and d2 need not be the same.
#' @param num_neighbor number of nearest neighbours. Defaults to 30
#' @param dist_power exponent to raise distance to. Either 1 or 2
#' @param chunkSize size of partition for the dataset to be run in parallel
#' @param nCores number of cores used in parallel running
#'
#' @return vector of LS score
#'
#' @author Jui Wan Loh
#'
#' @import data.table RANN parallel
#'
#' @examples
#' lsScore <- calcLS(inpPCA, inpUMAP)
#'
#' @export

calcLS <- function(embed1, embed2, num_neighbor = 30, dist_power = 2,
                   chunkSize = 1000, nCores = 4) {
  # Do checks
  if (nrow(embed1) != nrow(embed2)){
    stop("embed1 and embed2 does not have same number of rows!")
  }
  if (num_neighbor < 15){stop("num_neighbor must be >= 15!")}
  if (num_neighbor > 500){stop("num_neighbor must be <= 500!")}
  if (!dist_power %in% c(1,2)){stop("dist_power should be either 1 or 2!")}
  if (num_neighbor >= nrow(embed1)){num_neighbor <- nrow(embed1)-1}

  nChunk = nrow(embed1) %/% chunkSize
  if(nrow(embed1) %% chunkSize > 5){nChunk = nChunk + 1}
  # Merge last chunk if too small
  tmpChunk1 = seq(1, nChunk * chunkSize, chunkSize)
  tmpChunk2 = tmpChunk1 + chunkSize - 1
  tmpChunk2[nChunk] = nrow(embed1)
  avgDist.mat <- do.call(
    rbind, mclapply(seq(nChunk), function(iChunk){
      # Start calculation (plus1 because nn2 finds itself)
      nearest.embed1 <- nn2(embed1, k = num_neighbor+1,
                            query = embed1[tmpChunk1[iChunk]:tmpChunk2[iChunk],])
      nearest.embed2 <- nn2(embed2, k = num_neighbor+1,
                            query = embed2[tmpChunk1[iChunk]:tmpChunk2[iChunk],])
      avgDist.embed1.vec <- numeric(length(tmpChunk1[iChunk]:tmpChunk2[iChunk]))
      avgDist.embed2.vec <- numeric(length(tmpChunk1[iChunk]:tmpChunk2[iChunk]))
      for (i in 1:length(avgDist.embed1.vec)){
        k <- c(tmpChunk1[iChunk]:tmpChunk2[iChunk])[i]
        # up to 4-5 decimal place correct with the nn.dists from nn2 function
        avgDist.embed1.vec[i] <- mean(1/pdist(embed1[k,],embed1[
          nearest.embed1$nn.idx[i,2:ncol(nearest.embed1$nn.idx)],])@dist**dist_power)
        # calculating distance in embed1 space for embed2 neighbors
        avgDist.embed2.vec[i] <- mean(1/pdist(embed1[k,],embed1[
          nearest.embed2$nn.idx[i,2:ncol(nearest.embed2$nn.idx)],])@dist**dist_power)
      }
    return(cbind(avgDist.embed1.vec,avgDist.embed2.vec))
  }, mc.cores = nCores))

  out.vec <- avgDist.mat[,2]/avgDist.mat[,1]
  names(out.vec) <- rownames(embed1)
  return (out.vec)
}


