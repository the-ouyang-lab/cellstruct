#' discrete colors
#'
#' Generate a list of discrete colors
#'
#' @param id_vec a vector to be associated with different color for each element
#'
#' @return a color vector
#'
#' @author Jui Wan Loh
#'
#' @import RColorBrewer
#'
#' @examples
#' colorVec <- discrete_colors(clusVec)
#'
#' @export

discreteColor <- function(id_vec){
  n <- length(id_vec)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vec = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))),n)
  dict <- new.env(hash = T)
  Add <- function(key,val) dict[[key]] <- val
  return (mapply(Add, as.character(id_vec), col_vec))
}

