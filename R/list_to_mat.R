#' @noRd
#'
list_to_mat <- function(l) { # Converts a list of values into a matrix by organizing them into a block-diagonal format.
  ## transform a list of values to a band matrix, similar to a series of
  ## NONMEM $OMEGA blocks
  mdim <- map_dbl(l,function(i) 0.5*(sqrt(1+8*length(i))-1))
  m <- matrix(0,sum(mdim),sum(mdim))
  cums <- cumsum(c(0,head(mdim,-1)))
  for (i in seq_along(mdim)) {
    indices <- cums[i]+seq_len(mdim[i])
    mtemp <- matrix(0,mdim[i],mdim[i])
    mtemp[upper.tri(mtemp,diag=TRUE)] <- l[[i]]
    m[indices,indices] <- uptolow(mtemp)
  }
  m
}
