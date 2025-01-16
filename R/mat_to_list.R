#' @noRd

mat_to_list <- function(x) { # Converts a matrix back into a list, extracting diagonal and off-diagonal blocks.
  ## transform a sparse matrix to a list of values, similar to list_to_mat
  if (nrow(x)!=ncol(x)) stop("Need square matrix!")
  if (nrow(x)==1) return(list(c(x)))
  check <- apply(x,2,function(i) sum(rle(i!=0)$values)>1)
  if (any(check)) stop("Matrix needs to be either full, band or diagonal; all other types are unsupported!")
  i <- 1
  res <- list()
  while(i<=ncol(x)) {
    l <- sum(x[,i]!=0)
    if (l==0)  l <- 1 ## for the special case of zero omega
    ind <- i-1+seq_len(l)
    newcomp <- c(x[ind,ind][upper.tri(x[ind,ind],diag=TRUE)])
    res <- c(res,list(newcomp))
    i <- i+l
  }
  res
}
