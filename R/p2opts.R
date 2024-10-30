#' @noRd

p2opts <- function(opts,pp) {
  if (!is.list(pp)) {
    opts$pt <- pp
    opts$p <- opts$ptrans(pp)
  } else {
    opts$p <- pp
    opts$pt <- p_to_optim(pp)$values
  }
  return(opts)
}
