#' @noRd

derivtransfunchooser <- function(fs) { # Returns appropriate back-transformation or derivative functions based on a list of transformation names.
  translist <- c(log=exp,
                 covlogit=function(x) 2*exp(x)/(1+exp(x))^2,
                 logit=function(x) exp(x)/(1+exp(x))^2,
                 identity=function(x) 1)
  res <- translist[fs]
  for (i in which(map_lgl(res,is.null))) {
    names(res)[i]  <- fs[i]
    res[[i]] <- function(x) grad(fs[i],x)
  }
  res
}
