#' @noRd
#'

backtransfunchooser <- function(fs) { # Returns appropriate back-transformation or derivative functions based on a list of transformation function.
  translist <- c(log=exp,
                 covlogit=anticovlogit,
                 logit=function(x) exp(x)/(1+exp(x)),
                 identity=identity)
  res <- translist[fs]
  for (i in which(map_lgl(res,is.null))) {
    names(res)[i]  <- fs[i]
    res[[i]] <- numtransfun(fs[i])
  }
  res
}
