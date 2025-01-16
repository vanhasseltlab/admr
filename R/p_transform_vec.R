#' @noRd


p_transform_vec <- function(x) { # Transforms a parameter vector for optimization, handling transformations like log, covlogit, and fixed parameters.
  ### Transforms a vector of input values into something that can be fed to optim
  ## params with "fix" are removed
  origx <- x
  fixed <- str_detect(x,"fix") | str_detect(x,"Fix") | str_detect(x,"FIX")
  x <- x[!fixed]
  ## if other strings are present, then they are assumed to relate to the desired transformation
  xtransfuns <- rep("log",length(x))
  othertransfuns <- str_extract(format(x,scientific=FALSE),"[[:alpha:]]+")
  xtransfuns[!is.na(othertransfuns)] <- na.omit(othertransfuns)[1]
  x2 <- str_extract(x,".*[[:digit:]]+") %>% as.numeric()
  values <- map_dbl(seq_along(x2),~do.call(xtransfuns[.],list(x2[.])))
  backtransseq <- backtransfunchooser(xtransfuns)
  derivseq <- derivtransfunchooser(xtransfuns)
  backtransformfunc <- function(p,numeric=TRUE) {
    if (missing(p)) p <- values
    p2 <- map_dbl(seq_along(p),~backtransseq[[.]](p[.]))
    res <- origx
    res[!fixed] <- p2
    if (numeric & !is.numeric(res)) as.numeric(str_extract(res,"[[:digit:] [:punct:]]+")) else res
  }
  d_psi_d_psitrans_long <- function(p) {
    if (missing(p)) p <- values
    res <- rep(0,length(origx))
    map(seq_along(p),function(i) {
      ii <- which(!fixed)[i]
      res[ii] <- derivseq[[i]](p[i])
      res
    })
  }
  d_psi_d_psitrans_short <- function(p) {
    if (missing(p)) p <- values
    map2(derivseq,p,~.x(.y)) %>% unlist()
  }
  list(values=values,backtransformfunc=backtransformfunc,transfunclist=xtransfuns,origx=origx,
       d_psi_d_psitrans_long=d_psi_d_psitrans_long,d_psi_d_psitrans_short=d_psi_d_psitrans_short)
}
