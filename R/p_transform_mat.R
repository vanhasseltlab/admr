#' @noRd
#'
p_transform_mat <- function(x) { # Transforms a matrix (typically representing variances and covariances) for optimization by converting it to a vector with appropriate transformations applied.
  ### Converts a matrix to a form that can be fed to optim
  ## matrix to list-format
  origx <- x
  if (is.matrix(x)) x <- mat_to_list(x)
  sames <- map_lgl(x,~all(. %in% c("same","Same","SAME")))
  if (length(sames)>0 && sames[1]) stop("Cannot start with the _same_ argument!")
  pn <- 0
  matmapping <- x
  for (i in seq_along(x)) {
    if (any(str_detect(x[[i]],"fix|Fix|FIX"))) {
      numbers <- parse_number(x[[i]])
      matmapping[[i]] <- paste0(numbers," fix")
      x[[i]] <- numbers
    } else if (!sames[i]) {
      pnseq <- pn+seq_along(x[[i]])
      pn <- tail(pnseq,1)
      matmapping[[i]] <- paste0("p",pnseq)
    } else {
      matmapping[[i]] <- matmapping[[i-1]]
      x[[i]] <- x[[i-1]]
    }
  }
  ## x2 <- list_to_mat(x)
  x2 <- list_to_mat(map(x,as.numeric))
  x3 <- cov2cor(x2);diag(x3) <- diag(x2) ## retain diagonals, offdiagonals to correlations
  matmapping <- list_to_mat(matmapping)
  parlocs <- which(grepl("p",c(matmapping)) & !duplicated(c(matmapping)))
  transformlist <- rep("log",length(parlocs))
  for (i in 2:length(transformlist)) {
    matind <- which(matmapping==paste0("p",i),arr.ind=TRUE)[1,]
    if (diff(matind)!=0) ## if non-diagonal
      transformlist[i] <- "covlogit"
  }
  values <- map_dbl(seq_along(parlocs),~do.call(transformlist[.],
                                                list(as.numeric(x3[parlocs[.]]))))
  backtransformfunc <- function(p) {
    if (missing(p)) p <- values
    res <- x2
    backfuns <- backtransfunchooser(transformlist)
    seqnow <- seq_along(p)
    pfill <- map_dbl(seqnow,~backfuns[[.]](p[.]))
    for (i in seqnow)
      res[paste0("p",i) ==matmapping] <- pfill[i]
    sdmat <- diag(sqrt(diag(res)))
    cormat <- res;diag(cormat) <- 1
    res2 <- sdmat %*% cormat %*% sdmat
    fixedvals <- grepl("fix",matmapping)
    res2[fixedvals] <- x2[fixedvals]
    res2
  }
  d_psi_d_psitrans_long <- function(p) {
    if (missing(p)) p <- values
    actvals <- backtransformfunc(p)
    ## then just return by map(p) actvals(p) for diagonals, and something a bit more elaborate for offdiags!
    res0 <- matrix(0,nrow(x2),ncol(x2))
    map(seq_along(p),function(i) {
      res <- res0
      indices <- which(matmapping==paste0("p",i),arr.ind=TRUE)
      if (transformlist[i]=="log") {
        res[indices] <- actvals[indices]
      } else {
        corr <- 2*exp(p[i])/(exp(p[i]) + 1)^2
        res[indices] <- corr*sqrt(actvals[indices[1,1],indices[1,1]]*actvals[indices[1,2],indices[1,2]])
      }
      res
    })
  }
  derivseq <- derivtransfunchooser(transformlist)
  d_psi_d_psitrans_short <- function(p) {
    if (missing(p)) p <- values
    d_psi_d_psitrans_long(p) %>% map(~unique(.[.!=0]))
  }
  d_bi_d_omega <- function(biseqt,pt) {
    if (missing(pt)) pt <- values
    m <- backtransformfunc(pt)
    map(seq_along(values),function(i) {
      pind <- paste0("p",i)
      if (pind %in% diag(matmapping)) {
        x <- rep(0,nrow(m))
        x[diag(matmapping)==pind] <- m[matmapping==pind]
        ifelse(x==0,0,biseqt * 1/sqrt(x)/2)
      } else {
        # gen_bi_here <- function(p) {
        #   ptp <- pt
        #   ptp[i] <- ptp[i]+p/exp(ptp[i]) ## get derivatives on normal scale
        #   Omega <- backtransformfunc(ptp)
        #   if (any(diag(Omega)==0))
        #     diag(Omega)[diag(Omega)==0] <- 1e-10
        #   bi %*% chol(Omega)
        # }
        gen_bi_here <- function(p) {
          Omega <- m
          Omega[matmapping==pind] <- Omega[matmapping==pind]+p
          if (any(diag(Omega)==0))
            diag(Omega)[diag(Omega)==0] <- 1e-10
          biseqt %*% chol(Omega)
        }
        c(jacobiann(gen_bi_here,0))
      }
    })
  }
  mm <- matrix(0,nrow(x2),ncol(x2))
  d_omega_d_Omega <- map(paste0("p",seq_along(values)),function(i) {
    mm[matmapping==i] <- 1
    mm})
  list(values=values,backtransformfunc=backtransformfunc,transfunclist=transformlist,mapping=matmapping,
       d_psi_d_psitrans_long=d_psi_d_psitrans_long,d_psi_d_psitrans_short=d_psi_d_psitrans_short,
       d_bi_d_omega=d_bi_d_omega,d_omega_d_Omega=d_omega_d_Omega)
}
