#' @noRd
anticovlogit <- function(x) { # Performs logit and inverse logit transformations for values constrained to the [-1,1] domain.
  2*exp(x)/(exp(x)+1)-1
}

#' @noRd
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

#' @noRd
covlogit <- function(x) { # Performs logit and inverse logit transformations for values constrained to the [-1,1] domain.
  x2 <- (x+1)/2
  log(x2/(1-x2))
}

#' @noRd
derivtransfunchooser <- function(fs) { # Returns appropriate back-transformation or derivative functions based on a list of transformation names.
  translist <- c(log=exp,
                 covlogit=function(x) 2*exp(x)/(1+exp(x))^2,
                 logit=function(x) exp(x)/(1+exp(x))^2,
                 identity=function(x) 1)
  res <- translist[fs]
  for (i in which(map_lgl(res,is.null))) {
    names(res)[i] <- fs[i]
    res[[i]] <- function(x) grad(fs[i],x)
  }
  res
  }

#' @noRd
foceapprEV_single <- function(opts, bi, biseq) {
  if (missing(bi)) bi <- matrix(0, 1, ncol(opts$p$Omega))
  theta_i <- g_iter(opts, bi)
  Eind <- opts$f(opts$time, theta_i)
  Etyp <- opts$f(opts$time, opts$p$beta)

  d_f_d_bi <- jacobiann_vec_fast_cpp(function(et) opts$f(opts$time, g_iter(opts, et)), bi)

  E <- t(c(Eind) - d_f_d_bi %*% t(bi))
  V <- compute_variance_cpp(d_f_d_bi, opts$p$Omega)

  if (opts$interact) {
    V <- opts$h(list(V = V, E = Eind), opts$p)$V
  } else {
    V <- opts$h(list(V = V, E = Etyp), opts$p)$V
  }

  list(E = E, V = V, Etyp = Etyp, Eind = Eind, d_f_d_bi = d_f_d_bi)
}


#' @noRd
grad <- function(f,var) gradient(f,unname(var)) # Gradient of a function without names.

#' @noRd
g_iter <- function(opts, bi) {
  if (!is.matrix(bi)) {
    opts$g(opts$p$beta, bi, opts$ai[1])
  } else {
    g_iter_generic_cpp(opts$g, opts$p$beta, bi, opts$ai[, 1])
  }
}

#' @noRd
gen_bi <- function(opts,expanded=TRUE) { # Generates random effects bi based on the covariance matrix in the options, with optional expansion.
  Omega <- opts$p$Omega
  if (expanded) Omega <- Omega*opts$omega_expansion
  gen_bi2(Omega,opts$biseq)
}


#' @noRd
gen_bi2 <- function(Omega,biseq) { # Generates random effects bi based on the covariance matrix in the options, with optional expansion.
  ## error handling for the  special case of zero omegas
  gen_bi2_cpp(Omega, biseq)
}


#' @noRd
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

#' @noRd
logdens2wt <- function(x) { # Converts log-densities into normalized weights for importance sampling.
  # Use C++ version for better performance
  as.vector(logdens2wt_cpp(x))
}

#' @noRd
lowtoup <- function(m) { # Copies the upper triangular part of a matrix to the lower triangular
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}

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

#' @noRd
maxfunc_old <- function(opts) {
  bi <- gen_bi(opts)
  origbeta <- opts$p$beta
  rawpreds <- opts$f(opts$time,g_iter(opts,bi))
  NAfilters <- apply(rawpreds,1,function(i) !any(is.na(i)))
  adjust <- !any(is.na(opts$single_betas))
  #adjust = F
  propdens <- mnorm::dmnorm(bi,mean=rep(0,nrow(opts$p$Omega)),sigma=opts$p$Omega*opts$omega_expansion,log=TRUE)$den

  function(pnew,returnEV=FALSE) {
    pneww <- opts$ptrans(pnew)
    newdens <- opts$p_thetai(pneww,origbeta,bi)
    wttot <- logdens2wt(newdens[NAfilters]-propdens[NAfilters])
    EVnow <- with(cov.wt(rawpreds[NAfilters,],wttot,method="ML"),list(E=center,V=cov))
    EVnow <- opts$h(EVnow,pneww)
    if (adjust) {
      kappa <- opts$f(opts$time, t(as.matrix(opts$g(ifelse(opts$single_betas, pneww$beta, origbeta))))) - rawpreds[1,]
      EVnow$E <- EVnow$E + kappa
    }

    if (returnEV) {
      EVnow
    } else if (opts$no_cov) {
      nllfun_var(opts$obs, EVnow, n = opts$n)
    } else {
      nllfun(opts$obs, EVnow, n = opts$n)
    }
  }
}

#' @noRd
maxfunc <- function(opts) {
  bi <- gen_bi(opts)
  origbeta <- opts$p$beta
  rawpreds <- opts$f(opts$time,g_iter(opts,bi))
  NAfilters <- apply(rawpreds,1,function(i) !any(is.na(i)))
  adjust <- !any(is.na(opts$single_betas))
  #adjust = F
  propdens <- mnorm::dmnorm(bi,mean=rep(0,nrow(opts$p$Omega)),sigma=opts$p$Omega*opts$omega_expansion,log=TRUE)$den

  function(pnew) {
    pneww <- opts$ptrans(pnew)
    newdens <- opts$p_thetai(pneww,origbeta,bi)
    wttot <- logdens2wt(newdens[NAfilters]-propdens[NAfilters])
    EVnow <- with(cov.wt(rawpreds[NAfilters,],wttot,method="ML"),list(E=center,V=cov))
    EVnow <- opts$h(EVnow,pneww)
    if (adjust) {
      kappa <- opts$f(opts$time, t(as.matrix(opts$g(ifelse(opts$single_betas, pneww$beta, origbeta))))) - rawpreds[1,]
      EVnow$E <- EVnow$E + kappa
    }
    if  (opts$no_cov) {
      nllfun_var(opts$obs, EVnow, n = opts$n)
    } else {
      nllfun(opts$obs, EVnow, n = opts$n)
    }
  }
}

#' @noRd
MCapprEV <- function(opts) {
  p <- opts$p
  ## p is normal-scale parameters
  bi <- gen_bi(opts)
  theta_i <- g_iter(opts, bi)
  m <- opts$f(opts$time, theta_i)

  if (opts$omega_expansion == 1) {
    # Uniform weights if no expansion
    w <- rep(1, nrow(m))
    r <- meancov_cpp(m, w)
  } else {
    # Weighted case
    wt <- samplogdensfun_cpp(bi, p, opts$omega_expansion)
    r <- meancov_cpp(m, logdens2wt(wt))
  }

  opts$h(r, opts$p)
}



#' @noRd
nllfun <- function(obsEV,predEV,invpredV,n=1) {   # Computes the negative log-likelihood given observed and predicted expected values and their covariance.
  # Use C++ version for better performance
  nllfun_cpp(obsEV$E, obsEV$V, predEV$E, predEV$V, n)
}

#' @noRd
nllfun_var <- function(obsEV,predEV,n=1) {   # Computes the negative log-likelihood given observed and predicted expected values and their covariance.
  # Use C++ version for better performance
  nllfun_var_cpp(obsEV$E, diag(obsEV$V), predEV$E, diag(predEV$V), n)
}


#' @noRd
numtransfun <- function(origfun) # Creates a back-transformation function numerically based on an original transformation function.
  ## Numeric back-transformation function, inefficient.
  function(x) optimize(function(xx) (xx-do.call(origfun,list(x)))^2,c(-1000,1000))$minimum
translist_string <- c(log="exp",
                      covlogit="2*exp(x)/(exp(x)+1)-1",
                      logit="exp(x)/(1+exp(x))",
                      identity="x")

#' @noRd
obs2opts <- function(opts,obs) {opts$obs <- obs;opts} # Adds observed data to the model options.

#' @noRd
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
        gen_bi_here <- function(p) {
          Omega <- m
          Omega[matmapping==pind] <- Omega[matmapping==pind]+p
          if (any(diag(Omega)==0))
            diag(Omega)[diag(Omega)==0] <- 1e-10
          biseqt %*% chol(Omega)
        }
        c(jacobiann_vec_fast_cpp(gen_bi_here, 0))
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

#' @noRd
upd_opts <- function(opts0,optslist) { # Initializes and generates core options and settings for modeling and optimization, including random effects, simulation settings, and likelihood approximations.
  ## note: this function always resets biseq and adist, and takes away obs!!
  if (!all(names(optslist) %in% names(formals(genopts))))
    stop("Argument not understood by genopts() function!")
  update_biseq <- TRUE
  if ("biseq" %in% names(optslist)) update_biseq <- FALSE
  for (i in names(optslist)) opts0[[i]] <- optslist[[i]]
  if (update_biseq) opts0$biseq <- NA

  opts0$adist <- NULL
  opts0$ptrans <- NULL
  opts0$pderiv <- NULL
  opts0$ai <- NULL
  opts0$pt <- NULL
  opts0$obs <- NULL
  opts0$d_g_d_beta <- NULL
  opts0$d_g_d_bi <- NULL
  opts0$d_bi_d_omega <- NULL
  opts0$d_omega_d_Omega <- NULL
  do.call(genopts,opts0)
}

#' @noRd
uptolow <- function(m) { # Copies the upper triangular part of a matrix to the lower triangular and vice versa.
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

#' @noRd
calculate_bsv <- function(back_transformed_params) {
  tryCatch({
    # Extract diagonal elements from the Omega matrix
    omega_diag <- diag(back_transformed_params$Omega)

    # Calculate BSV (CV%) for (log-)normal parameters
    #bsv_cv <- sqrt(exp(omega_diag) - 1) * 100
    bsv_cv <- sqrt(omega_diag) * 100

    # Handle potential NA cases or missing values in Omega
    c(bsv_cv, rep(NA, length(omega_diag) - length(bsv_cv)))
  }, error = function(e) {
    # In case of any error, return NA for all
    rep(NA, length(diag(back_transformed_params$Omega)))
  })
}

#' @noRd
perturb_init <- function(init, perturbation) {
  init + rnorm(length(init), mean = 0, sd = perturbation * abs(init))  # Slight perturbation
}

#' @noRd
compute_nll <- function(opts, params, nomap) {
  if (nomap) {
    maxfunc(p2opts(opts, params))(params)
  } else {
    Reduce('+', map(opts, ~ maxfunc(p2opts(., params))(params)))
  }
}

#' @noRd
run_chain <- function(chain, opts, obs, init, maxiter, phase_fractions, convcrit_nll,
                      max_worse_iterations, perturbation, nomap, use_grad) {
  chain_start_time <- Sys.time()  # Start time for the chain
  chain_init <- if (chain == 1) init else perturb_init(init, perturbation)

  if (nomap) {
    opts <- opts %>% p2opts(chain_init) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts), function(i) opts[[i]] %>% p2opts(chain_init) %>% obs2opts(obs[[i]]))
  }

  res <- tibble(
    iter = 1,
    parameters = vector("list", maxiter),
    nll = NA_real_,
    approx_nll = NA_real_,
    iteration_time = 0
  )

  res$parameters[[1]] <- chain_init

  if (nomap) {
    res$nll[1] <- maxfunc(opts)(chain_init)
  } else {
    res$nll[1] <- Reduce('+', map(opts, ~ maxfunc(.)(chain_init)))
  }

  res$approx_nll[1] <- res$nll[1]

  if (chain == 1) {  # Print live output for single-chain mode
    cat(sprintf("Chain %d:\nIter | %-12s\n", chain, sprintf("NLL and Parameters (%d values)", length(init))))
    cat(sprintf("%s\n", strrep("-", 80)))
    cat(sprintf("%4d: %s\n", 1, paste(sprintf("%8.3f", c(res$nll[1], chain_init)), collapse = " ")))
  }

  phase_params <- list(
    list(bounds = 2, maxeval = 5000, phase_name = "Wide Search Phase"),
    list(bounds = 1, maxeval = 5000, phase_name = "Focussed Search Phase"),
    list(bounds = 0.5, maxeval = 5000, phase_name = "Fine-Tuning Phase"),
    list(bounds = 0.01, maxeval = 5000, phase_name = "Precision Phase")
  )

  if (length(phase_fractions) != length(phase_params)) {
    stop("Number of phase_fractions must match the number of phases.")
  }

  cumulative_iters <- cumsum(c(floor(phase_fractions[-length(phase_fractions)] * maxiter), maxiter))
  phase_start_iters <- c(2, head(cumulative_iters, -1) + 1)
  phase_end_iters <- cumulative_iters

  best_nll <- res$nll[1]
  best_params <- res$parameters[[1]]
  current_iter <- 2

  for (phase_idx in seq_along(phase_params)) {
    current_phase <- phase_params[[phase_idx]]

    if (chain == 1) {
      cat(sprintf("\n### %s ###\n", current_phase$phase_name))
    }

    # Call handle_phase
    phase_result <- handle_phase(
      current_phase = current_phase,
      opts = opts,
      chain_init = best_params,
      best_params = best_params,
      best_nll = best_nll,
      current_iter = current_iter,
      phase_end_iter = phase_end_iters[phase_idx],
      nomap = nomap,
      convcrit_nll = convcrit_nll,
      max_worse_iterations = max_worse_iterations,
      res = res,
      maxiter = maxiter,
      chain = chain,  # Pass the chain index for printing control
      use_grad = use_grad
    )

    # Update variables for the next phase
    best_params <- phase_result$best_params
    best_nll <- phase_result$best_nll
    res <- phase_result$res
    current_iter <- phase_result$current_iter

    if (phase_result$phase_converged) {
      cat(sprintf("Phase %s converged at iteration %d.\n", current_phase$phase_name, current_iter - 1))
    }
  }

  res <- res[!is.na(res$nll), ]
  chain_time <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))

  cat(sprintf("\nChain %d Complete: Final NLL = %.3f, Time Elapsed = %.2f seconds\n \n",
              chain, best_nll, chain_time))

  return(list(res = res, best_nll = best_nll, best_params = best_params, time = chain_time))
}

#' @noRd
run_chainMC_old <- function(chain, opts, obs, init, maxiter, convcrit_nll, perturbation, nomap) {
  chain_start_time <- Sys.time()  # Start time for the chain
  chain_init <- if (chain == 1) init else perturb_init(init, perturbation)

  if (nomap) {
    opts <- opts %>% p2opts(chain_init) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts), function(i) opts[[i]] %>% p2opts(chain_init) %>% obs2opts(obs[[i]]))
  }

  res <- tibble(
    iter = 1,
    parameters = vector("list", maxiter),
    nll = NA_real_,
    iteration_time = 0
  )

  i <- 1  # Initialize iteration counter for fitfun2


  # Generate objective function based on number of models
  if (nomap) {
    fitfun <- genfitfunc(opts, obs)
  } else {
    fitfuns <- map(seq_along(opts), ~ genfitfunc(opts[[.x]], obs[[.x]]))
    fitfun <- function(p) Reduce('+', map(fitfuns, ~ .(p)))
  }

  # Create wrapper function to track optimization progress
  fitfun2 <- function(p) {
    nllNow <- fitfun(p)
    res$parameters[[i]] <<- p
    res$iteration_time[i] <<- Sys.time()
    res$nll[i] <<- nllNow
    res$iter[i] <<- i

    # Print progress every 50 iterations
    if (i %% 50 == 0) {
      cat("Iteration:", i, "- NLL:", nllNow, "\n")
    }
    i <<- i + 1
    nllNow
  }

  # Run BOBYQA optimization with bounds
  est <- nloptr::nloptr(chain_init,
                        fitfun2,
                        lb=chain_init-2,  # Lower bounds: 2 units below initial values
                        ub=chain_init+2,  # Upper bounds: 2 units above initial values
                        opts=list(
                          algorithm="NLOPT_LN_BOBYQA",
                          ftol_rel=.Machine$double.eps^2,  # Relative function tolerance
                          maxeval = maxiter-2,                    # max number of evaluations
                          check_derivatives = F))  # Skip derivative checking for performance


  # Update variables for the next phase
  best_params <- est$solution
  best_nll <- fitfun(best_params)
  res <- res[!is.na(res$nll), ]
  chain_time <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))

  cat(sprintf("\nChain %d Complete: Final NLL = %.3f, Time Elapsed = %.2f seconds\n \n",
              chain, best_nll, chain_time))

  return(list(res = res, best_nll = best_nll, best_params = best_params, time = chain_time))
}

run_chainMC <- function(chain, opts, obs, init, maxiter, convcrit_nll, perturbation, nomap, use_grad) {

  chain_start_time <- Sys.time()  # Start time for the chain
  chain_init <- if (chain == 1) init else perturb_init(init, perturbation)
  n_params <- length(chain_init)

  # Precompute opts for this chain

  if (nomap) {
    opts_chain <- opts %>% p2opts(chain_init) %>% obs2opts(obs)
  } else {
    opts_chain <- map(seq_along(opts), function(i) opts[[i]] %>% p2opts(chain_init) %>% obs2opts(obs[[i]]))
  }

  # Precompute fit function

  if (nomap) {
    fitfun <- genfitfunc(opts_chain, obs)
  } else {
    fitfuns <- map(seq_along(opts_chain), ~ genfitfunc(opts_chain[[.x]], obs[[.x]]))
    fitfun <- function(p) Reduce('+', map(fitfuns, ~ .(p)))
  }

  # Gradient factory function



  if (use_grad) {
    gradfun <- make_gradfun_mc(fitfun)
  } else {
    gradfun <- NULL
  }

  # Preallocate storage

  res_params <- matrix(NA_real_, nrow = maxiter, ncol = n_params)
  res_nll <- numeric(maxiter)
  res_time <- numeric(maxiter)

  i <- 1
  iter_start_time <- Sys.time()

  # Wrapper to track optimization progress

  fitfun2 <- function(p) {
    nllNow <- fitfun(p)
    res_params[i, ] <<- p
    res_nll[i] <<- nllNow
    res_time[i] <<- as.numeric(difftime(Sys.time(), iter_start_time, units = "secs"))
    print_freq <- if (use_grad) 5 else 50
    if (i %% print_freq == 0) {
      cat("Iteration:", i, "- NLL:", nllNow, "\n")
    }
    i <<- i + 1
    nllNow
  }

  # Choose algorithm

  algo <- if (use_grad) "NLOPT_LD_LBFGS" else "NLOPT_LN_BOBYQA"

  # Run optimization

  est <- nloptr::nloptr(
    x0 = chain_init,
    eval_f = fitfun2,
    eval_grad_f = if (use_grad) gradfun else NULL,
    lb = chain_init - 2,
    ub = chain_init + 2,
    opts = list(
      algorithm = algo,
      ftol_rel = .Machine$double.eps^2,
      maxeval = maxiter,
      check_derivatives = FALSE
    )
  )

  best_params <- est$solution
  best_nll <- fitfun(best_params)

  # Build tibble result

  res <- tibble(
    iter = seq_len(i-1),
    parameters = split(res_params[1:(i-1), , drop = FALSE], seq_len(i-1)),
    nll = res_nll[1:(i-1)],
    iteration_time = res_time[1:(i-1)]
  )

  chain_time <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))

  cat(sprintf("\nChain %d Complete: Final NLL = %.3f, Time Elapsed = %.2f seconds\n\n",
              chain, best_nll, chain_time))

  return(list(res = res, best_nll = best_nll, best_params = best_params, time = chain_time))
}


#' @noRd
gradfun <- function(p) {
  eps <- 1e-6
  f0 <- fitfun(p)
  g <- numeric(length(p))
  for (k in seq_along(p)) {
    p_eps <- p
    p_eps[k] <- p_eps[k] + eps
    g[k] <- (fitfun(p_eps) - f0) / eps
  }
  g
}


#' @noRd
handle_phase_old <- function(current_phase, opts, chain_init, best_params, best_nll, current_iter, phase_end_iter,
                         nomap, convcrit_nll, max_worse_iterations, res, maxiter, chain) {
  phase_converged <- FALSE
  worse_counter <- 0

  for (i in current_iter:phase_end_iter) {
    if (current_iter > maxiter) {
      cat("Maximum iterations reached. Terminating optimization.\n")
      break
    }

    iter_start_time <- Sys.time()

    if (nomap) {
      ff <- maxfunc(p2opts(opts, chain_init))
    } else {
      ffs <- map(opts, ~ maxfunc(p2opts(., chain_init)))
      ff <- function(p) Reduce('+', map(ffs, ~ .(p)))
    }

    m0 <- nloptr::nloptr(
      x0 = chain_init,
      eval_f = function(params) ff(params),
      lb = chain_init - current_phase$bounds,
      ub = chain_init + current_phase$bounds,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",
        ftol_rel = .Machine$double.eps,
        maxeval = current_phase$maxeval
      )
    )

    # Optimization step
    chain_init <- m0$solution
    res$parameters[[current_iter]] <- chain_init
    res$nll[current_iter] <- compute_nll(opts, chain_init, nomap)
    res$approx_nll[current_iter] <- m0$objective
    res$iteration_time[current_iter] <- as.numeric(difftime(Sys.time(), iter_start_time, units = "secs"))

    # Update best parameters if a better NLL is found
    if (res$nll[current_iter] < best_nll) {
      best_nll <- res$nll[current_iter]
      best_params <- chain_init
    }

    # Printing live updates for the chain
    if (chain == 1) {  # Only print for the first chain
      cat(sprintf("%4d: %s\n", current_iter, paste(sprintf("%8.3f", c(res$nll[current_iter], chain_init)), collapse = " ")))
    }

    # Convergence check
    if (abs(res$nll[current_iter] - res$approx_nll[current_iter]) < convcrit_nll) {
      phase_converged <- TRUE
      current_iter <- current_iter + 1
      break  # Skip the rest of this phase, continue to the next phase
    }

    # Early phase termination due to worse counter
    worse_counter <- if (res$nll[current_iter] > best_nll) worse_counter + 1 else 0
    if (worse_counter >= max_worse_iterations) {
      current_iter <- current_iter + 1
      break  # Skip the rest of this phase, continue to the next phase
    }

    current_iter <- current_iter + 1  # Increment after every iteration
  }

  return(list(
    best_params = best_params,
    best_nll = best_nll,
    res = res,
    phase_converged = phase_converged,
    current_iter = current_iter
  ))
}

handle_phase <- function(current_phase, opts, chain_init, best_params, best_nll, current_iter, phase_end_iter,
                         nomap, convcrit_nll, max_worse_iterations, res, maxiter, chain,
                         use_grad) {

  phase_converged <- FALSE
  worse_counter <- 0

  for (i in current_iter:phase_end_iter) {
    if (current_iter > maxiter) {
      cat("Maximum iterations reached. Terminating optimization.\n")
      break
    }

    iter_start_time <- Sys.time()

    # Build objective function
    if (nomap) {
      ff <- maxfunc(p2opts(opts, chain_init))
    } else {
      ffs <- map(opts, ~ maxfunc(p2opts(., chain_init)))
      ff <- function(p) Reduce('+', map(ffs, ~ .(p)))
    }

    # Only define gradient once
    if (use_grad) gradfun <- make_gradfun_irmc(ff)

    # Choose algorithm
    algo <- if (use_grad) "NLOPT_LD_LBFGS" else "NLOPT_LN_BOBYQA"

    m0 <- nloptr::nloptr(
      x0 = chain_init,
      eval_f = ff,
      eval_grad_f = if (use_grad) gradfun else NULL,
      lb = chain_init - current_phase$bounds,
      ub = chain_init + current_phase$bounds,
      opts = list(
        algorithm = algo,
        ftol_rel = .Machine$double.eps,
        maxeval = current_phase$maxeval,
        check_derivatives = FALSE
      )
    )

    # Optimization step
    chain_init <- m0$solution
    res$parameters[[current_iter]] <- chain_init
    res$nll[current_iter] <- compute_nll(opts, chain_init, nomap)
    res$approx_nll[current_iter] <- m0$objective
    res$iteration_time[current_iter] <- as.numeric(difftime(Sys.time(), iter_start_time, units = "secs"))

    # Update best parameters if a better NLL is found
    if (res$nll[current_iter] < best_nll) {
      best_nll <- res$nll[current_iter]
      best_params <- chain_init
    }

    # Printing live updates (same as original)
    if (chain == 1) {
      cat(sprintf("%4d: %s\n", current_iter,
                  paste(sprintf("%8.3f", c(res$nll[current_iter], chain_init)), collapse = " ")))
    }

    # Convergence check
    if (abs(res$nll[current_iter] - res$approx_nll[current_iter]) < convcrit_nll) {
      phase_converged <- TRUE
      current_iter <- current_iter + 1
      break
    }

    # Early termination due to worse counter
    worse_counter <- if (res$nll[current_iter] > best_nll) worse_counter + 1 else 0
    if (worse_counter >= max_worse_iterations) {
      current_iter <- current_iter + 1
      break
    }

    current_iter <- current_iter + 1
  }

  return(list(
    best_params = best_params,
    best_nll = best_nll,
    res = res,
    phase_converged = phase_converged,
    current_iter = current_iter
  ))
}

#' @noRd
prepare_boxplot_df <- function(E, V, time_points) { # Function to convert E and V to plot-ready data frame
  mu <- as.numeric(E)
  sd <- sqrt(diag(V))

  data.frame(
    time = time_points,
    mean = mu,
    lower_q1 = mu - 0.674 * sd,
    upper_q3 = mu + 0.674 * sd,
    lower_95 = mu - 1.96 * sd,
    upper_95 = mu + 1.96 * sd
  )
}

#' @noRd
prepare_E_df <- function(E, label, time_points) {
  data.frame(
    time = time_points,
    mean = as.numeric(E),
    source = label
  )
}

#' @noRd
prepare_V_df <- function(V, label, time_points) {
  V_mat <- as.matrix(V)

  # Extract time points from row/col names
  time_points_row <- rownames(V_mat)
  time_points_col <- colnames(V_mat)
  expand.grid(
    time1 = time_points,
    time2 = time_points
  ) |>
    transform(value = as.vector(V_mat), source = label)
}

#' @noRd
make_gradfun_irmc <- function(ff, eps = 1e-6) {
  function(p) {
    f0 <- ff(p)
    g <- numeric(length(p))
    for (k in seq_along(p)) {
      p_eps <- p
      p_eps[k] <- p_eps[k] + eps
      g[k] <- (ff(p_eps) - f0) / eps
    }
    g
  }
}

#' @noRd
make_gradfun_mc <- function(fitfun, eps = 1e-6) {
  function(p) {
    f0 <- fitfun(p)
    g <- numeric(length(p))
    for (k in seq_along(p)) {
      p_eps <- p
      p_eps[k] <- p_eps[k] + eps
      g[k] <- (fitfun(p_eps) - f0) / eps
    }
    g
  }
}
