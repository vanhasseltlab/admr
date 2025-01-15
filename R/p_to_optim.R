#' Parameters to optimizable form
#'
#'[p_to_optim()] Converts a parameter list into a form that can be optimized, returning transformed values, back-transformation functions, and derivatives.
#'
#' @param p The parameter list
#'
#' @returns optimized parameter list
#' @export
#' @examples
#' #test

p_to_optim <- function(p) {
  ### a master function to feed parameters to optim
  if (any(names(p) == "Sigma_exp")) {
    message("Message from p_to_optim:\n Sigma_exp will be transformed to Sigma_prop!")
    if (any(names(p) == "Sigma_prop")) {
      p$Sigma_prop <- p$Sigma_prop + p$Sigma_exp
    } else {
      p$Sigma_prop <- p$Sigma_exp
    }
    p$Sigma_exp <- NULL
  }

  # Ensure names are preserved during transformations
  p2 <- map(p, ~if (is.list(.) || is.matrix(.)) {
    p_transform_mat(.)
  } else {
    transformed <- p_transform_vec(.)
    if (is.vector(.)) names(transformed[[1]]) <- names(.)  # Preserve names
    transformed
  })

  # Extract transformed values while keeping names
  p3 <- map(p2, ~{
    res <- .[[1]]
    if (!is.null(names(res))) names(res) <- names(.[[1]])  # Ensure names persist
    res
  })

  # Clean up names to remove 'beta.' prefix
  valuess <- do.call(c, p3)
  names(valuess) <- gsub("^beta\\.", "", names(valuess))  # Remove 'beta.' prefix

  backtransformfunc <- function(pp) {
    if (any(is.na(pp))) stop("Parameters contain NA values!!")
    if (missing(pp)) pp <- valuess
    p5 <- relist(pp, p3)
    p5 <- map(seq_along(p5), ~p2[[.]][[2]](p5[[.]]))
    names(p5) <- names(p)
    p5
  }

  d_psi_d_psitrans_long <- function(pp) {
    if (missing(pp)) pp <- valuess
    res <- do.call(c, map(names(p3), function(n) {
      r <- p2[[n]]$d_psi_d_psitrans_long(pp[grepl(n, names(valuess))])
      r
    }))
    names(res) <- rep(names(p3), map_dbl(p3, length))
    res
  }

  tottransfunclist <- do.call(c, map(p2, ~.$transfunclist))
  d_psi_d_psitrans_short <- function(pp) {  ## new
    if (missing(pp)) pp <- valuess
    pp <- split(unlist(pp), factor(str_extract(names(tottransfunclist), ".*[^[:digit:]]")))
    res <- map2(p2, pp, ~.x$d_psi_d_psitrans_short(.y)) %>% unlist()
    names(res) <- names(tottransfunclist)
    res
  }

  list(
    values = valuess,
    backtransformfunc = backtransformfunc,
    d_psi_d_psitrans_long = d_psi_d_psitrans_long,
    d_psi_d_psitrans_short = d_psi_d_psitrans_short,
    d_bi_d_omega = p2$Omega$d_bi_d_omega,
    d_omega_d_Omega = p2$Omega$d_omega_d_Omega
  )
}
