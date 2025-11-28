#' @docType package
#' @name admr-package
#' @aliases admr-package admr
#' @title admr: Aggregate data modelling in R
#' @description
#' A novel method to aggregate data and model it using a non-linear estimation techniques.
#' @author \itemize{
#' **Maintainer:** H. van de Beek \email{h.van.de.beek@lacdr.leidenuniv.nl}
#'
#' Authors:
#'
#' \item H. van de Beek
#' \item P.A.J. Välitalo}
#' @keywords internal
"_PACKAGE"
## usethis namespace: start
#' @import ggplot2
#' @importFrom calculus gradient jacobian
#' @importFrom magrittr %>%
#' @importFrom mnorm dmnorm
#' @importFrom nloptr nloptr
#' @importFrom numDeriv hessian
#' @importFrom optimx optimx
#' @importFrom purrr map map2 map_dbl map_lgl
#' @importFrom randtoolbox sobol
#' @importFrom Rcpp sourceCpp
#' @importFrom readr parse_number
#' @importFrom stats coef cov.wt lm D qnorm cov2cor na.omit rnorm weights optimize median
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom tibble tibble
#' @importFrom utils relist tail head
#' @useDynLib admr, .registration = TRUE
## usethis namespace: end
NULL
