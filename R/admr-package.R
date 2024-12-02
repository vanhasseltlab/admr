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
#' @importFrom stats coef cov.wt lm D qnorm cov2cor na.omit rnorm weights optimize
#' @importFrom utils relist tail head
#' @importFrom purrr map map2 map_dbl map_lgl
#' @importFrom stringr str_extract
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom mnorm dmnorm
#' @importFrom randtoolbox sobol
#' @importFrom calculus gradient jacobian
#' @importFrom stringr str_detect
#' @importFrom readr parse_number
#' @importFrom nloptr nloptr
#' @importFrom optimx optimx
#' @importFrom numDeriv hessian
## usethis namespace: end
NULL
