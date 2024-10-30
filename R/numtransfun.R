#' @noRd

numtransfun <- function(origfun) # Creates a back-transformation function numerically based on an original transformation function.
  ## Numeric back-transformation function, inefficient.
  function(x) optimize(function(xx) (xx-do.call(origfun,list(x)))^2,c(-1000,1000))$minimum
translist_string <- c(log="exp",
                      covlogit="2*exp(x)/(exp(x)+1)-1",
                      logit="exp(x)/(1+exp(x))",
                      identity="x")
