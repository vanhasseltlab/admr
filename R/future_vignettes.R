library(rxode2)
library(nlmixr2)
library(dplyr)
library(tidyr)
library(mnorm)
library(MASS)
library(data.table)
library(tidyverse)

generate_data <- function(n, times, seed = 1) {
  set.seed(seed)

  mod <- rxode({
    # Parameters
    ke = cl / v1             # Elimination rate constant
    k12 = q / v1             # Rate constant for central to peripheral transfer
    k21 = q / v2             # Rate constant for peripheral to central transfer

    # Differential equations for drug amount in compartments
    d/dt(depot)    = -ka * depot
    d/dt(central)  = ka * depot - ke * central -
      k12 * central +
      k21 * peripheral
    d/dt(peripheral) = k12 * central -
      k21 * peripheral

    # Concentration in the central compartment
    cp = central / v1
  })

  # Generate covariates
  theta <- c(cl=5, v1 = 10, v2 =30, q =10, ka = 1)

  omegaCor <- matrix(c(1,  0,  0,  0,  0,
                       0,  1,  0,  0,  0,
                       0,  0,  1,  0,  0,
                       0,  0,  0,  1,  0,
                       0,  0,  0,  0,  1),
                     dimnames=list(NULL,c("eta.cl","eta.v1","eta.v2", "eta.q",
                                          "eta.ka")), nrow=5)

  iiv.sd <- c(0.3, 0.3, 0.3, 0.3, 0.3) ## SDs of model parameters

  iiv <- iiv.sd %*% t(iiv.sd)
  omega <- iiv * omegaCor  # covariance matrix

  mv <- mvrnorm(n, rep(0, dim(omega)[1]), omega)

  params.all <-
    data.table(
      "ID" = seq(1:n),
      "cl" = theta['cl'] * exp(mv[, 1]),
      "v1" = theta['v1'] * exp(mv[, 2]),
      "v2" = theta['v2'] * exp(mv[, 3]),
      "q"  = theta['q']  * exp(mv[, 4]),
      "ka" = theta['ka'] * exp(mv[, 5]),
      "COV1" = round(rnorm(n, 70, 15))
    )

  params.all$cl <- params.all$cl * (params.all$COV1 / 70)^0.75

  # Event table
  ev <- et() %>%
    et(amt = 100) %>%  # Add single dose
    et(0) %>%  # Add initial time observation
    et(times) %>%  # Sampling schedule
    et(ID = seq(1, n)) %>%  # Assign unique IDs for all subjects
    as.data.frame()

  # Solve the model
  sim <- rxSolve(mod, events = ev, iCov = params.all, cores = 0, addCov = T) %>%
    mutate(ID = as.integer(id), TIME = as.numeric(time)) %>%
    dplyr::select(-c(id, time)) %>%
    mutate(AMT =  ifelse(TIME == 0, 100, 0)) %>%
    mutate(EVID = ifelse(TIME == 0, 101, 0)) %>%
    mutate(CMT = ifelse(TIME == 0, 1, 2))

  # Simulate residual error and combine
  sim$rv <- rnorm(nrow(sim), 0, 0.2)
  sim$DV <- round(sim$cp * (1 + sim$rv), 3)
  sim <- merge(sim, params.all)

  dat <- sim %>%
    dplyr::select("ID", "TIME", "DV", "AMT", "EVID", "CMT", COV1)

  return(dat)
}

examplomycin <- generate_data(n = 1000, times = c(.1,.25,.5,1,2,3,5,8,12),
                              seed = 1)
head(examplomycin)
