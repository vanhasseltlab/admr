# Step 1: Create a wide format dataset by filtering out EVID == 101
examplomycin_wide <- examplomycin %>%
  filter(EVID != 101) %>%
  dplyr::select(ID, TIME, DV) %>%
  pivot_wider(names_from = TIME, values_from = DV) %>%
  dplyr::select(-c(1))  # Remove the ID column

# Step 2: Add COV1 to the data
examplomycin_wide_with_cov1 <- cbind(examplomycin_wide, COV1 = examplomycin$COV1)

# Step 3: Apply `meancov()` function from the `admr` package
examplomycin_aggregated <- admr::meancov(examplomycin_wide)
examplomycin_aggregated <- admr::meancov(examplomycin_wide_with_cov1)
head(examplomycin_aggregated)


rxModel <- RxODE({
  cl <- cl * (COV1 / 70)^0.75
  cp = linCmt(    # Solved one- or two-compartment model
    cl,           # Clearance
    v1,           # Volume of the central compartment
    v2,           # Volume of the peripheral compartment
    q,            # Inter-compartmental clearance
    ka            # Absorption rate constant
  )
})

predder <- function(time, theta_i, dose = 100) {
  n_individuals <- nrow(theta_i)

  if (is.null(n_individuals)) {
    n_individuals <- 1
  }

  cov1_values <- as.data.frame(rnorm(n_individuals, mean = 70, sd = 15))

  # Create the event table for dosing and sampling
  ev <- eventTable(amount.units="mg", time.units="hours")
  ev$add.dosing(dose = dose, nbr.doses = 1, start.time = 0)  # Assume same dose
  ev$add.sampling(time)  # Time points to simulate

  # Solve the RxODE model for all individuals
  out <- rxSolve(rxModel, params = theta_i, iCov = cov1_values, events = ev, cores = 0)

  # Extract the predicted concentrations (cp) and format them as a matrix
  cp_matrix <- matrix(out$cp, nrow = n_individuals, ncol = length(time),
                      byrow = TRUE)

  # Return matrix of predicted concentrations
  return(cp_matrix)
}

opts <- genopts(time=c(.1,.25,.5,1,2,3,5,8,12),
                  p=list(beta=c(cl = 5, v1 = 10, v2 = 30,q = 10, ka = 1, COV1eff = 0.75),Omega=omegas(.09,0,5),Sigma_prop=0.01),
                  nsim=1000,n=1000,
                  adist=NULL, interact=TRUE,fo_appr=FALSE,biseq=NA,omega_expansion=1,
                  f=predder)

#fit <- timedbobyqa_cov(opts$pt, opts, examplomycin_aggregated)


predder(c(.1,.25,.5,1,2,3,5,8,12), c(cl = 5, v1 = 10, v2 = 30,q = 10, ka = 1, COV1eff = 0.75), dose = 100)
