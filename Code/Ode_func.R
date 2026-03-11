# ODE model functions ----------------------------------

#' Create the initial state vector for the ODE model.
#'
#' @param population: total population size
#' @param initial_seeds: the number of infected cases in the beginning
#' @param Base_Immu: the proportion of population that is initially immune
#' @param n_virus: the number of viruses in the model
#'
#' @return a named numeric vector representing the initial state of the system
Get.InitState <- function(
  population = 100000,
  initial_seeds = 10,
  Base_Immu = 0,
  n_virus = 4
) {
  state <- numeric(3 * n_virus)
  names(state) <- as.vector(rbind(
    paste0("S_", 1:n_virus),
    paste0("I_", 1:n_virus),
    paste0("R_", 1:n_virus)
  ))

  R0 <- round(population * Base_Immu)
  I0 <- initial_seeds
  S0 <- population - R0 - I0
  if (S0 < 0) {
    stop("Initial S < 0. Reduce Base_Immu or initial_seeds.")
  }

  for (i in 1:n_virus) {
    state[paste0("S_", i)] <- S0
    state[paste0("I_", i)] <- I0
    state[paste0("R_", i)] <- R0
  }
  return(state)
}


#' Post-process the raw output from the ODE simulation to extract the incidence data.
#'
#' @param Sim: the raw output from the ODE simulation (a matrix)
#' @param virus_names: optional vector of virus names to rename the incidence columns
#' @param after: only include data after this date (default "2019-01-01")
#'
#' @return a data.table containing the time and incidence columns for each virus
Model.GetI <- function(Sim, virus_names = NULL, after = as.Date("2019-01-01")) {
  dt <- Sim %>% as.data.table()
  dt[, Time := as.Date(time, origin = "1970-01-01")][,
    ISOweek := strftime(Time, "%G-W%V")
  ]
  inc_cols <- grep("^Inc_\\d+$", names(dt), value = TRUE)

  out <- dt[Time >= after, c("Time", "ISOweek", inc_cols), with = FALSE]

  if (!is.null(virus_names)) {
    stopifnot(length(virus_names) == length(inc_cols))
    setnames(out, inc_cols, virus_names)
  } else {
    setnames(out, inc_cols, paste0("Virus_", seq_along(inc_cols)))
  }
  return(out)
}


#' Run the ODE simulation with the given parameters and return the results.
#'
#' @param Parm: a list of parameters for the ODE model
#'
#' @return a list containing the raw simulation output and the processed incidence data
Model.RunSim.ode <- function(
  Parm,
  Plot = TRUE,
  virus_names = NULL,
  after = as.Date("2019-01-01")
) {
  # times <- seq(from = 1, to = 365 * Parm[["years"]])
  times <- as.numeric(seq(
    from = as.Date(Parm[["year_start"]]),
    to = as.Date(Parm[["year_end"]]),
    by = 1
  ))
  state <- Get.InitState(
    population = Parm[["num_of_agent"]],
    initial_seeds = Parm[["initial_seeds"]],
    Base_Immu = Parm[["base_immune"]],
    n_virus = length(Parm[["beta0"]])
  )

  SimResult <- ode(
    y = state,
    times = times,
    func = ParmInferenceCpp,
    parms = Parm,
    method = "rk4"
  )

  SimResult_Inc <- Model.GetI(
    SimResult,
    virus_names = virus_names,
    after = after
  )

  ResultList <- list(
    RawData = SimResult,
    Data = SimResult_Inc
  )

  if (Plot) {
    # Plotting results for each virus
    fig1 <- SimResult_Inc %>%
      as.data.frame() %>%
      pivot_longer(
        cols = !c(Time, ISOweek),
        names_to = "virus",
        values_to = "cases"
      ) %>%
      ggplot(., aes(x = Time, y = cases)) +
      geom_line() +
      scale_x_date(date_labels = "%Y-%b", date_breaks = "3 months") +
      # scale_y_continuous(limits = c(0, 10000)) +
      theme_bw() +
      facet_wrap(vars(virus), nrow = 4, scales = "free_y") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # fig2 <- SimResult_Inc %>%
    #   as.data.frame() %>%
    #   pivot_longer(cols = !Time, names_to = "virus", values_to = "cases") %>%
    #   # filter(Time > 53*4) %>%
    #   ggplot(., aes(x = Time, y = cases, colour = virus)) +
    #   geom_line(alpha = 0.7) +
    #   scale_x_continuous(breaks = seq(1, length(SimResult_Inc[, 1]), by = 52)) +
    #   # scale_y_continuous(limits = c(0, 10000)) +
    #   theme_bw()

    ResultList$fig1 <- fig1
    # ResultList$fig2 <- fig2
  }
  return(ResultList)
}


# Likelihood --------------------------------------------
Utility.Clamp01 <- function(p, epsilon = 1e-9) {
  pmin(pmax(p, epsilon), 1 - epsilon)
}

Utility.Logit <- function(p) log(p / (1 - p))

`%||%` <- function(x, y) if (is.null(x)) y else x


#' Beta-Binomial PMF log-likelihood
#' @param y: observed counts
#' @param n: total counts (population size)
#' @param p: predicted probability of infection from the model
#' @param rho: overdispersion parameter (0 < rho < 1)
#' @param epsilon: small value to avoid numerical issues with p and rho
LLH.BetaBinomial.PMF <- function(y, n, p, rho, epsilon = 1e-9) {
  p <- Utility.Clamp01(p, epsilon)
  rho <- pmin(pmax(rho, epsilon), 1 - epsilon)
  theta <- (1 - rho) / rho
  a <- p * theta
  b <- (1 - p) * theta
  # 这里用 rho in (0,1) 表示 overdispersion（类内相关）
  # theta = (1-rho)/rho, alpha = p*theta, beta = (1-p)*theta
  return(lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b))
}


#' Calculate the total log-likelihood of the observed data given the model predictions using a Beta-Binomial likelihood.
#' @param Dat: a data.table containing columns 'y' (observed counts), 'N' (population size), and 'p_sim' (predicted probabilities from the model)
#' @param rho: overdispersion parameter (0 < rho < 1)
#' @param epsilon: small value to avoid numerical issues with p and rho
LLH.BetaBinomial <- function(
  Dat,
  rho = 0.02,
  epsilon = 1e-9
) {
  return(sum(
    LLH.BetaBinomial.PMF(Dat$y, Dat$N, Dat$p_sim, rho = rho, epsilon = epsilon),
    na.rm = TRUE
  ))
}


#' Calculate the log-likelihood based on the "Ratio + Absolute Anchor" method.
#' This method uses the log-ratio of predicted vs observed probabilities for each virus relative to a reference virus, as well as the absolute logit of the reference virus.
#' @param Dat: a data.table containing columns 'Time', 'ISOweek', 'Virus', 'p_sim' (predicted probabilities), and 'p_obs' (observed probabilities)
#' @param VirNames: optional vector of virus names to include in the likelihood calculation (default is all viruses in the data)
#' @param sigma_ratio: standard deviation for the log-ratio component of the likelihood
#' @param sigma_abs: standard deviation for the absolute anchor component of the likelihood
LLH.RatioAbs <- function(
  Dat,
  VirNames = paste0("Virus_", 1:4),
  sigma_ratio = 0.7,
  sigma_abs = 0.7
) {
  DatWide <- dcast(
    Dat[, .(Time, ISOweek, Virus, p_sim, p_obs)],
    Time + ISOweek ~ Virus,
    value.var = c("p_obs", "p_sim"),
    sep = "__"
  )

  ref_obs_col <- paste0("p_obs__", "Virus_1")
  ref_sim_col <- paste0("p_sim__", "Virus_1")
  if (
    !(ref_obs_col %in% names(DatWide)) || !(ref_sim_col %in% names(DatWide))
  ) {
    return(-Inf)
  }
  Ref_Obs <- Utility.Clamp01(DatWide[[ref_obs_col]])
  Ref_Sim <- Utility.Clamp01(DatWide[[ref_sim_col]])

  ok_ref <- is.finite(Ref_Obs) & is.finite(Ref_Sim)
  if (!any(ok_ref)) {
    return(-Inf)
  }

  DatWide <- DatWide[ok_ref]
  Ref_Obs <- Ref_Obs[ok_ref]
  Ref_Sim <- Ref_Sim[ok_ref]

  # absolute anchor for ref: a_t = logit(p_ref,t)
  Abs_Obs <- Utility.Logit(Ref_Obs)
  Abs_Sim <- Utility.Logit(Ref_Sim)
  LLH_Abs <- sum(
    dnorm(Abs_Obs, mean = Abs_Sim, sd = sigma_abs, log = TRUE),
    na.rm = TRUE
  )

  # ratios for other viruses: z_{i,t} = log(p_i,t) - log(p_ref,t)
  others <- setdiff(VirNames, "Virus_1")
  LLH_Ratios <- 0
  for (v in others) {
    pObs <- DatWide[[paste0("p_obs__", v)]]
    pSim <- DatWide[[paste0("p_sim__", v)]]
    ok <- !is.na(pObs) & !is.na(pSim)
    if (!any(ok)) {
      next
    }
    z_obs <- log(pObs[ok]) - log(Ref_Obs[ok])
    z_sim <- log(pSim[ok]) - log(Ref_Sim[ok])
    LLH_Ratios <- LLH_Ratios +
      sum(
        dnorm(z_obs, mean = z_sim, sd = sigma_ratio, log = TRUE),
        na.rm = TRUE
      )
  }
  return(LLH_Abs + LLH_Ratios)
}


#' Calculate the log-likelihood based on the "Log Difference" method, which compares the log-transformed predicted probabilities to the log-transformed observed probabilities.
#' @param Dat: a data.table containing columns 'Time', 'ISOweek', 'Virus', 'p_sim' (predicted probabilities), and 'p_obs' (observed probabilities)
#' @param transform: the transformation to apply to the probabilities before calculating the difference ("logit" or "log")
#' @param sigma: standard deviation for the normal distribution of the differences
#' @param epsilon: small value to avoid numerical issues with probabilities
LLH.LogDiff <- function(
  Dat,
  transform = c("logit", "log"),
  sigma = 0.7,
  epsilon = 1e-9
) {
  transform <- match.arg(transform)

  # x_obs = transform(p_obs), x_Sim = transform(p_sim)
  # x_obs ~ Normal(x_Sim, sigma)
  if (transform == "logit") {
    x_obs <- Utility.Logit(Dat$p_obs)
    x_Sim <- Utility.Logit(Dat$p_sim)
  } else {
    x_obs <- log(Dat$p_obs)
    x_Sim <- log(Dat$p_sim)
  }

  LLH <- sum(dnorm(x_obs, mean = x_Sim, sd = sigma, log = TRUE), na.rm = TRUE)
  return(LLH)
}


#' Calculate the log PMF of a Dirichlet distribution for a given observed vector X and concentration parameters Alpha.
#' @param X: a matrix of observed proportions (T x K) where T is the number of time points and K is the number of viruses
#' @param Alpha: a matrix of concentration parameters (T x K) for the Dirichlet distribution, where Alpha[t, i] corresponds to the concentration parameter for virus i at time t
LLH.Dirichlet.PMF <- function(X, Alpha) {
  lgamma(rowSums(Alpha)) -
    rowSums(lgamma(Alpha)) +
    rowSums((Alpha - 1) * log(X))
}


#' Calculate the log-likelihood based on a Dirichlet likelihood, which compares the observed proportions of each virus to the predicted proportions from the model, accounting for overdispersion with a concentration parameter kappa.
#' @param Dat: a data.table containing columns 'Time', 'ISOweek', 'Virus', 'p_sim' (predicted probabilities), and 'p_obs' (observed probabilities)
#' @param kappa: concentration parameter for the Dirichlet distribution
#' @param c: small constant added to observed counts to avoid zero counts
LLH.Dirichlet <- function(Dat, kappa, c = 1e-6) {
  VirNames <- sort(unique(as.character(Dat$Virus)))

  DatWide <- dcast(
    Dat[, .(Time, ISOweek, Virus, p_sim, p_obs)],
    Time + ISOweek ~ Virus,
    value.var = c("p_obs", "p_sim"),
    sep = "__"
  )

  obs_cols <- paste0("p_obs__", VirNames)
  sim_cols <- paste0("p_sim__", VirNames)

  keep <- obs_cols %in% names(DatWide) & sim_cols %in% names(DatWide)
  obs_cols <- obs_cols[keep]
  sim_cols <- sim_cols[keep]

  S_obs <- as.matrix(DatWide[, obs_cols, with = FALSE])
  P_sim <- as.matrix(DatWide[, sim_cols, with = FALSE])

  S_obs <- pmax(S_obs, 1e-12)
  S_obs <- S_obs / rowSums(S_obs)

  P_total <- rowSums(P_sim)
  valid <- which(P_total > 1e-12)
  if (length(valid) == 0) {
    return(-Inf)
  }

  S_sim <- P_sim[valid, , drop = FALSE] / rowSums(P_sim[valid, , drop = FALSE])

  Alpha <- kappa * S_sim + c
  Alpha <- pmax(Alpha, 1e-10)

  X <- S_obs[valid, , drop = FALSE]
  X <- pmax(X, 1e-12)
  X <- X / rowSums(X)

  LLH <- sum(LLH.Dirichlet.PMF(X, Alpha))
  return(LLH)
}


#' Calculate the log-likelihood based on the specified method
#' @param Parm: a list of parameters for the ODE model
#' @param TargetDat: a data.table containing the observed data to compare against
#' @param Method: the method to use for calculating the likelihood
#' @return the calculated log-likelihood value
Model.RunSim.LLH <- function(
  Parm,
  after = as.Date("2015-08-31"),
  TargetDat,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  LLHArg = list()
) {
  Method <- match.arg(Method)
  Tar <- copy(TargetDat)

  # Run simulation
  SimData <- Model.RunSim.ode(Parm, Plot = FALSE, after = after)[["Data"]]
  SimData[, Time := Time - (as.integer(strftime(Time, "%u")) - 1L)] # align to Monday
  SimData <- SimData[,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(Time, ISOweek),
    .SDcols = setdiff(names(SimData), c("Time", "ISOweek"))
  ][order(Time)]

  SimData_long <- melt(
    SimData,
    id.vars = c("Time", "ISOweek"),
    variable.name = "Virus",
    value.name = "p_sim"
  )
  SimData_long[, p_sim := p_sim / Parm[["num_of_agent"]]]
  SimData_long[, p_sim := Utility.Clamp01(p_sim)]

  if (Method %in% c("RatioAbs", "LogDiff", "Dirichlet")) {
    pseudo <- LLHArg$pseudo %||% 0.5
    Tar[, p_obs := (y + pseudo) / (N + 2 * pseudo)]
    Tar[, p_obs := Utility.Clamp01(p_obs)]
  }

  MergeData <- SimData_long[
    Tar,
    on = c("ISOweek", "Virus"),
    nomatch = 0
  ]
  setorder(MergeData, ISOweek, Virus)

  # Calculate log-likelihood based on the specified method
  if (Method == "BetaBinomial") {
    rho <- LLHArg$rho %||% 0.02
    LLH <- LLH.BetaBinomial(MergeData, rho = rho)
  } else if (Method == "RatioAbs") {
    sigma_ratio <- LLHArg$sigma_ratio %||% 0.7
    sigma_abs <- LLHArg$sigma_abs %||% 0.7

    LLH <- LLH.RatioAbs(
      MergeData,
      VirNames = paste0("Virus_", 1:4),
      sigma_ratio = sigma_ratio,
      sigma_abs = sigma_abs
    )
  } else if (Method == "LogDiff") {
    sigma <- LLHArg$sigma %||% 0.7
    transform <- LLHArg$transform %||% "logit"

    LLH <- LLH.LogDiff(
      MergeData,
      transform = transform,
      sigma = sigma
    )
  } else if (Method == "Dirichlet") {
    kappa <- LLHArg$kappa %||% 100
    c_add <- LLHArg$c %||% 1e-6

    LLH <- LLH.Dirichlet(
      Dat = MergeData,
      kappa = kappa,
      c = c_add
    )
  } else {
    stop("Invalid Method")
  }
  return(LLH)
}

# Markov Chain Monte Carlo sampling ----------------------
#' Log-prior for the competition parameters, assuming they follow a log-normal distribution with mean 0 and specified standard deviation on the log scale.
#' @param comp: a vector of competition parameters for the viruses
#' @param sdlog: standard deviation of the log-normal distribution (default is 1)
MCMC.LogPrior.Comp <- function(comp, sdlog = 1) {
  z <- log(comp)
  z <- z - mean(z) # identifiability
  return(sum(dnorm(z, mean = 0, sd = sdlog, log = TRUE)))
}


#' Generate a new proposal for the competition parameters by adding a random perturbation to the log-transformed parameters, and then normalizing to ensure identifiability.
#' @param Parm: a vector of current parameter values (only the first 4 elements are used for the competition parameters)
#' @param step1: the maximum absolute value of the random perturbation added to the log-transformed parameters (default is 1)
MCMC.Proposal <- function(Parm, step1 = 0.1) {
  ParmReal <- log(Parm)
  ParmUpdate <- ParmReal + runif(length(Parm), -step1, step1)

  # identifiability: geometric mean(comp) = 1
  ParmUpdate <- ParmUpdate - mean(ParmUpdate)

  NewParm <- exp(ParmUpdate)
  return(NewParm)
}


#' Run MCMC sampling using the Metropolis-Hastings algorithm to estimate the competition parameters of the ODE model.
#' @param Initial: a vector of initial values for the competition parameters
#' @param n_iterations: the number of MCMC iterations to run
#' @param step: the step size for the proposal distribution
#' @param TargetDat: a data.table containing the observed data to compare against
#' @param Method: the method to use for calculating the likelihood
#' @param BaseParm: a list of base parameters for the ODE model, which will be updated with the proposed competition parameters in each iteration
#' @param after: only include data after this date when calculating the likelihood
#' @param LLHArg: a list of additional arguments to pass to the likelihood function
#' @param LogPriorFun: a function that takes a vector of competition parameters and returns the log-prior probability of those parameters
#' @return a matrix containing the MCMC samples of the competition parameters, with additional attributes for the log-likelihood, log-prior, log-posterior, acceptance indicators, and acceptance rate
MCMC.MH <- function(
  Initial,
  n_iterations,
  step = 0.1,
  TargetDat = TargetDat,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  BaseParm = Parameter.Create(),
  after = as.Date("2015-08-31"),
  LLHArg = list(),
  LogPriorFun = NULL,
  verbose = TRUE
) {
  Method <- match.arg(Method)

  if (is.null(LogPriorFun)) {
    LogPriorFun <- function(comp) MCMC.LogPrior.Comp(comp, sdlog = 1)
  }

  Initial <- as.numeric(Initial)
  # normalize initial comp
  Initial <- Initial / exp(mean(log(Initial)))
  n_comp <- length(Initial)
  chain <- matrix(NA_real_, nrow = n_iterations, ncol = n_comp)
  colnames(chain) <- paste0("comp_", seq_len(n_comp))

  llh_trace <- rep(NA_real_, n_iterations)
  lpr_trace <- rep(NA_real_, n_iterations)
  lpo_trace <- rep(NA_real_, n_iterations)
  accepted <- rep(FALSE, n_iterations)

  chain[1, ] <- Initial

  Parm_cur <- BaseParm
  Parm_cur[["comp"]] <- chain[1, ]

  llh_trace[1] <- Model.RunSim.LLH(
    Parm = Parm_cur,
    after = after,
    TargetDat = TargetDat,
    Method = Method,
    LLHArg = LLHArg
  )
  lpr_trace[1] <- LogPriorFun(chain[1, ])
  lpo_trace[1] <- llh_trace[1] + lpr_trace[1]

  pb <- progress_bar$new(
    total = n_iterations,
    clear = TRUE,
    format = "  [:bar] :percent :etas"
  )
  pb$tick()

  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal(Parm = chain[i - 1, ], step = step)
    Parm_prop <- BaseParm
    Parm_prop[["comp"]] <- proposal

    llh_prop <- Model.RunSim.LLH(
      Parm = Parm_prop,
      after = after,
      TargetDat = TargetDat,
      Method = Method,
      LLHArg = LLHArg
    )
    lpr_prop <- LogPriorFun(proposal)
    lpo_prop <- llh_prop + lpr_prop
    log_alpha <- lpo_prop - lpo_trace[i - 1]

    log_alpha <- lpo_prop - lpo_trace[i - 1]

    if (log(runif(1)) < log_alpha) {
      chain[i, ] <- proposal
      llh_trace[i] <- llh_prop
      lpr_trace[i] <- lpr_prop
      lpo_trace[i] <- lpo_prop
      accepted[i] <- TRUE
    } else {
      chain[i, ] <- chain[i - 1, ]
      llh_trace[i] <- llh_trace[i - 1]
      lpr_trace[i] <- lpr_trace[i - 1]
      lpo_trace[i] <- lpo_trace[i - 1]
    }
    pb$tick()
  }

  attr(chain, "LLH") <- llh_trace
  attr(chain, "LogPrior") <- lpr_trace
  attr(chain, "LogPost") <- lpo_trace
  attr(chain, "Accepted") <- accepted
  attr(chain, "AcceptanceRate") <- mean(accepted[-1])
  attr(chain, "Method") <- Method
  attr(chain, "LLHArg") <- LLHArg

  return(chain)
}
