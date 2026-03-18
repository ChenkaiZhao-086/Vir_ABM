# ODE model functions

#' Create the initial state vector for the ODE model.
#' @param population: total population size
#' @param initial_seeds: the number of infected cases in the beginning
#' @param Base_Immu: the proportion of population that is initially immune
#' @param n_virus: the number of viruses in the model
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
#' @param Sim: the raw output from the ODE simulation (a matrix)
#' @param virus_names: optional vector of virus names to rename the incidence columns
#' @param after: only include data after this date (default "2019-01-01")
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
#' @param Parm: a list of parameters for the ODE model
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

# Likelihood / Inference / MCMC
Utility.Clamp01 <- function(p, epsilon = 1e-9) {
  pmin(pmax(p, epsilon), 1 - epsilon)
}

Utility.Logit <- function(p) log(p / (1 - p))

`%||%` <- function(x, y) if (is.null(x)) y else x


# Fast preparation of parameters for the ODE simulation -----------------------------------

#' Take the base parameters and prepare the initial state and week mapping for inference.
#' @param TargetDat: a data.table containing the target data
#' @param viruses: a character vector of virus names to include
#' @return a list containing the prepared target data
Target.Prepare <- function(TargetDat, viruses = NULL) {
  DT <- copy(TargetDat)

  if (is.null(viruses)) {
    viruses <- sort(unique(as.character(DT$Virus)))
  }
  DT <- DT[Virus %in% viruses]

  if (nrow(DT) == 0L) {
    stop("TargetDat has no rows after virus filtering.")
  }

  # 将周数据预处理成矩阵，避免 MCMC 中反复 melt/dcast
  DT_y <- DT[, .(y = sum(y, na.rm = TRUE)), by = .(ISOweek, Virus)]
  DT_N <- DT[,
    .(
      N = if (all(is.na(N))) NA_real_ else max(N, na.rm = TRUE)
    ),
    by = .(ISOweek, Virus)
  ]

  y_wide <- dcast(
    DT_y,
    ISOweek ~ Virus,
    value.var = "y",
    fill = NA_real_
  )
  N_wide <- dcast(
    DT_N,
    ISOweek ~ Virus,
    value.var = "N",
    fill = NA_real_
  )

  for (v in viruses) {
    if (!v %in% names(y_wide)) {
      y_wide[[v]] <- NA_real_
    }
    if (!v %in% names(N_wide)) N_wide[[v]] <- NA_real_
  }

  setcolorder(y_wide, c("ISOweek", viruses))
  setcolorder(N_wide, c("ISOweek", viruses))

  ord <- order(y_wide$ISOweek)
  y_wide <- y_wide[ord]
  N_wide <- N_wide[match(y_wide$ISOweek, N_wide$ISOweek)]

  y_mat <- as.matrix(y_wide[, ..viruses])
  N_mat <- as.matrix(N_wide[, ..viruses])

  sameN <- apply(N_mat, 1, function(z) {
    z <- z[is.finite(z)]
    length(z) <= 1L || all(abs(z - z[1]) < 1e-8)
  })
  if (!all(sameN)) {
    warning(
      "Some weeks have non-identical N across viruses; Dirichlet weight will use the first finite N in that week."
    )
  }

  N_week <- apply(N_mat, 1, function(z) {
    z <- z[is.finite(z)]
    if (length(z) == 0L) NA_real_ else z[1]
  })

  return(list(
    ISOweek = y_wide$ISOweek,
    viruses = viruses,
    y = y_mat,
    N = N_mat,
    N_week = N_week,
    n_virus = length(viruses)
  ))
}


#' Create a mapping from dates to ISO weeks for the simulation period.
#' This is used to aggregate daily simulation output into weekly data for comparison with targets.
#' @param Parm: the list of parameters containing the simulation start and end dates
#' @param after: only include dates after this date (default "2015-08-31")
#' @return a list containing the times, keep vector, group factor, ISOweek labels, and number of days in each week
Model.MakeWeekMap <- function(Parm, after = as.Date("2015-08-31")) {
  times <- as.numeric(seq(
    from = as.Date(Parm[["year_start"]]),
    to = as.Date(Parm[["year_end"]]),
    by = 1
  ))
  dates <- as.Date(times, origin = "1970-01-01")
  monday <- dates - (as.integer(strftime(dates, "%u")) - 1L)

  keep <- dates >= after
  monday_keep <- monday[keep]
  week_levels <- unique(monday_keep)
  group <- factor(monday_keep, levels = week_levels, ordered = TRUE)

  return(list(
    times = times,
    keep = keep,
    group = group,
    ISOweek = strftime(week_levels, "%G-W%V"),
    n_day_in_week = as.integer(table(group))
  ))
}


#' Take the base parameters and prepare the initial state and week mapping for inference.
#' @param TargetDat: a data.table containing the target data
#' @param BaseParm: a list of base parameters for the ODE model
#' @param after: only include data after this date (default "2015-08-31")
#' @param viruses: a character vector of virus names to include (default NULL, meaning all viruses in TargetDat)
#' @return a list containing the prepared target data, week mapping, simulation index, initial state, and after date
Inference.Setup <- function(
  TargetDat,
  BaseParm,
  after = as.Date("2015-08-31"),
  viruses = NULL
) {
  if (is.null(viruses)) {
    viruses <- paste0("Virus_", seq_along(BaseParm[["beta0"]]))
  }

  target <- Target.Prepare(TargetDat, viruses = viruses)
  week_map <- Model.MakeWeekMap(BaseParm, after = after)

  sim_index <- match(target$ISOweek, week_map$ISOweek)
  keep <- !is.na(sim_index)

  if (!all(keep)) {
    warning(
      sum(!keep),
      " target weeks are outside the simulation window and were dropped."
    )
  }
  if (!any(keep)) {
    stop(
      "No overlapping ISO weeks between TargetDat and the simulation window."
    )
  }

  target$ISOweek <- target$ISOweek[keep]
  target$y <- target$y[keep, , drop = FALSE]
  target$N <- target$N[keep, , drop = FALSE]
  target$N_week <- target$N_week[keep]

  init_state <- Get.InitState(
    population = BaseParm[["num_of_agent"]],
    initial_seeds = BaseParm[["initial_seeds"]],
    Base_Immu = BaseParm[["base_immune"]],
    n_virus = length(BaseParm[["beta0"]])
  )

  return(list(
    target = target,
    week_map = week_map,
    sim_index = sim_index[keep],
    init_state = init_state,
    after = after
  ))
}


#' Run the ODE simulation and directly aggregate the daily output into weekly incidence for likelihood calculation.
#' @param Parm: a list of parameters for the ODE model
#' @param WeekMap a list containing the times, keep vector, group factor, ISOweek labels, and number of days in each week (output of Model.MakeWeekMap)
#' @param metric the metric to use for likelihood calculation ("Inc_sum" for total weekly incidence, "I_mean" for average weekly prevalence)
#' @param InitState the initial state vector for the ODE simulation (optional)
#' @return a list containing the raw simulation output and the processed weekly incidence data
Model.RunSim.Weekly <- function(
  Parm,
  WeekMap,
  metric = c("Inc_sum", "I_mean"),
  InitState = NULL
) {
  metric <- match.arg(metric)
  n_virus <- length(Parm[["beta0"]])

  if (is.null(InitState) || length(InitState) != 3 * n_virus) {
    InitState <- Get.InitState(
      population = Parm[["num_of_agent"]],
      initial_seeds = Parm[["initial_seeds"]],
      Base_Immu = Parm[["base_immune"]],
      n_virus = n_virus
    )
  }

  sim <- ode(
    y = InitState,
    times = WeekMap$times,
    func = ParmInferenceCpp,
    parms = Parm,
    method = "rk4"
  )
  sim <- as.matrix(sim)

  if (metric == "Inc_sum") {
    cols <- grep("^Inc_\\d+$", colnames(sim))
    if (length(cols) == 0L) {
      stop("No incidence columns matching '^Inc_\\\\d+$' found in ODE output.")
    }
    daily <- sim[WeekMap$keep, cols, drop = FALSE]
    weekly <- rowsum(daily, group = WeekMap$group, reorder = FALSE)
  } else {
    cols <- grep("^I_\\d+$", colnames(sim))
    if (length(cols) == 0L) {
      stop("No state columns matching '^I_\\\\d+$' found in ODE output.")
    }
    daily <- sim[WeekMap$keep, cols, drop = FALSE]
    weekly <- rowsum(daily, group = WeekMap$group, reorder = FALSE)
    weekly <- sweep(weekly, 1, WeekMap$n_day_in_week, "/")
  }

  colnames(weekly) <- paste0("Virus_", seq_len(ncol(weekly)))

  return(list(
    ISOweek = WeekMap$ISOweek,
    mat = weekly
  ))
}


# Likelihoods ---------------------------------------------------------------

#' Calculate the log-likelihood using a Beta-Binomial model.
#' @param y vector of observed counts
#' @param n vector of total counts
#' @param p vector of probabilities
#' @param rho overdispersion parameter
#' @param epsilon small value to avoid numerical issues
LLH.BetaBinomial.PMF <- function(y, n, p, rho, epsilon = 1e-9) {
  p <- Utility.Clamp01(p, epsilon)
  rho <- pmin(pmax(rho, epsilon), 1 - epsilon)
  theta <- (1 - rho) / rho
  a <- p * theta
  b <- (1 - p) * theta
  out <- lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)
  bad <- !is.finite(y) |
    !is.finite(n) |
    !is.finite(p) |
    (y < 0) |
    (n < 0) |
    (y > n)
  out[bad] <- NA_real_
  return(out)
}


#' Calculate the log-likelihood for the Beta-Binomial model.
#' This function sums the log-likelihood contributions from all valid entries in the matrices.
#' @param y_mat matrix of observed counts
#' @param N_mat matrix of total counts
#' @param p_sim_mat matrix of simulated probabilities
#' @param rho overdispersion parameter
#' @param epsilon small value to avoid numerical issues
LLH.BetaBinomial <- function(
  y_mat,
  N_mat,
  p_sim_mat,
  rho = 0.02,
  epsilon = 1e-9
) {
  valid <- is.finite(y_mat) &
    is.finite(N_mat) &
    is.finite(p_sim_mat) &
    (y_mat >= 0) &
    (N_mat >= 0) &
    (y_mat <= N_mat)

  if (!any(valid)) {
    return(-Inf)
  }

  return(sum(
    LLH.BetaBinomial.PMF(
      y = y_mat[valid],
      n = N_mat[valid],
      p = p_sim_mat[valid],
      rho = rho,
      epsilon = epsilon
    ),
    na.rm = TRUE
  ))
}


#' Calculate the log-likelihood for the Ratio-Absolute method, which combines the likelihood of the absolute incidence and the ratios between viruses.
#' This method uses the first virus as a reference and models the log ratios of other viruses to the reference, as well as the absolute logit-transformed incidence of the reference virus.
#' @param p_obs_mat matrix of observed probabilities (incidence proportions)
#' @param p_sim_mat matrix of simulated probabilities
#' @param VirNames vector of virus names corresponding to the columns of the matrices
#' @param sigma_ratio standard deviation for the ratio component of the likelihood
#' @param sigma_abs standard deviation for the absolute component of the likelihood
#' @param epsilon small value to avoid numerical issues
LLH.RatioAbs <- function(
  p_obs_mat,
  p_sim_mat,
  VirNames = NULL,
  sigma_ratio = 0.7,
  sigma_abs = 0.7,
  epsilon = 1e-9
) {
  p_obs_mat <- Utility.Clamp01(p_obs_mat, epsilon)
  p_sim_mat <- Utility.Clamp01(p_sim_mat, epsilon)

  if (is.null(VirNames)) {
    VirNames <- paste0("Virus_", seq_len(ncol(p_obs_mat)))
  }

  ref_idx <- match("Virus_1", VirNames)
  if (is.na(ref_idx)) {
    return(-Inf)
  }

  Ref_Obs <- p_obs_mat[, ref_idx]
  Ref_Sim <- p_sim_mat[, ref_idx]

  ok_ref <- is.finite(Ref_Obs) & is.finite(Ref_Sim)
  if (!any(ok_ref)) {
    return(-Inf)
  }

  p_obs_mat <- p_obs_mat[ok_ref, , drop = FALSE]
  p_sim_mat <- p_sim_mat[ok_ref, , drop = FALSE]
  Ref_Obs <- Ref_Obs[ok_ref]
  Ref_Sim <- Ref_Sim[ok_ref]

  Abs_Obs <- Utility.Logit(Ref_Obs)
  Abs_Sim <- Utility.Logit(Ref_Sim)
  LLH_Abs <- sum(
    dnorm(Abs_Obs, mean = Abs_Sim, sd = sigma_abs, log = TRUE),
    na.rm = TRUE
  )

  others <- setdiff(seq_along(VirNames), ref_idx)
  LLH_Ratios <- 0

  for (j in others) {
    pObs <- p_obs_mat[, j]
    pSim <- p_sim_mat[, j]
    ok <- is.finite(pObs) & is.finite(pSim)
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


#' Calculate the log-likelihood for the Log-Difference method, which models the log-transformed probabilities of each virus independently.
#' This method applies a transformation (logit or log) to the observed and simulated probabilities and then models the transformed values with a normal distribution.
#' @param p_obs_mat matrix of observed probabilities (incidence proportions)
#' @param p_sim_mat matrix of simulated probabilities
#' @param transform the transformation to apply to the probabilities ("logit" or "log")
#' @param sigma standard deviation for the likelihood
#' @param epsilon small value to avoid numerical issues
LLH.LogDiff <- function(
  p_obs_mat,
  p_sim_mat,
  transform = c("logit", "log"),
  sigma = 0.7,
  epsilon = 1e-9
) {
  transform <- match.arg(transform)

  p_obs_mat <- Utility.Clamp01(p_obs_mat, epsilon)
  p_sim_mat <- Utility.Clamp01(p_sim_mat, epsilon)

  valid <- is.finite(p_obs_mat) & is.finite(p_sim_mat)
  if (!any(valid)) {
    return(-Inf)
  }

  if (transform == "logit") {
    x_obs <- Utility.Logit(p_obs_mat[valid])
    x_sim <- Utility.Logit(p_sim_mat[valid])
  } else {
    x_obs <- log(p_obs_mat[valid])
    x_sim <- log(p_sim_mat[valid])
  }

  return(sum(dnorm(x_obs, mean = x_sim, sd = sigma, log = TRUE), na.rm = TRUE))
}


#' Calculate the log-likelihood for the Dirichlet method.
#' @param X matrix of observed proportions
#' @param Alpha matrix of concentration parameters
LLH.Dirichlet.PMF <- function(X, Alpha) {
  lgamma(rowSums(Alpha)) -
    rowSums(lgamma(Alpha)) +
    rowSums((Alpha - 1) * log(X))
}


#' Calculate the log-likelihood for the Dirichlet method.
#' This method first converts the observed counts into proportions, then uses the simulated proportions to construct the concentration parameters of the Dirichlet distribution, and finally calculates the log-likelihood using the Dirichlet PMF.
#' @param y_mat matrix of observed counts
#' @param sim_mat matrix of simulated counts or proportions
#' @param kappa concentration parameter for the Dirichlet distribution
#' @param N_vec optional vector of total counts for weighting the concentration parameter (used if weight is "N" or "sum_y_capN")
#' @param c small value to avoid numerical issues when converting counts to proportions
#' @param weight method for weighting the concentration parameter ("sum_y_capN" uses the sum of observed counts capped by N_vec, "N"
LLH.Dirichlet <- function(
  y_mat,
  sim_mat,
  kappa,
  N_vec = NULL,
  c = 1e-8,
  weight = c("sum_y_capN", "N", "none")
) {
  # 修改后的 Dirichlet：
  # 1) 对 share 建模：y_i / sum_j y_j
  # 2) sim 端用各病毒 share
  # 3) 跳过总和为 0 的周
  # 4) concentration 用 kappa * weight_t
  weight <- match.arg(weight)

  if (!is.finite(kappa) || kappa <= 0) {
    return(-Inf)
  }

  complete_obs <- apply(is.finite(y_mat), 1, all)
  complete_sim <- apply(is.finite(sim_mat), 1, all)

  obs_total <- rowSums(y_mat, na.rm = TRUE)
  sim_total <- rowSums(sim_mat, na.rm = TRUE)

  keep <- complete_obs & complete_sim & (obs_total > 0) & (sim_total > 0)
  if (!any(keep)) {
    return(-Inf)
  }

  Y <- y_mat[keep, , drop = FALSE]
  M <- sim_mat[keep, , drop = FALSE]

  S_obs <- Y + c
  S_obs <- sweep(S_obs, 1, rowSums(S_obs), "/")

  S_sim <- pmax(M, 0) + c
  S_sim <- sweep(S_sim, 1, rowSums(S_sim), "/")

  if (weight == "none") {
    w <- rep(1, nrow(Y))
  } else if (weight == "N") {
    if (is.null(N_vec)) {
      stop("weight = 'N' requires N_vec.")
    }
    N_use <- N_vec[keep]
    if (any(!is.finite(N_use))) {
      stop("weight = 'N' but some N_vec are missing.")
    }
    w <- N_use
  } else {
    # sum_y_capN
    w <- rowSums(Y)
    if (!is.null(N_vec)) {
      N_use <- N_vec[keep]
      ok <- is.finite(N_use)
      w[ok] <- pmin(w[ok], N_use[ok])
    }
  }

  w <- pmax(w, 1)

  Alpha <- sweep(S_sim, 1, kappa * w, "*")
  Alpha <- pmax(Alpha, 1e-12)

  return(sum(LLH.Dirichlet.PMF(S_obs, Alpha)))
}


# Main log-likelihood wrapper -----------------------------------------------

#' Run the ODE simulation and calculate the log-likelihood based on the specified method.
#' @param Parm list of parameters for the ODE simulation
#' @param after date after which to consider the data
#' @param TargetDat target data for the simulation
#' @param Prep preprocessed data for the simulation
#' @param Method inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param LLHArg list of additional arguments for the log-likelihood calculation, including:
#'   - viruses: vector of virus names to include in the likelihood calculation
#'   - proxy: metric to use for the likelihood calculation ("Inc_sum" or "I_mean")
#'   - epsilon: small value to avoid numerical issues
#'   - rho: overdispersion parameter for the Beta-Binomial method
#'   - pseudo: pseudo-count for calculating observed probabilities in RatioAbs and LogDiff methods
#'   - sigma_ratio: standard deviation for the ratio component in RatioAbs method
#'   - sigma_abs: standard deviation for the absolute component in RatioAbs method
#'   - sigma: standard deviation for the LogDiff method
#'   - transform: transformation to apply in LogDiff method ("logit" or "log")
#'   - kappa: concentration parameter for the Dirichlet method
#'   - c: small value to avoid numerical issues in the Dirichlet method
#'   - weight: weighting method for the Dirichlet method ("sum_y_capN", "N", "none")
Model.RunSim.LLH <- function(
  Parm,
  after = as.Date("2015-08-31"),
  TargetDat = NULL,
  Prep = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  LLHArg = list()
) {
  Method <- match.arg(Method)

  if (is.null(Prep)) {
    if (is.null(TargetDat)) {
      stop("Provide either Prep or TargetDat.")
    }
    Prep <- Inference.Setup(
      TargetDat = TargetDat,
      BaseParm = Parm,
      after = after,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_along(Parm[["beta0"]]))
    )
  }

  proxy <- LLHArg$proxy %||% "Inc_sum"
  epsilon <- LLHArg$epsilon %||% 1e-9

  sim <- Model.RunSim.Weekly(
    Parm = Parm,
    WeekMap = Prep$week_map,
    metric = proxy,
    InitState = Prep$init_state
  )

  target <- Prep$target
  sim_mat <- sim$mat[Prep$sim_index, target$viruses, drop = FALSE]
  p_sim_mat <- Utility.Clamp01(
    sim_mat / Parm[["num_of_agent"]],
    epsilon = epsilon
  )

  if (Method == "BetaBinomial") {
    rho <- LLHArg$rho %||% 0.02
    return(LLH.BetaBinomial(
      y_mat = target$y,
      N_mat = target$N,
      p_sim_mat = p_sim_mat,
      rho = rho,
      epsilon = epsilon
    ))
  }

  if (Method %in% c("RatioAbs", "LogDiff")) {
    pseudo <- LLHArg$pseudo %||% 0.5
    p_obs_mat <- (target$y + pseudo) / (target$N + 2 * pseudo)
    p_obs_mat[!is.finite(target$N) | target$N < 0] <- NA_real_
    p_obs_mat <- Utility.Clamp01(p_obs_mat, epsilon = epsilon)
  }

  if (Method == "RatioAbs") {
    sigma_ratio <- LLHArg$sigma_ratio %||% 0.7
    sigma_abs <- LLHArg$sigma_abs %||% 0.7

    return(LLH.RatioAbs(
      p_obs_mat = p_obs_mat,
      p_sim_mat = p_sim_mat,
      VirNames = target$viruses,
      sigma_ratio = sigma_ratio,
      sigma_abs = sigma_abs,
      epsilon = epsilon
    ))
  }

  if (Method == "LogDiff") {
    sigma <- LLHArg$sigma %||% 0.7
    transform <- LLHArg$transform %||% "logit"

    return(LLH.LogDiff(
      p_obs_mat = p_obs_mat,
      p_sim_mat = p_sim_mat,
      transform = transform,
      sigma = sigma,
      epsilon = epsilon
    ))
  }

  # Dirichlet
  kappa <- LLHArg$kappa %||% 100
  c_add <- LLHArg$c %||% 1e-8
  weight <- LLHArg$weight %||% "sum_y_capN"

  LLH.Dirichlet(
    y_mat = target$y,
    sim_mat = sim_mat,
    kappa = kappa,
    N_vec = target$N_week,
    c = c_add,
    weight = weight
  )
}


# Parameterization of comp --------------------------------------------------

# 采用 centered log-comp:
# eta_full = (eta_1, ..., eta_{K-1}, -sum eta)
# comp = exp(eta_full)
# 好处：
# 1) comp 始终 > 0
# 2) 几何均值固定为 1
# 3) 保留原始 comp 的正参数解释

#' Convert the free eta parameters to the full eta vector by appending the negative sum of the free parameters.
#' This ensures that the geometric mean of the resulting comp values will be 1, which helps with identifiability and interpretability.
#' @param eta_free vector of free eta parameters (length K-1)
#' @param n_virus total number of viruses (K)
Comp.EtaFull <- function(eta_free, n_virus = length(eta_free) + 1L) {
  stopifnot(length(eta_free) == n_virus - 1L)
  c(eta_free, -sum(eta_free))
}

Comp.FromEta <- function(eta_free, n_virus = length(eta_free) + 1L) {
  exp(Comp.EtaFull(eta_free, n_virus = n_virus))
}


# 兼容旧思路：comp 的 log-normal 型先验
MCMC.LogPrior.Comp <- function(comp = NULL, eta = NULL, sd = 1) {
  if (!is.null(comp)) {
    if (any(!is.finite(comp)) || any(comp <= 0)) {
      return(-Inf)
    }
    eta <- log(comp)
  }
  if (is.null(eta)) {
    stop("Provide either comp or eta.")
  }
  if (any(!is.finite(eta))) {
    return(-Inf)
  }

  eta <- eta - mean(eta)
  sum(dnorm(eta, mean = 0, sd = sd, log = TRUE))
}


# Method-specific theta spec ------------------------------------------------

#' Create a specification for the parameters of the inference model based on the method and number of viruses.
#' @param Method: the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus: the number of viruses
#' @return a list containing the method, number of viruses, component names, auxiliary names, parameter names, and number of parameters
Inference.ParamSpec <- function(
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L
) {
  Method <- match.arg(Method)

  comp_names <- paste0("eta_", seq_len(n_virus - 1L))
  aux_names <- switch(
    Method,
    BetaBinomial = "z_rho",
    Dirichlet = "log_kappa",
    RatioAbs = c("log_sigma_ratio", "log_sigma_abs"),
    LogDiff = "log_sigma"
  )

  return(list(
    Method = Method,
    n_virus = n_virus,
    comp_names = comp_names,
    aux_names = aux_names,
    par_names = c(comp_names, aux_names),
    n_par = length(c(comp_names, aux_names))
  ))
}


#' Decode the parameter vector theta into the components and auxiliary parameters.
#' This function takes the free eta parameters, reconstructs the full eta vector, and then calculates the comp values. It also extracts the auxiliary parameters based on the specified method.
#' @param theta vector of parameters to decode
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
#' @return a list containing the full eta vector, comp values, and auxiliary parameters
Inference.DecodeTheta <- function(
  theta,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)

  if (length(theta) != spec$n_par) {
    stop("theta length does not match method specification.")
  }

  eta_free <- theta[seq_len(n_virus - 1L)]
  eta_full <- Comp.EtaFull(eta_free, n_virus = n_virus)
  comp <- exp(eta_full)

  out <- list(
    eta = eta_full,
    comp = comp
  )

  idx <- n_virus
  if (Method == "BetaBinomial") {
    out$rho <- plogis(theta[idx])
  } else if (Method == "Dirichlet") {
    out$kappa <- exp(theta[idx])
  } else if (Method == "RatioAbs") {
    out$sigma_ratio <- exp(theta[idx])
    out$sigma_abs <- exp(theta[idx + 1L])
  } else if (Method == "LogDiff") {
    out$sigma <- exp(theta[idx])
  }

  return(out)
}

# Initial values ------------------------------------------------------------

#' Create the initial parameter vector for MCMC based on the specified method and initial values.
#' @param Initial: a list of initial values for the parameters, or a numeric vector of the same length as the number of parameters
#' @param Method: the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus: the number of viruses
#' @return a numeric vector of initial parameter values for MCMC
MCMC.MakeInitial <- function(
  Initial = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L
) {
  Method <- match.arg(Method)

  aux_default <- switch(
    Method,
    BetaBinomial = qlogis(0.02),
    Dirichlet = log(100),
    RatioAbs = c(log(0.7), log(0.7)),
    LogDiff = log(0.7)
  )

  theta0 <- c(rep(0, n_virus - 1L), aux_default)

  if (is.null(Initial)) {
    return(theta0)
  }

  if (is.list(Initial)) {
    if (!is.null(Initial$eta)) {
      eta0 <- as.numeric(Initial$eta)
      if (length(eta0) == n_virus) {
        # Center the eta values to ensure identifiability
        eta0 <- eta0 - mean(eta0)
        theta0[seq_len(n_virus - 1L)] <- eta0[seq_len(n_virus - 1L)]
      } else if (length(eta0) == n_virus - 1L) {
        theta0[seq_len(n_virus - 1L)] <- eta0
      } else {
        stop("Initial$eta has wrong length.")
      }
    }

    if (!is.null(Initial$comp)) {
      comp0 <- as.numeric(Initial$comp)
      if (
        length(comp0) != n_virus || any(comp0 <= 0) || any(!is.finite(comp0))
      ) {
        stop("Initial$comp must be a positive vector of length n_virus.")
      }
      # Center the log-comp values to ensure identifiability
      eta0 <- log(comp0)
      eta0 <- eta0 - mean(eta0)
      theta0[seq_len(n_virus - 1L)] <- eta0[seq_len(n_virus - 1L)]
    }

    if (Method == "BetaBinomial" && !is.null(Initial$rho)) {
      theta0[n_virus] <- qlogis(Initial$rho)
    }
    if (Method == "Dirichlet" && !is.null(Initial$kappa)) {
      theta0[n_virus] <- log(Initial$kappa)
    }
    if (Method == "RatioAbs") {
      if (!is.null(Initial$sigma_ratio)) {
        theta0[n_virus] <- log(Initial$sigma_ratio)
      }
      if (!is.null(Initial$sigma_abs)) {
        theta0[n_virus + 1L] <- log(Initial$sigma_abs)
      }
    }
    if (Method == "LogDiff" && !is.null(Initial$sigma)) {
      theta0[n_virus] <- log(Initial$sigma)
    }

    return(theta0)
  }

  Initial <- as.numeric(Initial)
  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)

  if (length(Initial) == spec$n_par) {
    return(Initial)
  }

  if (length(Initial) == n_virus) {
    if (all(is.finite(Initial)) && all(Initial > 0)) {
      eta0 <- log(Initial)
    } else {
      eta0 <- Initial
    }
    eta0 <- eta0 - mean(eta0)
    theta0[seq_len(n_virus - 1L)] <- eta0[seq_len(n_virus - 1L)]
    return(theta0)
  }

  if (length(Initial) == n_virus - 1L) {
    theta0[seq_len(n_virus - 1L)] <- Initial
    return(theta0)
  }

  stop("Initial has unsupported length.")
}


# Proposal ------------------------------------------------------------------

#' Resolve the step sizes for the proposal distribution in MCMC.
#' @param step a scalar or a list specifying the standard deviation of the proposal distribution for the component parameters and auxiliary parameters. If a scalar is provided, it will be used for all parameters
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
MCMC.ResolveStep <- function(
  step = 0.1,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)

  if (is.numeric(step) && length(step) == 1L) {
    return(list(
      comp = rep(step, n_virus - 1L),
      aux = rep(step, length(spec$aux_names))
    ))
  }

  if (!is.list(step)) {
    stop("step must be either a scalar or a list.")
  }

  step_comp <- rep(step$comp %||% 0.10, n_virus - 1L)
  step_aux <- rep(step$aux %||% 0.10, length(spec$aux_names))

  if (length(spec$aux_names) > 0L) {
    for (j in seq_along(spec$aux_names)) {
      nm <- spec$aux_names[j]
      if (!is.null(step[[nm]])) {
        step_aux[j] <- step[[nm]]
      }
    }
  }

  return(list(comp = step_comp, aux = step_aux))
}


#' Generate a proposal parameter vector for MCMC by adding normally distributed random noise to the current parameter vector Theta.
#' The noise is added separately to the component parameters and the auxiliary parameters based on the specified method and step sizes.
#' @param Theta a numeric vector of the current parameter values (length should match the number of parameters for the specified method)
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
#' @param step a scalar or a list specifying the standard deviation of the proposal distribution for the component parameters and auxiliary parameters. If a scalar is provided, it will be used for all parameters
MCMC.Proposal <- function(
  Theta,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  step = 0.1
) {
  Method <- match.arg(Method)
  steps <- MCMC.ResolveStep(step = step, Method = Method, n_virus = n_virus)

  out <- Theta
  idx_comp <- seq_len(n_virus - 1L)
  out[idx_comp] <- out[idx_comp] + rnorm(length(idx_comp), 0, steps$comp)

  if (length(steps$aux) > 0L) {
    idx_aux <- (n_virus):(length(Theta))
    out[idx_aux] <- out[idx_aux] + rnorm(length(idx_aux), 0, steps$aux)
  }

  return(out)
}


# Adaptative proposal -------------------------------------------------------------

#' Resolve the step sizes for the adaptive proposal distribution in MCMC and return them as a single vector.
#' This function calls MCMC.ResolveStep to get the step sizes for the component parameters and auxiliary parameters, and then concatenates them into a single vector.
#' @param step a scalar or a list specifying the standard deviation of the proposal distribution for the component parameters and auxiliary parameters. If a scalar is provided, it will be used for all parameters
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
MCMC.ResolveStepVec <- function(
  step = 0.1,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L
) {
  Method <- match.arg(Method)
  st <- MCMC.ResolveStep(step = step, Method = Method, n_virus = n_virus)
  return(c(st$comp, st$aux))
}


#' Initialize the covariance estimation for the adaptive MCMC proposal distribution.
#' This function creates a list to store the number of samples (n), the current mean vector, and the sum of squared deviations matrix (S) for the parameters.
#' This will be updated iteratively as new samples are accepted during the MCMC process.
#' @param NParams the number of parameters in the model (length of the parameter vector theta)
MCMC.CovInit <- function(NParams) {
  list(
    n = 0L,
    mean = numeric(NParams),
    S = matrix(0, nrow = NParams, ncol = NParams)
  )
}


#' Update the covariance matrix for the adaptive MCMC proposal distribution using Welford's online algorithm.
#' This function updates the mean and the sum of squared deviations (S) based on the new parameter vector x. The covariance matrix can be calculated from S and n when needed.
#' @param obj a list containing the current state of the covariance estimation, including n (number of samples), mean (current mean vector), and S (current sum of squared deviations matrix)
#' @param x the new parameter vector to incorporate into the covariance estimation
MCMC.CovUpdate <- function(obj, x) {
  x <- as.numeric(x)
  obj$n <- obj$n + 1L

  if (obj$n == 1L) {
    obj$mean <- x
    return(obj)
  }

  delta <- x - obj$mean # difference between the new sample and the current mean
  obj$mean <- obj$mean + delta / obj$n # update the mean with the new sample
  obj$S <- obj$S + tcrossprod(delta, x - obj$mean) # update the sum of squared deviations

  return(obj)
}


#' Generate a proposal parameter vector for MCMC using an adaptive multivariate normal distribution.
#' The covariance matrix of the proposal distribution is updated based on the history of accepted samples.
#' @param Theta a numeric vector of the current parameter values (length should match the number of parameters for the specified method)
#' @param Sigma the covariance matrix for the multivariate normal proposal distribution
MCMC.Proposal.MVN <- function(Theta, Sigma) {
  U <- MCMC.SafeChol(Sigma)
  # this function attempts to ensure that the covariance matrix is positive definite by adding jitter if necessary.
  # It returns the upper triangular Cholesky factor of the covariance matrix,
  # which is then used to generate the proposal by multiplying it with a standard normal random vector and adding it to the current parameter vector Theta.
  z <- rnorm(length(Theta))
  return(as.numeric(Theta + drop(z %*% U)))
}


#' Calculate the covariance matrix for the adaptive MCMC proposal distribution based on the accumulated sum of squared deviations (S) and the number of samples (n).
#' This function returns the covariance matrix by dividing S by (n - 1) to get the sample covariance, and adds a small diagonal matrix to ensure positive definiteness.
#' @param obj a list containing the current state of the covariance estimation, including n (number of samples), mean (current mean vector), and S (current sum of squared deviations matrix)
#' @param eps a small value to add to the diagonal of the covariance matrix to ensure positive definiteness
MCMC.CovGet <- function(obj, eps = 1e-6) {
  NParams <- length(obj$mean)
  if (obj$n < 2L) {
    return(diag(eps, NParams))
  }
  return(obj$S / (obj$n - 1L) + diag(eps, NParams))
}


#' Safely compute the Cholesky decomposition of a covariance matrix, adding jitter if necessary to ensure positive definiteness.
#' This function attempts to compute the Cholesky decomposition of Sigma, adding jitter if necessary to ensure it is positive definite.
#' @param Sigma the covariance matrix to decompose
#' @param jitter the initial amount of jitter to add to the diagonal if Sigma is not positive definite
#' @param max_try the maximum number of attempts to add jitter and compute the Cholesky decomposition before giving up
MCMC.SafeChol <- function(Sigma, jitter = 1e-8, max_try = 8L) {
  NParams <- nrow(Sigma)
  for (k in 0:max_try) {
    U <- tryCatch(
      chol(Sigma + diag(jitter * 10^k, NParams)),
      error = function(e) NULL
    )
    if (!is.null(U)) {
      return(U)
    }
  }
  stop("Proposal covariance is not positive definite.")
}


# Priors --------------------------------------------------------------------

#' Calculate the log-prior probability of the parameter using a log-normal distribution on the eta values, and also includes priors for the auxiliary parameters based on the specified method.
#' @param theta vector of parameters (free eta parameters followed by auxiliary parameters)
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
#' @param PriorArg a list of additional arguments for the prior calculation, such as standard deviations and means for the priors, including:
#'  - sd_comp: standard deviation for the log-normal prior on the comp parameters
#'  - mu_logit_rho: mean for the logit-transformed rho parameter in the Beta-Binomial method
#'  - sd_logit_rho: standard deviation for the logit-transformed rho parameter in the Beta-Binomial method
#'  - mu_log_kappa: mean for the log-transformed kappa parameter in the Dirichlet method
#'  - sd_log_kappa: standard deviation for the log-transformed kappa parameter in the Dirichlet method
#'  - mu_log_sigma_ratio: mean for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - sd_log_sigma_ratio: standard deviation for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - mu_log_sigma_abs: mean for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - sd_log_sigma_abs: standard deviation for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - mu_log_sigma: mean for the log-transformed sigma parameter in the LogDiff method
#'  - sd_log_sigma: standard deviation for the log-transformed sigma parameter in the LogDiff method
MCMC.LogPrior.Theta <- function(
  theta,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  PriorArg = list()
) {
  Method <- match.arg(Method)

  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)
  if (length(theta) != spec$n_par) {
    return(-Inf)
  }
  if (any(!is.finite(theta))) {
    return(-Inf)
  }

  sd_comp <- PriorArg$sd_comp %||% 1
  eta_full <- Comp.EtaFull(theta[seq_len(n_virus - 1L)], n_virus = n_virus)
  lp <- sum(dnorm(eta_full, mean = 0, sd = sd_comp, log = TRUE))

  idx <- n_virus
  if (Method == "BetaBinomial") {
    mu <- PriorArg$mu_logit_rho %||% qlogis(0.02)
    sd <- PriorArg$sd_logit_rho %||% 1.5
    lp <- lp + dnorm(theta[idx], mean = mu, sd = sd, log = TRUE)
  } else if (Method == "Dirichlet") {
    mu <- PriorArg$mu_log_kappa %||% log(100)
    sd <- PriorArg$sd_log_kappa %||% 1
    lp <- lp + dnorm(theta[idx], mean = mu, sd = sd, log = TRUE)
  } else if (Method == "RatioAbs") {
    mu1 <- PriorArg$mu_log_sigma_ratio %||% log(0.7)
    sd1 <- PriorArg$sd_log_sigma_ratio %||% 1
    mu2 <- PriorArg$mu_log_sigma_abs %||% log(0.7)
    sd2 <- PriorArg$sd_log_sigma_abs %||% 1

    lp <- lp +
      dnorm(theta[idx], mean = mu1, sd = sd1, log = TRUE) +
      dnorm(theta[idx + 1L], mean = mu2, sd = sd2, log = TRUE)
  } else if (Method == "LogDiff") {
    mu <- PriorArg$mu_log_sigma %||% log(0.7)
    sd <- PriorArg$sd_log_sigma %||% 1
    lp <- lp + dnorm(theta[idx], mean = mu, sd = sd, log = TRUE)
  }

  return(lp)
}


# Evaluate one theta --------------------------------------------------------

#' Evaluate the log-likelihood, log-prior, and log-posterior for a given parameter vector theta.
#' This function decodes the theta vector into the model parameters, runs the ODE simulation, calculates the log-likelihood based on the specified method, and then calculates the log-prior.
#' @param theta vector of parameters to evaluate
#' @param BaseParm base parameters for the ODE simulation
#' @param Prep preprocessed data for the likelihood calculation
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param LLHArg additional arguments for the likelihood calculation
#' @param PriorArg additional arguments for the prior calculation
#' @return a list containing the log-likelihood (llh), log-prior (lpr), and log-posterior (lpo) for the given theta
MCMC.EvaluateTheta <- function(
  theta,
  BaseParm,
  Prep,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  LLHArg = list(),
  PriorArg = list()
) {
  Method <- match.arg(Method)
  n_virus <- length(BaseParm[["beta0"]])

  decoded <- Inference.DecodeTheta(
    theta = theta,
    Method = Method,
    n_virus = n_virus
  )

  if (any(!is.finite(decoded$comp)) || any(decoded$comp <= 0)) {
    return(list(llh = -Inf, lpr = -Inf, lpo = -Inf))
  }

  Parm <- BaseParm
  Parm[["comp"]] <- decoded$comp

  LLHArg_cur <- LLHArg
  if (Method == "BetaBinomial") {
    LLHArg_cur$rho <- decoded$rho
  } else if (Method == "Dirichlet") {
    LLHArg_cur$kappa <- decoded$kappa
  } else if (Method == "RatioAbs") {
    LLHArg_cur$sigma_ratio <- decoded$sigma_ratio
    LLHArg_cur$sigma_abs <- decoded$sigma_abs
  } else if (Method == "LogDiff") {
    LLHArg_cur$sigma <- decoded$sigma
  }

  llh <-
    Model.RunSim.LLH(
      Parm = Parm,
      Prep = Prep,
      Method = Method,
      LLHArg = LLHArg_cur
    )

  lpr <- MCMC.LogPrior.Theta(
    theta = theta,
    Method = Method,
    n_virus = n_virus,
    PriorArg = PriorArg
  )

  lpo <- if (is.finite(llh) && is.finite(lpr)) llh + lpr else -Inf

  return(list(
    llh = llh,
    lpr = lpr,
    lpo = lpo
  ))
}


# Metropolis-Hastings Markov Chain Monte Carlo ----------------------------------------------

#' Run Metropolis-Hastings MCMC to sample from the posterior distribution of parameters.
#' @param Initial initial theta vector or list of parameter values (comp / rho / kappa / sigma)
#' @param n_iterations number of MCMC iterations
#' @param step proposal step size (scalar or list with comp and aux)
#' @param TargetDat target data for likelihood calculation (if Prep not provided)
#' @param Prep preprocessed data and mappings for likelihood calculation (if TargetDat not provided
#' @param Method likelihood method to use ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param BaseParm base parameters for the model simulation
#' @param after date to define the start of the simulation window
#' @param LLHArg list of additional arguments for the log-likelihood calculation, including:
#'   - viruses: vector of virus names to include in the likelihood calculation
#'   - proxy: metric to use for the likelihood calculation ("Inc_sum" or "I_mean")
#'   - epsilon: small value to avoid numerical issues
#'   - rho: overdispersion parameter for the Beta-Binomial method
#'   - pseudo: pseudo-count for calculating observed probabilities in RatioAbs and LogDiff methods
#'   - sigma_ratio: standard deviation for the ratio component in RatioAbs method
#'   - sigma_abs: standard deviation for the absolute component in RatioAbs method
#'   - sigma: standard deviation for the LogDiff method
#'   - transform: transformation to apply in LogDiff method ("logit" or "log")
#'   - kappa: concentration parameter for the Dirichlet method
#'   - c: small value to avoid numerical issues in the Dirichlet method
#'   - weight: weighting method for the Dirichlet method ("sum_y_capN", "N", "none")
#' @param PriorArg a list of additional arguments for the prior calculation, such as standard deviations and means for the priors, including:
#'  - sd_comp: standard deviation for the log-normal prior on the comp parameters
#'  - mu_logit_rho: mean for the logit-transformed rho parameter in the Beta-Binomial method
#'  - sd_logit_rho: standard deviation for the logit-transformed rho parameter in the Beta-Binomial method
#'  - mu_log_kappa: mean for the log-transformed kappa parameter in the Dirichlet method
#'  - sd_log_kappa: standard deviation for the log-transformed kappa parameter in the Dirichlet method
#'  - mu_log_sigma_ratio: mean for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - sd_log_sigma_ratio: standard deviation for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - mu_log_sigma_abs: mean for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - sd_log_sigma_abs: standard deviation for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - mu_log_sigma: mean for the log-transformed sigma parameter in the LogDiff method
#'  - sd_log_sigma: standard deviation for the log-transformed sigma parameter in the LogDiff method
#' @param verbose whether to show a progress bar
#' @return a matrix of MCMC samples with attributes for log-likelihood, log-prior, log-posterior, acceptance, and method details.
#' The parameter vector theta is method-specific and includes both the compositional parameters (via eta) and any auxiliary parameters (e.g., rho for BetaBinomial, kappa for Dirichlet).
#' The chain is returned on the transformed scale (eta for comp, logit for rho, log for kappa/sigma).
#' Users can decode the chain back to natural parameters using MCMC.DecodeChain and Inference.DecodeTheta.
MCMC.MH <- function(
  Initial = NULL,
  n_iterations,
  step = 0.1,
  TargetDat = NULL,
  Prep = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  BaseParm = Parameter.Create(),
  after = as.Date("2015-08-31"),
  LLHArg = list(),
  PriorArg = list(),
  verbose = TRUE
) {
  Method <- match.arg(Method)
  n_virus <- length(BaseParm[["beta0"]])

  if (is.null(Prep)) {
    if (is.null(TargetDat)) {
      stop("Provide either Prep or TargetDat.")
    }
    Prep <- Inference.Setup(
      TargetDat = TargetDat,
      BaseParm = BaseParm,
      after = after,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
    )
  }

  # spec defines the structure of the parameter vector theta based on the method and number of viruses
  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)
  theta_init <- MCMC.MakeInitial(
    Initial = Initial,
    Method = Method,
    n_virus = n_virus
  ) # initial theta vector on the transformed scale

  chain <- matrix(NA_real_, nrow = n_iterations, ncol = spec$n_par)
  colnames(chain) <- spec$par_names

  llh_trace <- rep(NA_real_, n_iterations)
  lpr_trace <- rep(NA_real_, n_iterations)
  lpo_trace <- rep(NA_real_, n_iterations)
  accepted <- rep(FALSE, n_iterations)

  eval0 <- MCMC.EvaluateTheta(
    theta = theta_init,
    BaseParm = BaseParm,
    Prep = Prep,
    Method = Method,
    LLHArg = LLHArg,
    PriorArg = PriorArg
  )

  if (!is.finite(eval0$lpo)) {
    stop(
      "Initial log-posterior is not finite. Please adjust Initial / LLHArg / PriorArg."
    )
  }

  chain[1, ] <- theta_init
  llh_trace[1] <- eval0$llh
  lpr_trace[1] <- eval0$lpr
  lpo_trace[1] <- eval0$lpo

  if (verbose) {
    pb <- progress::progress_bar$new(
      total = n_iterations,
      clear = TRUE,
      format = "  [:bar] :percent :etas"
    )
    pb$tick()
  }

  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal(
      Theta = chain[i - 1, ],
      Method = Method,
      n_virus = n_virus,
      step = step
    )

    eval_prop <- MCMC.EvaluateTheta(
      theta = proposal,
      BaseParm = BaseParm,
      Prep = Prep,
      Method = Method,
      LLHArg = LLHArg,
      PriorArg = PriorArg
    )

    log_alpha <- eval_prop$lpo - lpo_trace[i - 1]

    if (is.finite(eval_prop$lpo) && log(runif(1)) < log_alpha) {
      chain[i, ] <- proposal
      llh_trace[i] <- eval_prop$llh
      lpr_trace[i] <- eval_prop$lpr
      lpo_trace[i] <- eval_prop$lpo
      accepted[i] <- TRUE
    } else {
      chain[i, ] <- chain[i - 1, ]
      llh_trace[i] <- llh_trace[i - 1]
      lpr_trace[i] <- lpr_trace[i - 1]
      lpo_trace[i] <- lpo_trace[i - 1]
    }

    if (verbose) pb$tick()
  }

  attr(chain, "LLH") <- llh_trace
  attr(chain, "LogPrior") <- lpr_trace
  attr(chain, "LogPost") <- lpo_trace
  attr(chain, "Accepted") <- accepted
  attr(chain, "AcceptanceRate") <- mean(accepted[-1])
  attr(chain, "Method") <- Method
  attr(chain, "LLHArg") <- LLHArg
  attr(chain, "PriorArg") <- PriorArg
  attr(chain, "Spec") <- spec

  return(chain)
}


#' Run adaptive Metropolis-Hastings MCMC to sample from the posterior distribution of parameters, where the proposal covariance is adapted based on the history of the chain.
#' @param Initial initial theta vector or list of parameter values (comp / rho / kappa / sigma)
#' @param n_iterations number of MCMC iterations
#' @param step initial proposal step size (scalar or list with comp and aux)
#' @param TargetDat target data for likelihood calculation (if Prep not provided)
#' @param Prep preprocessed data and mappings for likelihood calculation (if TargetDat not provided)
#' @param Method likelihood method to use ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param BaseParm base parameters for the model simulation
#' @param after date to define the start of the simulation window
#' @param LLHArg list of additional arguments for the log-likelihood calculation, including:
#'   - viruses: vector of virus names to include in the likelihood calculation
#'   - proxy: metric to use for the likelihood calculation ("Inc_sum" or "I_mean")
#'   - epsilon: small value to avoid numerical issues
#'   - rho: overdispersion parameter for the Beta-Binomial method
#'   - pseudo: pseudo-count for calculating observed probabilities in RatioAbs and LogDiff methods
#'   - sigma_ratio: standard deviation for the ratio component in RatioAbs method
#'   - sigma_abs: standard deviation for the absolute component in RatioAbs method
#'   - sigma: standard deviation for the LogDiff method
#'   - transform: transformation to apply in LogDiff method ("logit" or "log")
#'   - kappa: concentration parameter for the Dirichlet method
#'   - c: small value to avoid numerical issues in the Dirichlet method
#'   - weight: weighting method for the Dirichlet method ("sum_y_capN", "N", "none")
#' @param PriorArg a list of additional arguments for the prior calculation, such as standard deviations and means for the priors, including:
#'  - sd_comp: standard deviation for the log-normal prior on the comp parameters
#'  - mu_logit_rho: mean for the logit-transformed rho parameter in the Beta-Binomial method
#'  - sd_logit_rho: standard deviation for the logit-transformed rho parameter in the Beta-Binomial method
#'  - mu_log_kappa: mean for the log-transformed kappa parameter in the Dirichlet method
#'  - sd_log_kappa: standard deviation for the log-transformed kappa parameter in the Dirichlet method
#'  - mu_log_sigma_ratio: mean for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - sd_log_sigma_ratio: standard deviation for the log-transformed sigma_ratio parameter in the RatioAbs method
#'  - mu_log_sigma_abs: mean for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - sd_log_sigma_abs: standard deviation for the log-transformed sigma_abs parameter in the RatioAbs method
#'  - mu_log_sigma: mean for the log-transformed sigma parameter in the LogDiff method
#'  - sd_log_sigma: standard deviation for the log-transformed sigma parameter in the LogDiff method
#' @param Adapt a list of arguments for the adaptive proposal, including:
#'  - use: whether to use adaptive proposal
#'  - start: iteration to start adapting the proposal covariance
#'  - end: iteration to stop adapting the proposal covariance
#'  - every: how often to update the proposal covariance
#'  - target_accept: target acceptance rate for the adaptive proposal (default is 0.44 for 1D and 0.234 for higher dimensions)
#'  - gamma: exponent for the adaptation step size (default is 0.6)
#'  - eps_cov: small value to add to the diagonal of the covariance matrix to ensure positive definiteness (default is 1e-6)
#' @param verbose whether to show a progress bar
#' @return a matrix of MCMC samples with attributes for log-likelihood, log-prior, log-posterior, acceptance, and method details.
#' The parameter vector theta is on the transformed scale (eta for comp, logit for rho, log for kappa/sigma).
#' Users can decode the chain back to natural parameters using MCMC.DecodeChain and Inference.DecodeTheta.
MCMC.MH.Adaptive <- function(
  Initial = NULL,
  n_iterations,
  step = 0.1,
  TargetDat = NULL,
  Prep = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  BaseParm = Parameter.Create(),
  after = as.Date("2015-08-31"),
  LLHArg = list(),
  PriorArg = list(),
  Adapt = list(
    use = TRUE,
    start = 100L,
    end = NULL,
    every = 25L,
    target_accept = NULL,
    gamma = 0.6,
    eps_cov = 1e-6
  ),
  verbose = TRUE
) {
  Method <- match.arg(Method)
  n_virus <- length(BaseParm[["beta0"]])

  if (is.null(Prep)) {
    if (is.null(TargetDat)) {
      stop("Provide either Prep or TargetDat.")
    }
    Prep <- Inference.Setup(
      TargetDat = TargetDat,
      BaseParm = BaseParm,
      after = after,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
    )
  }

  # spec defines the structure of the parameter vector theta based on the method and number of viruses
  spec <- Inference.ParamSpec(Method = Method, n_virus = n_virus)
  theta_init <- MCMC.MakeInitial(
    Initial = Initial,
    Method = Method,
    n_virus = n_virus
  ) # initial theta vector on the transformed scale

  chain <- matrix(NA_real_, nrow = n_iterations, ncol = spec$n_par)
  colnames(chain) <- spec$par_names

  llh_trace <- rep(NA_real_, n_iterations)
  lpr_trace <- rep(NA_real_, n_iterations)
  lpo_trace <- rep(NA_real_, n_iterations)
  accepted <- rep(FALSE, n_iterations)

  eval0 <- MCMC.EvaluateTheta(
    theta = theta_init,
    BaseParm = BaseParm,
    Prep = Prep,
    Method = Method,
    LLHArg = LLHArg,
    PriorArg = PriorArg
  )

  if (!is.finite(eval0$lpo)) {
    stop("Initial log-posterior is not finite.")
  }

  chain[1, ] <- theta_init
  llh_trace[1] <- eval0$llh
  lpr_trace[1] <- eval0$lpr
  lpo_trace[1] <- eval0$lpo

  step_vec <- MCMC.ResolveStepVec(
    step = step,
    Method = Method,
    n_virus = n_virus
  )

  Sigma_prop <- diag(step_vec^2, nrow = spec$n_par, ncol = spec$n_par)

  adapt_use <- isTRUE(Adapt$use)
  adapt_start <- Adapt$start %||% 100L
  adapt_end <- Adapt$end %||% floor(0.5 * n_iterations)
  adapt_every <- Adapt$every %||% 25L
  gamma_exp <- Adapt$gamma %||% 0.6
  eps_cov <- Adapt$eps_cov %||% 1e-6
  target_accept <- Adapt$target_accept %||%
    if (spec$n_par == 1L) 0.44 else 0.234

  log_lambda <- 0
  scale_hist <- rep(NA_real_, n_iterations)
  scale_hist[1] <- exp(log_lambda)

  cov_state <- MCMC.CovInit(spec$n_par)
  cov_state <- MCMC.CovUpdate(cov_state, theta_init)

  acc_batch <- 0L

  if (verbose) {
    pb <- progress::progress_bar$new(
      total = n_iterations,
      clear = TRUE,
      format = "  [:bar] :percent :etas"
    )
    pb$tick()
  }

  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal.MVN(
      Theta = chain[i - 1, ],
      Sigma = Sigma_prop
    )

    eval_prop <- MCMC.EvaluateTheta(
      theta = proposal,
      BaseParm = BaseParm,
      Prep = Prep,
      Method = Method,
      LLHArg = LLHArg,
      PriorArg = PriorArg
    )

    log_alpha <- eval_prop$lpo - lpo_trace[i - 1]

    if (is.finite(eval_prop$lpo) && log(runif(1)) < log_alpha) {
      chain[i, ] <- proposal
      llh_trace[i] <- eval_prop$llh
      lpr_trace[i] <- eval_prop$lpr
      lpo_trace[i] <- eval_prop$lpo
      accepted[i] <- TRUE
    } else {
      chain[i, ] <- chain[i - 1, ]
      llh_trace[i] <- llh_trace[i - 1]
      lpr_trace[i] <- lpr_trace[i - 1]
      lpo_trace[i] <- lpo_trace[i - 1]
    }

    # Update covariance state with the current sample (regardless of acceptance) to ensure the proposal distribution adapts based on the history of the chain.
    cov_state <- MCMC.CovUpdate(cov_state, chain[i, ])

    if (adapt_use && i >= adapt_start && i <= adapt_end) {
      acc_batch <- acc_batch + accepted[i]

      if (((i - adapt_start + 1L) %% adapt_every) == 0L) {
        batch_id <- (i - adapt_start + 1L) / adapt_every
        gamma_i <- 1 / ((batch_id + 10)^gamma_exp) # diminishing adaptation factor

        acc_rate_batch <- acc_batch / adapt_every
        log_lambda <- log_lambda + gamma_i * (acc_rate_batch - target_accept)

        emp_cov <- MCMC.CovGet(cov_state, eps = eps_cov)

        Sigma_prop <- exp(2 * log_lambda) *
          ((2.38^2 / spec$n_par) * emp_cov + diag(eps_cov, spec$n_par))

        acc_batch <- 0L
      }
    }

    scale_hist[i] <- exp(log_lambda)

    if (verbose) {
      pb$tick()
    }
  }

  attr(chain, "LLH") <- llh_trace
  attr(chain, "LogPrior") <- lpr_trace
  attr(chain, "LogPost") <- lpo_trace
  attr(chain, "Accepted") <- accepted
  attr(chain, "AcceptanceRate") <- mean(accepted[-1])
  attr(chain, "Method") <- Method
  attr(chain, "LLHArg") <- LLHArg
  attr(chain, "PriorArg") <- PriorArg
  attr(chain, "Spec") <- spec
  attr(chain, "Adaptive") <- list(
    use = adapt_use,
    start = adapt_start,
    end = adapt_end,
    every = adapt_every,
    target_accept = target_accept,
    final_scale = exp(log_lambda),
    final_cov = Sigma_prop,
    scale_history = scale_hist
  )

  return(chain)
}


# Post-processing ----------------------------------------------------------------

#' Decode the MCMC chain from the transformed parameter space back to the natural parameter space, including the compositional parameters (comp) and any auxiliary parameters (e.g., rho for BetaBinomial, kappa for Dirichlet).
#' @param chain the MCMC chain matrix with parameters on the transformed scale (eta for comp, logit for rho, log for kappa/sigma)
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses (optional, inferred from chain if not provided)
#' @param include_eta whether to include the eta parameters in the output (default is FALSE, which returns only the decoded comp and auxiliary parameters)
MCMC.DecodeChain <- function(
  chain,
  Method = attr(chain, "Method"),
  n_virus = NULL,
  include_eta = FALSE
) {
  chain <- as.matrix(chain)

  if (is.null(n_virus)) {
    spec <- attr(chain, "Spec")
    if (!is.null(spec)) {
      n_virus <- spec$n_virus
    } else {
      stop("Please provide n_virus.")
    }
  }

  eta_free <- chain[, seq_len(n_virus - 1L), drop = FALSE]
  eta_full <- cbind(eta_free, -rowSums(eta_free))
  colnames(eta_full) <- paste0("eta_", seq_len(n_virus))

  comp_full <- exp(eta_full)
  colnames(comp_full) <- paste0("comp_", seq_len(n_virus))

  aux <- switch(
    Method,
    BetaBinomial = matrix(
      plogis(chain[, n_virus]),
      ncol = 1,
      dimnames = list(NULL, "rho")
    ),
    Dirichlet = matrix(
      exp(chain[, n_virus]),
      ncol = 1,
      dimnames = list(NULL, "kappa")
    ),
    RatioAbs = cbind(
      sigma_ratio = exp(chain[, n_virus]),
      sigma_abs = exp(chain[, n_virus + 1L])
    ),
    LogDiff = matrix(
      exp(chain[, n_virus]),
      ncol = 1,
      dimnames = list(NULL, "sigma")
    )
  )

  if (include_eta) {
    return(cbind(eta_full, comp_full, aux))
  } else {
    return(cbind(comp_full, aux))
  }
}


# Parallelization ----------------------------------------------------------------

MCMC.RunChains.Future <- function(
  n_chain = 4L,
  InitialList = NULL,
  n_iterations,
  step = 0.1,
  TargetDat = NULL,
  Prep = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  BaseParm = Parameter.Create(),
  after = as.Date("2015-08-31"),
  LLHArg = list(),
  PriorArg = list(),
  Adapt = list(),
  workers = min(n_chain, future::availableCores()),
  seed = 123,
  source_r = NULL,
  source_cpp = NULL,
  wait = TRUE
) {
  Method <- match.arg(Method)
  n_virus <- length(BaseParm[["beta0"]])

  if (is.null(Prep)) {
    if (is.null(TargetDat)) {
      stop("Provide either Prep or TargetDat.")
    }
    Prep <- Inference.Setup(
      TargetDat = TargetDat,
      BaseParm = BaseParm,
      after = after,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
    )
  }

  if (is.null(InitialList)) {
    InitialList <- MCMC.MakeInitialList(
      n_chain = n_chain,
      Initial = NULL,
      Method = Method,
      n_virus = n_virus,
      jitter = 0.3,
      seed = seed
    )
  }
  stopifnot(length(InitialList) == n_chain)

  old_plan <- plan()
  plan(multisession, workers = workers)

  # 仅同步模式下自动恢复 plan
  if (wait) {
    on.exit(plan(old_plan), add = TRUE)
  }

  fs <- lapply(seq_len(n_chain), function(k) {
    future(
      {
        if (!is.null(source_r)) {
          source(source_r)
        }

        if (!is.null(source_cpp)) {
          # 避免并发 sourceCpp 冲突：每个 worker 用独立 cache 目录
          cache_dir <- file.path(tempdir(), paste0("rcpp_cache_", Sys.getpid()))
          dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
          sourceCpp(source_cpp, cacheDir = cache_dir, rebuild = FALSE)
        }

        MCMC.MH.Adaptive(
          Initial = InitialList[[k]],
          n_iterations = n_iterations,
          step = step,
          Prep = Prep,
          Method = Method,
          BaseParm = BaseParm,
          after = after,
          LLHArg = LLHArg,
          PriorArg = PriorArg,
          Adapt = Adapt,
          verbose = FALSE
        )
      },
      seed = seed + k
    )
  })
  names(fs) <- paste0("chain_", seq_len(n_chain))

  if (wait) {
    chains <- lapply(fs, future::value)
    names(chains) <- names(fs)
    return(chains)
  }

  structure(
    list(futures = fs, old_plan = old_plan),
    class = "mcmc_future_job"
  )
}
