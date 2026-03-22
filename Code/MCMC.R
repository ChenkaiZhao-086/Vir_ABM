Utility.Clamp01 <- function(p, epsilon = 1e-9) {
  pmin(pmax(p, epsilon), 1 - epsilon)
}

Utility.Logit <- function(p) log(p / (1 - p))

`%||%` <- function(x, y) if (is.null(x)) y else x


# Fast preparation of parameters for the ODE simulation -----------------------------------

#' Aggregate target data into matrices indexed by ISO week and virus.
#' @param TargetDat a data.table with columns Virus, ISOweek, y, N
#' @param viruses   character vector of virus names to include (NULL = all)
#' @return list: ISOweek, viruses, y (matrix), N (matrix), N_week, n_virus
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
    .(N = if (all(is.na(N))) NA_real_ else max(N, na.rm = TRUE)),
    by = .(ISOweek, Virus)
  ]

  y_wide <- dcast(DT_y, ISOweek ~ Virus, value.var = "y", fill = NA_real_)
  N_wide <- dcast(DT_N, ISOweek ~ Virus, value.var = "N", fill = NA_real_)

  for (v in viruses) {
    if (!v %in% names(y_wide)) {
      y_wide[[v]] <- NA_real_
    }
    if (!v %in% names(N_wide)) N_wide[[v]] <- NA_real_
  }

  setcolorder(y_wide, c("ISOweek", viruses))
  setcolorder(N_wide, c("ISOweek", viruses))

  y_wide <- y_wide[order(y_wide$ISOweek)]
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

  list(
    ISOweek = y_wide$ISOweek,
    viruses = viruses,
    y = y_mat,
    N = N_mat,
    N_week = N_week,
    n_virus = length(viruses)
  )
}


#' Build a date-to-ISO-week mapping for the simulation period.
#' @param Parm  parameter list with year_start and year_end
#' @param after only include dates on or after this date
#' @param before only include dates on or before this date
#' @return list: times, keep, group, ISOweek, n_day_in_week
Model.MakeWeekMap <- function(
  Parm,
  after = as.Date("2015-08-31"),
  before = as.Date("2020-01-06")
) {
  times <- as.numeric(seq(
    from = as.Date(Parm[["year_start"]]),
    to = as.Date(Parm[["year_end"]]),
    by = 1
  ))
  dates <- as.Date(times, origin = "1970-01-01")
  monday <- dates - (as.integer(strftime(dates, "%u")) - 1L)

  keep <- dates >= after & dates <= before
  monday_keep <- monday[keep]
  week_levels <- unique(monday_keep)
  group <- factor(monday_keep, levels = week_levels, ordered = TRUE)

  list(
    times = times,
    keep = keep,
    group = group,
    ISOweek = strftime(week_levels, "%G-W%V"),
    n_day_in_week = as.integer(table(group))
  )
}


#' Prepare all data structures needed for MCMC inference.
#' @param TargetDat data.table of target observations
#' @param BaseParm  base parameter list for the ODE model
#' @param after     only use data after this date
#' @param viruses   virus names to include (NULL defaults to Virus_1, Virus_2, ...)
#' @return list: target, week_map, sim_index, init_state, after
Inference.Setup <- function(
  TargetDat,
  BaseParm,
  after = as.Date("2015-08-31"),
  before = as.Date("2020-01-06"),
  viruses = NULL
) {
  if (is.null(viruses)) {
    viruses <- paste0("Virus_", seq_along(BaseParm[["beta0"]]))
  }

  target <- Target.Prepare(TargetDat, viruses = viruses)
  week_map <- Model.MakeWeekMap(BaseParm, after = after, before = before)

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

  list(
    target = target,
    week_map = week_map,
    sim_index = sim_index[keep],
    init_state = init_state,
    after = after,
    before = before
  )
}


#' Resolve and validate inference options.
#' @param InferArg list of inference arguments
#' @param n_virus  number of viruses
Inference.ResolveInferArg <- function(InferArg = list(), n_virus = 4L) {
  out <- modifyList(
    list(
      infer_penal = FALSE,
      penal_center_with_comp = FALSE,
      infer_NPISes = FALSE,
      NPISes_shared = FALSE
    ),
    InferArg
  )

  out$infer_penal <- isTRUE(out$infer_penal)
  out$penal_center_with_comp <- isTRUE(out$penal_center_with_comp) &&
    out$infer_penal
  out$infer_NPISes <- isTRUE(out$infer_NPISes)
  out$NPISes_shared <- isTRUE(out$NPISes_shared)

  out$n_npises <- if (!out$infer_NPISes) {
    0L
  } else if (out$NPISes_shared) {
    1L
  } else {
    as.integer(n_virus)
  }

  out
}


#' Run the ODE and aggregate daily output to weekly incidence or prevalence.
#' Simulation columns are always named Virus_1, Virus_2, ... in order.
#' @param Parm      ODE parameter list
#' @param WeekMap   output of Model.MakeWeekMap
#' @param metric    "Inc_sum" (weekly incidence) or "I_mean" (mean prevalence)
#' @param InitState initial state vector (auto-generated if NULL or wrong length)
#' @return list: ISOweek (character), mat (matrix with columns Virus_1, ...)
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

  sim <- as.matrix(ode(
    y = InitState,
    times = WeekMap$times,
    func = ParmInferenceCpp,
    parms = Parm,
    method = "rk4"
  ))

  if (metric == "Inc_sum") {
    cols <- grep("^Inc_\\d+$", colnames(sim))
    if (length(cols) == 0L) {
      stop("No incidence columns matching '^Inc_\\d+$' found in ODE output.")
    }
    daily <- sim[WeekMap$keep, cols, drop = FALSE]
    weekly <- rowsum(daily, group = WeekMap$group, reorder = FALSE)
  } else {
    cols <- grep("^I_\\d+$", colnames(sim))
    if (length(cols) == 0L) {
      stop("No state columns matching '^I_\\d+$' found in ODE output.")
    }
    daily <- sim[WeekMap$keep, cols, drop = FALSE]
    weekly <- rowsum(daily, group = WeekMap$group, reorder = FALSE)
    weekly <- sweep(weekly, 1, WeekMap$n_day_in_week, "/")
  }

  colnames(weekly) <- paste0("Virus_", seq_len(ncol(weekly)))
  list(ISOweek = WeekMap$ISOweek, mat = weekly)
}


# Likelihoods ---------------------------------------------------------------

#' Beta-Binomial log-PMF. Inputs are assumed pre-validated by the caller.
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
  lchoose(n, y) + lbeta(y + a, n - y + b) - lbeta(a, b)
}


#' Aggregate Beta-Binomial log-likelihood over all valid matrix entries.
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

  sum(
    LLH.BetaBinomial.PMF(
      y = y_mat[valid],
      n = N_mat[valid],
      p = p_sim_mat[valid],
      rho = rho,
      epsilon = epsilon
    ),
    na.rm = TRUE
  )
}


#' Log-likelihood combining absolute incidence (logit scale, Virus_1 as reference)
#' with log-ratios of all other viruses to the reference.
#' @param p_obs_mat matrix of observed probabilities (incidence proportions)
#' @param p_sim_mat matrix of simulated probabilities
#' @param VirNames vector of virus names corresponding to the columns of the matrices
#' @param sigma_ratio standard deviation for the ratio component of the likelihood
#' @param sigma_abs standard deviation for the absolute component of the likelihood
#' @param epsilon small value to avoid numerical issues
LLH.RatioAbs <- function(
  p_obs_mat,
  p_sim_mat,
  viruses = NULL,
  sigma_ratio = 0.7,
  sigma_abs = 0.7,
  epsilon = 1e-9
) {
  p_obs_mat <- Utility.Clamp01(p_obs_mat, epsilon)
  p_sim_mat <- Utility.Clamp01(p_sim_mat, epsilon)

  if (is.null(viruses)) {
    viruses <- paste0("Virus_", seq_len(ncol(p_obs_mat)))
  }

  ref_idx <- match("Virus_1", viruses)
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

  LLH_Abs <- sum(
    dnorm(
      Utility.Logit(Ref_Obs),
      mean = Utility.Logit(Ref_Sim),
      sd = sigma_abs,
      log = TRUE
    ),
    na.rm = TRUE
  )

  LLH_Ratios <- 0
  for (j in setdiff(seq_along(viruses), ref_idx)) {
    pObs <- p_obs_mat[, j]
    pSim <- p_sim_mat[, j]
    ok <- is.finite(pObs) & is.finite(pSim)
    if (!any(ok)) {
      next
    }
    LLH_Ratios <- LLH_Ratios +
      sum(
        dnorm(
          log(pObs[ok]) - log(Ref_Obs[ok]),
          mean = log(pSim[ok]) - log(Ref_Sim[ok]),
          sd = sigma_ratio,
          log = TRUE
        ),
        na.rm = TRUE
      )
  }

  LLH_Abs + LLH_Ratios
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

  sum(dnorm(x_obs, mean = x_sim, sd = sigma, log = TRUE), na.rm = TRUE)
}


#'
LLH.RollingMaxPoisson <- function(
  y_mat,
  sim_mat,
  roll_n = 52L,
  roll_align = "center",
  lambda_floor = 0
) {
  Obs <- as.matrix(y_mat)
  Sim <- as.matrix(sim_mat)

  if (!all(dim(Obs) == dim(Sim))) {
    stop("y_mat and sim_mat must have the same dimensions.")
  }

  nr <- nrow(Obs)
  nc <- ncol(Obs)

  Sim_roll <- matrix(NA_real_, nr, nc)
  Obs_roll <- matrix(NA_real_, nr, nc)

  for (j in seq_len(nc)) {
    Sim_roll[, j] <- data.table::frollmax(
      as.numeric(Sim[, j]),
      n = roll_n,
      align = roll_align,
      na.rm = TRUE
    )
    Obs_roll[, j] <- data.table::frollmax(
      as.numeric(Obs[, j]),
      n = roll_n,
      align = roll_align,
      na.rm = TRUE
    )
  }

  Prop <- Sim / Sim_roll
  Lambda <- Prop * Obs_roll

  if (is.finite(lambda_floor) && lambda_floor > 0) {
    Lambda <- pmax(Lambda, lambda_floor)
  }

  Obs_use <- round(Obs)

  ok <- is.finite(Obs_use) & is.finite(Lambda) & (Obs_use >= 0) & (Lambda >= 0)
  if (!any(ok)) {
    return(0)
  }

  ll <- matrix(NA_real_, nr, nc)
  ll[ok] <- dpois(x = Obs_use[ok], lambda = Lambda[ok], log = TRUE)
  sum(ll, na.rm = TRUE)
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
#' @param weight method for weighting the concentration parameter
#' - "sum_y_capN" uses the sum of observed counts capped by N_vec,
#' - "N" uses N_vec directly,
#' - "none" uses no weighting
LLH.Dirichlet <- function(
  y_mat,
  sim_mat,
  kappa,
  N_vec = NULL,
  c = 1e-8,
  weight = c("sum_y_capN", "N", "none"),
  add_roll_poisson = TRUE,
  roll_n = 52L,
  roll_align = "center",
  poisson_weight = 1,
  lambda_floor = 0
) {
  weight <- match.arg(weight)
  if (!is.finite(kappa) || kappa <= 0) {
    return(-Inf)
  }

  keep <- apply(is.finite(y_mat), 1, all) &
    apply(is.finite(sim_mat), 1, all) &
    rowSums(y_mat, na.rm = TRUE) > 0 &
    rowSums(sim_mat, na.rm = TRUE) > 0

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
    w <- rowSums(Y)
    if (!is.null(N_vec)) {
      N_use <- N_vec[keep]
      ok <- is.finite(N_use)
      w[ok] <- pmin(w[ok], N_use[ok])
    }
  }

  Alpha <- sweep(S_sim, 1, kappa * pmax(w, 1), "*")
  Alpha <- pmax(Alpha, 1e-12)

  llh_dir <- sum(LLH.Dirichlet.PMF(S_obs, Alpha))

  if (
    !isTRUE(add_roll_poisson) ||
      !is.finite(poisson_weight) ||
      poisson_weight == 0
  ) {
    return(llh_dir)
  }

  llh_roll <- LLH.RollingMaxPoisson(
    y_mat = y_mat,
    sim_mat = sim_mat,
    roll_n = roll_n,
    roll_align = roll_align,
    lambda_floor = lambda_floor
  )

  llh_dir + poisson_weight * llh_roll
}


# Main log-likelihood wrapper -----------------------------------------------

#' Run the ODE and compute the log-likelihood under the chosen method.
#' @param Parm      ODE parameter list
#' @param after     simulation start date (used only when Prep is NULL)
#' @param TargetDat raw target data (used only when Prep is NULL)
#' @param Prep      pre-built inference setup (preferred; avoids re-computation)
#' @param Method    likelihood method: "BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"
#' @param LLHArg    list of additional arguments for the log-likelihood calculation, including:
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
  before = as.Date("2020-01-06"),
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
      before = before,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_along(Parm[["beta0"]]))
    )
  }

  proxy <- LLHArg$proxy %||% "Inc_sum"
  epsilon <- LLHArg$epsilon %||% 1e-9
  target <- Prep$target

  sim <- Model.RunSim.Weekly(
    Parm = Parm,
    WeekMap = Prep$week_map,
    metric = proxy,
    InitState = Prep$init_state
  )

  sim_mat <- sim$mat[Prep$sim_index, target$viruses, drop = FALSE]
  p_sim_mat <- Utility.Clamp01(
    sim_mat / Parm[["num_of_agent"]],
    epsilon = epsilon
  )

  if (Method == "BetaBinomial") {
    return(LLH.BetaBinomial(
      y_mat = target$y,
      N_mat = target$N,
      p_sim_mat = p_sim_mat,
      rho = LLHArg$rho %||% 0.02,
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
    return(LLH.RatioAbs(
      p_obs_mat = p_obs_mat,
      p_sim_mat = p_sim_mat,
      viruses = target$viruses,
      sigma_ratio = LLHArg$sigma_ratio %||% 0.7,
      sigma_abs = LLHArg$sigma_abs %||% 0.7,
      epsilon = epsilon
    ))
  }

  if (Method == "LogDiff") {
    return(LLH.LogDiff(
      p_obs_mat = p_obs_mat,
      p_sim_mat = p_sim_mat,
      transform = LLHArg$transform %||% "logit",
      sigma = LLHArg$sigma %||% 0.7,
      epsilon = epsilon
    ))
  }

  # Dirichlet
  LLH.Dirichlet(
    y_mat = target$y,
    sim_mat = sim_mat,
    kappa = LLHArg$kappa %||% 100,
    N_vec = target$N_week,
    c = LLHArg$c %||% 1e-8,
    weight = LLHArg$weight %||% "sum_y_capN",
    add_roll_poisson = LLHArg$add_roll_poisson %||% TRUE,
    roll_n = LLHArg$roll_n %||% 52L,
    roll_align = LLHArg$roll_align %||% "center",
    poisson_weight = LLHArg$poisson_weight %||% 1,
    lambda_floor = LLHArg$lambda_floor %||% 0
  )
}


# Parameterization of comp --------------------------------------------------

#' Append -sum(eta_free) to form the full K-length centered eta vector.
Comp.EtaFull <- function(eta_free, n_virus = length(eta_free) + 1L) {
  stopifnot(length(eta_free) == n_virus - 1L)
  c(eta_free, -sum(eta_free))
}


#' Log-prior for centered compositional parameters.
#' @param comp full comp vector (will be mean-centred internally), or
#' @param eta  free eta vector (last element appended as -sum)
MCMC.LogPrior.Comp <- function(comp = NULL, eta = NULL, sd = 1) {
  if (!is.null(comp)) {
    comp <- as.numeric(comp)
    if (any(!is.finite(comp))) {
      return(-Inf)
    }
    comp <- comp - mean(comp)
    return(sum(dnorm(comp, mean = 0, sd = sd, log = TRUE)))
  }
  if (is.null(eta)) {
    stop("Provide either comp or eta.")
  }
  eta <- as.numeric(eta)
  if (any(!is.finite(eta))) {
    return(-Inf)
  }
  sum(dnorm(c(eta, -sum(eta)), mean = 0, sd = sd, log = TRUE))
}


# Method-specific theta spec ------------------------------------------------

#' Build the parameter specification for a given inference method.
#' @param Method   inference method
#' @param n_virus  number of viruses
#' @param InferArg: a list of inference arguments that may affect the parameter specification, including:
#'   - infer_penal: whether to include a penalization parameter
#'   - penal_center_with_comp: whether the penalization parameter is centered with the comp parameters
#'   - infer_NPISes: whether to include NPISes parameters
#'   - NPISes_shared: whether the NPISes parameters are shared across viruses
#' @return a list containing the method, number of viruses, component names, auxiliary names, parameter names, and number of parameters
Inference.ParamSpec <- function(
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  InferArg = list()
) {
  Method <- match.arg(Method)
  ia <- Inference.ResolveInferArg(InferArg = InferArg, n_virus = n_virus)

  comp_names <- if (ia$penal_center_with_comp) {
    paste0("cp_free_", seq_len(n_virus))
  } else {
    paste0("comp_free_", seq_len(n_virus - 1L))
  }

  method_aux_names <- switch(
    Method,
    BetaBinomial = "z_rho",
    Dirichlet = "log_kappa",
    RatioAbs = c("log_sigma_ratio", "log_sigma_abs"),
    LogDiff = "log_sigma"
  )

  base_aux_names <- c(
    if (ia$infer_penal && !ia$penal_center_with_comp) "penal" else character(0),
    if (ia$infer_NPISes) {
      if (ia$NPISes_shared) {
        "log_NPISes_shared"
      } else {
        paste0("log_NPISes_", seq_len(n_virus))
      }
    } else {
      character(0)
    }
  )

  par_names <- c(comp_names, method_aux_names, base_aux_names)
  ptr <- length(comp_names)

  idx_method_aux <- if (length(method_aux_names) > 0L) {
    (ptr + 1L):(ptr + length(method_aux_names))
  } else {
    integer(0)
  }
  ptr <- ptr + length(method_aux_names)

  idx_penal <- if (ia$infer_penal && !ia$penal_center_with_comp) {
    ptr + 1L
  } else {
    integer(0)
  }
  ptr <- ptr + length(idx_penal)

  idx_npises <- if (ia$n_npises > 0L) {
    (ptr + 1L):(ptr + ia$n_npises)
  } else {
    integer(0)
  }

  list(
    Method = Method,
    n_virus = n_virus,
    Infer = ia,
    comp_names = comp_names,
    method_aux_names = method_aux_names,
    base_aux_names = base_aux_names,
    aux_names = c(method_aux_names, base_aux_names),
    par_names = par_names,
    n_par = length(par_names),
    idx = list(
      comp = seq_along(comp_names),
      method_aux = idx_method_aux,
      penal = idx_penal,
      npises = idx_npises
    )
  )
}


#' Decode the theta vector into interpretable model parameters.
#' This function takes the free eta parameters, reconstructs the full eta vector, and then calculates the comp values. It also extracts the auxiliary parameters based on the specified method.
#' @param theta vector of parameters to decode
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
#' @return a list containing the full eta vector, comp values, and auxiliary parameters
Inference.DecodeTheta <- function(
  theta,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  InferArg = list()
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )

  if (length(theta) != spec$n_par) {
    stop("theta length does not match method specification.")
  }

  ia <- spec$Infer
  out <- list()
  comp_free <- theta[spec$idx$comp]

  if (ia$penal_center_with_comp) {
    cp_full <- c(comp_free, -sum(comp_free))
    out$comp <- cp_full[seq_len(n_virus)]
    out$Penal <- cp_full[n_virus + 1L]
    out$cp_full <- cp_full
  } else {
    comp_full <- c(comp_free, -sum(comp_free))
    out$comp <- comp_full
    out$eta <- comp_full
    if (length(spec$idx$penal) == 1L) out$Penal <- theta[spec$idx$penal]
  }

  idx_m <- spec$idx$method_aux
  if (Method == "BetaBinomial") {
    out$rho <- plogis(theta[idx_m[1L]])
  } else if (Method == "Dirichlet") {
    out$kappa <- exp(theta[idx_m[1L]])
  } else if (Method == "RatioAbs") {
    out$sigma_ratio <- exp(theta[idx_m[1L]])
    out$sigma_abs <- exp(theta[idx_m[2L]])
  } else if (Method == "LogDiff") {
    out$sigma <- exp(theta[idx_m[1L]])
  }

  if (length(spec$idx$npises) > 0L) {
    z <- theta[spec$idx$npises]
    out$NPISes <- if (ia$NPISes_shared) rep(exp(z[1L]), n_virus) else exp(z)
  }

  out
}

# Initial values ------------------------------------------------------------

#' Build the initial theta vector for MCMC.
#' @param Initial  NULL, a list of named parameter values, or a numeric vector
#' @param Method   inference method
#' @param n_virus  number of viruses
#' @param BaseParm base parameter list (used for defaults)
#' @param InferArg inference options
MCMC.MakeInitial <- function(
  Initial = NULL,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  BaseParm = NULL,
  InferArg = list()
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )
  ia <- spec$Infer

  aux_default <- switch(
    Method,
    BetaBinomial = qlogis(0.02),
    Dirichlet = log(100),
    RatioAbs = c(log(0.7), log(0.7)),
    LogDiff = log(0.7)
  )

  theta0 <- rep(0, spec$n_par)
  names(theta0) <- spec$par_names
  theta0[spec$idx$method_aux] <- aux_default

  if (length(spec$idx$penal) == 1L) {
    theta0[spec$idx$penal] <- BaseParm[["Penal"]] %||% 0
  }

  if (length(spec$idx$npises) > 0L) {
    np0 <- as.numeric(BaseParm[["NPISes"]] %||% rep(1, n_virus))
    if (ia$NPISes_shared) {
      theta0[spec$idx$npises] <- log(mean(np0))
    } else {
      if (length(np0) == 1L) {
        np0 <- rep(np0, n_virus)
      }
      if (length(np0) != n_virus) {
        stop("BaseParm$NPISes length mismatch.")
      }
      theta0[spec$idx$npises] <- log(np0)
    }
  }

  if (is.null(Initial)) {
    return(theta0)
  }

  if (is.list(Initial)) {
    if (!is.null(Initial$theta)) {
      th <- as.numeric(Initial$theta)
      if (length(th) != spec$n_par) {
        stop("Initial$theta has wrong length.")
      }
      return(th)
    }

    if (ia$penal_center_with_comp) {
      if (!is.null(Initial$cp_full)) {
        cp_full <- as.numeric(Initial$cp_full)
        if (length(cp_full) != n_virus + 1L) {
          stop("Initial$cp_full must have length n_virus + 1.")
        }
      } else {
        comp0 <- as.numeric(
          Initial$comp %||% (BaseParm[["comp"]] %||% rep(0, n_virus))
        )
        pen0 <- as.numeric(
          Initial$Penal %||% Initial$penal %||% (BaseParm[["Penal"]] %||% 0)
        )
        if (length(comp0) != n_virus) {
          stop("Initial$comp must have length n_virus.")
        }
        cp_full <- c(comp0, pen0)
      }
      if (any(!is.finite(cp_full))) {
        stop("Initial cp/penal contains non-finite values.")
      }
      cp_full <- cp_full - mean(cp_full)
      theta0[spec$idx$comp] <- cp_full[seq_len(n_virus)]
    } else {
      comp_in <- Initial$comp %||% Initial$eta
      if (!is.null(comp_in)) {
        comp_in <- as.numeric(comp_in)
        if (length(comp_in) == n_virus) {
          comp_in <- comp_in - mean(comp_in)
          theta0[spec$idx$comp] <- comp_in[seq_len(n_virus - 1L)]
        } else if (length(comp_in) == n_virus - 1L) {
          theta0[spec$idx$comp] <- comp_in
        } else {
          stop("Initial$comp / Initial$eta has wrong length.")
        }
      }
      if (length(spec$idx$penal) == 1L) {
        p0 <- Initial$Penal %||% Initial$penal
        if (!is.null(p0)) theta0[spec$idx$penal] <- as.numeric(p0)
      }
    }

    if (Method == "BetaBinomial" && !is.null(Initial$rho)) {
      theta0[spec$idx$method_aux[1L]] <- qlogis(Initial$rho)
    }
    if (Method == "Dirichlet" && !is.null(Initial$kappa)) {
      theta0[spec$idx$method_aux[1L]] <- log(Initial$kappa)
    }
    if (Method == "RatioAbs") {
      if (!is.null(Initial$sigma_ratio)) {
        theta0[spec$idx$method_aux[1L]] <- log(Initial$sigma_ratio)
      }
      if (!is.null(Initial$sigma_abs)) {
        theta0[spec$idx$method_aux[2L]] <- log(Initial$sigma_abs)
      }
    }
    if (Method == "LogDiff" && !is.null(Initial$sigma)) {
      theta0[spec$idx$method_aux[1L]] <- log(Initial$sigma)
    }

    if (length(spec$idx$npises) > 0L && !is.null(Initial$NPISes)) {
      np0 <- as.numeric(Initial$NPISes)
      if (ia$NPISes_shared) {
        theta0[spec$idx$npises] <- log(if (length(np0) > 1L) mean(np0) else np0)
      } else {
        if (length(np0) == 1L) {
          np0 <- rep(np0, n_virus)
        }
        if (length(np0) != n_virus) {
          stop("Initial$NPISes length mismatch.")
        }
        theta0[spec$idx$npises] <- log(np0)
      }
    }

    if (length(spec$idx$npises) > 0L && !is.null(Initial$log_NPISes)) {
      z <- as.numeric(Initial$log_NPISes)
      if (ia$NPISes_shared) {
        theta0[spec$idx$npises] <- if (length(z) > 1L) mean(z) else z
      } else {
        if (length(z) == 1L) {
          z <- rep(z, n_virus)
        }
        if (length(z) != n_virus) {
          stop("Initial$log_NPISes length mismatch.")
        }
        theta0[spec$idx$npises] <- z
      }
    }

    return(theta0)
  }

  # numeric Initial
  Initial <- as.numeric(Initial)
  if (length(Initial) == spec$n_par) {
    return(Initial)
  }
  if (length(Initial) == length(spec$idx$comp)) {
    theta0[spec$idx$comp] <- Initial
    return(theta0)
  }

  if (!ia$penal_center_with_comp && length(Initial) == n_virus) {
    x <- Initial - mean(Initial)
    theta0[spec$idx$comp] <- x[seq_len(n_virus - 1L)]
    return(theta0)
  }
  if (ia$penal_center_with_comp && length(Initial) == n_virus + 1L) {
    x <- Initial - mean(Initial)
    theta0[spec$idx$comp] <- x[seq_len(n_virus)]
    return(theta0)
  }

  stop("Initial has unsupported length.")
}


# Proposal ------------------------------------------------------------------

#' Resolve proposal step sizes into separate comp and aux vectors.
#' @param step a scalar or a list specifying the standard deviation of the proposal distribution for the component parameters and auxiliary parameters. If a scalar is provided, it will be used for all parameters
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
MCMC.ResolveStep <- function(
  step = 0.1,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  InferArg = list()
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )

  n_comp <- length(spec$idx$comp)
  n_aux <- length(spec$aux_names)

  if (is.numeric(step) && length(step) == 1L) {
    return(list(comp = rep(step, n_comp), aux = rep(step, n_aux)))
  }
  if (!is.list(step)) {
    stop("step must be either a scalar or a list.")
  }

  step_comp <- rep(step$comp %||% 0.10, n_comp)
  if (spec$Infer$penal_center_with_comp && !is.null(step$cp)) {
    step_comp <- rep(step$cp, n_comp)
  }

  step_aux <- rep(step$aux %||% 0.10, n_aux)
  for (j in seq_along(spec$aux_names)) {
    nm <- spec$aux_names[j]
    if (!is.null(step[[nm]])) step_aux[j] <- step[[nm]]
  }

  list(comp = step_comp, aux = step_aux)
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
  step = 0.1,
  InferArg = list()
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )
  steps <- MCMC.ResolveStep(
    step = step,
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )

  out <- Theta
  idx_comp <- spec$idx$comp
  if (length(idx_comp) > 0L) {
    out[idx_comp] <- out[idx_comp] + rnorm(length(idx_comp), 0, steps$comp)
  }
  idx_aux <- c(spec$idx$method_aux, spec$idx$penal, spec$idx$npises)
  if (length(idx_aux) > 0L) {
    out[idx_aux] <- out[idx_aux] + rnorm(length(idx_aux), 0, steps$aux)
  }
  out
}


# Adaptative proposal -------------------------------------------------------------

#' Return all step sizes as a single concatenated vector (comp then aux).
#' This function calls MCMC.ResolveStep to get the step sizes for the component parameters and auxiliary parameters, and then concatenates them into a single vector.
#' @param step a scalar or a list specifying the standard deviation of the proposal distribution for the component parameters and auxiliary parameters. If a scalar is provided, it will be used for all parameters
#' @param Method the inference method ("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet")
#' @param n_virus the number of viruses
MCMC.ResolveStepVec <- function(
  step = 0.1,
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  n_virus = 4L,
  InferArg = list()
) {
  st <- MCMC.ResolveStep(
    step = step,
    Method = match.arg(Method),
    n_virus = n_virus,
    InferArg = InferArg
  )
  c(st$comp, st$aux)
}


#' Initialise the Welford online-covariance accumulator.
#' This function creates a list to store the number of samples (n), the current mean vector, and the sum of squared deviations matrix (S) for the parameters.
#' This will be updated iteratively as new samples are accepted during the MCMC process.
#' @param NParams the number of parameters in the model (length of the parameter vector theta)
MCMC.CovInit <- function(n_params) {
  list(n = 0L, mean = numeric(n_params), S = matrix(0, n_params, n_params))
}


#' Update the Welford accumulator with a new sample x.
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
  n_params <- nrow(Sigma)
  for (k in 0:max_try) {
    U <- tryCatch(
      chol(Sigma + diag(jitter * 10^k, n_params)),
      error = function(e) NULL
    )
    if (!is.null(U)) return(U)
  }
  stop("Proposal covariance is not positive definite.")
}


# Priors --------------------------------------------------------------------

#' Log-prior for the full theta vector.
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
  PriorArg = list(),
  InferArg = list()
) {
  Method <- match.arg(Method)
  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )
  ia <- spec$Infer

  if (length(theta) != spec$n_par || any(!is.finite(theta))) {
    return(-Inf)
  }

  lp <- 0

  if (ia$penal_center_with_comp) {
    cp_full <- c(theta[spec$idx$comp], -sum(theta[spec$idx$comp]))
    mu <- PriorArg$mu_comp_penal %||% 0
    sd <- PriorArg$sd_comp_penal %||% 1
    if (length(mu) == 1L) {
      mu <- rep(mu, length(cp_full))
    }
    if (length(sd) == 1L) {
      sd <- rep(sd, length(cp_full))
    }
    if (length(mu) != length(cp_full) || length(sd) != length(cp_full)) {
      stop("mu_comp_penal / sd_comp_penal length mismatch.")
    }
    lp <- lp + sum(dnorm(cp_full, mu, sd, log = TRUE))
  } else {
    comp_full <- c(theta[spec$idx$comp], -sum(theta[spec$idx$comp]))
    mu <- PriorArg$mu_comp %||% 0
    sd <- PriorArg$sd_comp %||% 1
    if (length(mu) == 1L) {
      mu <- rep(mu, length(comp_full))
    }
    if (length(sd) == 1L) {
      sd <- rep(sd, length(comp_full))
    }
    if (length(mu) != length(comp_full) || length(sd) != length(comp_full)) {
      stop("mu_comp / sd_comp length mismatch.")
    }
    lp <- lp + sum(dnorm(comp_full, mu, sd, log = TRUE))

    if (length(spec$idx$penal) == 1L) {
      lp <- lp +
        dnorm(
          theta[spec$idx$penal],
          mean = PriorArg$mu_penal %||% 0,
          sd = PriorArg$sd_penal %||% 1,
          log = TRUE
        )
    }
  }

  idx <- spec$idx$method_aux
  if (Method == "BetaBinomial") {
    lp <- lp +
      dnorm(
        theta[idx[1L]],
        PriorArg$mu_logit_rho %||% qlogis(0.02),
        PriorArg$sd_logit_rho %||% 1.5,
        log = TRUE
      )
  } else if (Method == "Dirichlet") {
    lp <- lp +
      dnorm(
        theta[idx[1L]],
        PriorArg$mu_log_kappa %||% log(100),
        PriorArg$sd_log_kappa %||% 1,
        log = TRUE
      )
  } else if (Method == "RatioAbs") {
    lp <- lp +
      dnorm(
        theta[idx[1L]],
        PriorArg$mu_log_sigma_ratio %||% log(0.7),
        PriorArg$sd_log_sigma_ratio %||% 1,
        log = TRUE
      ) +
      dnorm(
        theta[idx[2L]],
        PriorArg$mu_log_sigma_abs %||% log(0.7),
        PriorArg$sd_log_sigma_abs %||% 1,
        log = TRUE
      )
  } else if (Method == "LogDiff") {
    lp <- lp +
      dnorm(
        theta[idx[1L]],
        PriorArg$mu_log_sigma %||% log(0.7),
        PriorArg$sd_log_sigma %||% 1,
        log = TRUE
      )
  }

  if (length(spec$idx$npises) > 0L) {
    z <- theta[spec$idx$npises]
    mu <- PriorArg$mu_log_NPISes %||% 0
    sd <- PriorArg$sd_log_NPISes %||% 1
    if (length(mu) == 1L) {
      mu <- rep(mu, length(z))
    }
    if (length(sd) == 1L) {
      sd <- rep(sd, length(z))
    }
    if (length(mu) != length(z) || length(sd) != length(z)) {
      stop("mu_log_NPISes / sd_log_NPISes length mismatch.")
    }
    lp <- lp + sum(dnorm(z, mu, sd, log = TRUE))
  }

  lp
}


# Evaluate one theta --------------------------------------------------------

#' Compute log-likelihood, log-prior, and log-posterior for a given theta.
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
  after = as.Date("2015-08-31"),
  before = as.Date("2020-01-06"),
  Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
  LLHArg = list(),
  PriorArg = list(),
  InferArg = list()
) {
  Method <- match.arg(Method)
  n_virus <- length(BaseParm[["beta0"]])

  decoded <- Inference.DecodeTheta(
    theta = theta,
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )

  if (any(!is.finite(decoded$comp))) {
    return(list(llh = -Inf, lpr = -Inf, lpo = -Inf))
  }

  Parm <- BaseParm
  Parm[["comp"]] <- decoded$comp
  if (!is.null(decoded$Penal)) {
    Parm[["Penal"]] <- decoded$Penal
  }
  if (!is.null(decoded$NPISes)) {
    Parm[["NPISes"]] <- decoded$NPISes
  }

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

  llh <- Model.RunSim.LLH(
    Parm = Parm,
    Prep = Prep,
    Method = Method,
    LLHArg = LLHArg_cur,
    after = after,
    before = before
  )
  lpr <- MCMC.LogPrior.Theta(
    theta = theta,
    Method = Method,
    n_virus = n_virus,
    PriorArg = PriorArg,
    InferArg = InferArg
  )
  lpo <- if (is.finite(llh) && is.finite(lpr)) llh + lpr else -Inf

  list(llh = llh, lpr = lpr, lpo = lpo)
}


# Metropolis-Hastings Markov Chain Monte Carlo ----------------------------------------------

#' Adaptive Metropolis-Hastings sampler with Welford covariance adaptation.
#' Returns a matrix of MCMC samples on the transformed scale with attached attributes
#' (LLH, LogPrior, LogPost, Accepted, AcceptanceRate, Method, Spec, Adaptive, ...).
#' Use MCMC.DecodeChain / Inference.DecodeTheta to recover natural-scale parameters.
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
  before = as.Date("2020-01-01"),
  LLHArg = list(),
  PriorArg = list(),
  InferArg = list(),
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
      before = before,
      viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
    )
  }

  spec <- Inference.ParamSpec(
    Method = Method,
    n_virus = n_virus,
    InferArg = InferArg
  )
  theta_init <- MCMC.MakeInitial(
    Initial = Initial,
    Method = Method,
    n_virus = n_virus,
    BaseParm = BaseParm,
    InferArg = InferArg
  )

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
    PriorArg = PriorArg,
    InferArg = InferArg
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
    n_virus = n_virus,
    InferArg = InferArg
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
  scale_hist[1] <- 1

  cov_state <- MCMC.CovUpdate(MCMC.CovInit(spec$n_par), theta_init)
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
    proposal <- MCMC.Proposal.MVN(Theta = chain[i - 1, ], Sigma = Sigma_prop)
    eval_prop <- MCMC.EvaluateTheta(
      theta = proposal,
      BaseParm = BaseParm,
      Prep = Prep,
      Method = Method,
      LLHArg = LLHArg,
      PriorArg = PriorArg,
      InferArg = InferArg
    )

    if (
      is.finite(eval_prop$lpo) &&
        log(runif(1)) < eval_prop$lpo - lpo_trace[i - 1]
    ) {
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

    cov_state <- MCMC.CovUpdate(cov_state, chain[i, ])

    if (adapt_use && i >= adapt_start && i <= adapt_end) {
      acc_batch <- acc_batch + accepted[i]

      if (((i - adapt_start + 1L) %% adapt_every) == 0L) {
        batch_id <- (i - adapt_start + 1L) / adapt_every
        gamma_i <- 1 / ((batch_id + 10)^gamma_exp)
        log_lambda <- log_lambda +
          gamma_i * (acc_batch / adapt_every - target_accept)
        Sigma_prop <- exp(2 * log_lambda) *
          ((2.38^2 / spec$n_par) *
            MCMC.CovGet(cov_state, eps = eps_cov) +
            diag(eps_cov, spec$n_par))
        acc_batch <- 0L
      }
    }

    scale_hist[i] <- exp(log_lambda)
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
  attr(chain, "InferArg") <- InferArg
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

  chain
}


# Post-processing ----------------------------------------------------------------

#' Decode the MCMC chain from the transformed scale back to the natural scale.
#' @param chain        MCMC chain matrix (with Spec attribute, or provide n_virus)
#' @param include_eta  if TRUE, prepend the raw free-parameter columns to the output
MCMC.DecodeChain <- function(
  chain,
  Method = attr(chain, "Method"),
  n_virus = NULL,
  InferArg = attr(chain, "InferArg"),
  include_eta = FALSE
) {
  chain <- as.matrix(chain)
  spec <- attr(chain, "Spec")

  if (is.null(spec)) {
    if (is.null(n_virus)) {
      stop("Please provide n_virus.")
    }
    spec <- Inference.ParamSpec(
      Method = Method,
      n_virus = n_virus,
      InferArg = InferArg %||% list()
    )
  } else {
    Method <- spec$Method
    n_virus <- spec$n_virus
  }

  ia <- spec$Infer
  parts <- list()
  free_block <- chain[, spec$idx$comp, drop = FALSE]

  if (include_eta) {
    eta_block <- free_block
    colnames(eta_block) <- spec$comp_names
    parts <- c(parts, list(eta_block))
  }

  if (ia$penal_center_with_comp) {
    cp_full <- cbind(free_block, -rowSums(free_block))
    colnames(cp_full) <- c(paste0("comp_", seq_len(n_virus)), "Penal")
    parts <- c(parts, list(cp_full))
  } else {
    comp_full <- cbind(free_block, -rowSums(free_block))
    colnames(comp_full) <- paste0("comp_", seq_len(n_virus))
    parts <- c(parts, list(comp_full))

    if (length(spec$idx$penal) == 1L) {
      pen <- matrix(
        chain[, spec$idx$penal],
        ncol = 1,
        dimnames = list(NULL, "Penal")
      )
      parts <- c(parts, list(pen))
    }
  }

  idx_m <- spec$idx$method_aux
  aux <- switch(
    Method,
    BetaBinomial = matrix(
      plogis(chain[, idx_m[1L]]),
      ncol = 1,
      dimnames = list(NULL, "rho")
    ),
    Dirichlet = matrix(
      exp(chain[, idx_m[1L]]),
      ncol = 1,
      dimnames = list(NULL, "kappa")
    ),
    RatioAbs = cbind(
      sigma_ratio = exp(chain[, idx_m[1L]]),
      sigma_abs = exp(chain[, idx_m[2L]])
    ),
    LogDiff = matrix(
      exp(chain[, idx_m[1L]]),
      ncol = 1,
      dimnames = list(NULL, "sigma")
    )
  )
  parts <- c(parts, list(aux))

  if (length(spec$idx$npises) > 0L) {
    z <- chain[, spec$idx$npises, drop = FALSE]
    if (ia$NPISes_shared) {
      npm <- matrix(
        exp(z[, 1L]),
        ncol = 1,
        dimnames = list(NULL, "NPISes_shared")
      )
    } else {
      npm <- exp(z)
      colnames(npm) <- paste0("NPISes_", seq_len(ncol(npm)))
    }
    parts <- c(parts, list(npm))
  }

  do.call(cbind, parts)
}

# Parallelization ----------------------------------------------------------------

# MCMC.RunChains.Future <- function(
#   n_chain = 4L,
#   InitialList = NULL,
#   n_iterations,
#   step = 0.1,
#   TargetDat = NULL,
#   Prep = NULL,
#   Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
#   BaseParm = Parameter.Create(),
#   after = as.Date("2015-08-31"),
#   before = as.Date("2020-01-01"),
#   LLHArg = list(),
#   PriorArg = list(),
#   Adapt = list(),
#   workers = min(n_chain, future::availableCores()),
#   seed = 123,
#   source_r = NULL,
#   source_cpp = NULL,
#   wait = TRUE
# ) {
#   Method <- match.arg(Method)
#   n_virus <- length(BaseParm[["beta0"]])

#   if (is.null(Prep)) {
#     if (is.null(TargetDat)) {
#       stop("Provide either Prep or TargetDat.")
#     }
#     Prep <- Inference.Setup(
#       TargetDat = TargetDat,
#       BaseParm = BaseParm,
#       after = after,
#       before = before,
#       viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
#     )
#   }

#   if (is.null(InitialList)) {
#     InitialList <- MCMC.MakeInitialList(
#       n_chain = n_chain,
#       Initial = NULL,
#       Method = Method,
#       n_virus = n_virus,
#       jitter = 0.3,
#       seed = seed
#     )
#   }
#   stopifnot(length(InitialList) == n_chain)

#   old_plan <- plan()
#   plan(multisession, workers = workers)

#   # 仅同步模式下自动恢复 plan
#   if (wait) {
#     on.exit(plan(old_plan), add = TRUE)
#   }

#   fs <- lapply(seq_len(n_chain), function(k) {
#     future(
#       {
#         if (!is.null(source_r)) {
#           source(source_r)
#         }

#         if (!is.null(source_cpp)) {
#           # 避免并发 sourceCpp 冲突：每个 worker 用独立 cache 目录
#           cache_dir <- file.path(tempdir(), paste0("rcpp_cache_", Sys.getpid()))
#           dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
#           sourceCpp(source_cpp, cacheDir = cache_dir, rebuild = FALSE)
#         }

#         MCMC.MH.Adaptive(
#           Initial = InitialList[[k]],
#           n_iterations = n_iterations,
#           step = step,
#           Prep = Prep,
#           Method = Method,
#           BaseParm = BaseParm,
#           after = after,
#           before = before,
#           LLHArg = LLHArg,
#           PriorArg = PriorArg,
#           Adapt = Adapt,
#           verbose = FALSE
#         )
#       },
#       seed = seed + k
#     )
#   })
#   names(fs) <- paste0("chain_", seq_len(n_chain))

#   if (wait) {
#     chains <- lapply(fs, future::value)
#     names(chains) <- names(fs)
#     return(chains)
#   }

#   structure(
#     list(futures = fs, old_plan = old_plan),
#     class = "mcmc_future_job"
#   )
# }

MCMC.ToLong <- function(chain) {
  Dat <- do.call(
    rbind,
    lapply(seq_along(chain), function(i) {
      ch <- chain[[i]]
      mat <- as.matrix(ch)
      n <- nrow(mat)

      mcpar <- attr(ch, "mcpar")
      iter <- seq(from = mcpar[1], to = mcpar[2], by = mcpar[3])

      data.frame(
        Iteration = iter,
        Chain = rep(paste0("Chain ", i), n),
        as.data.frame(mat, check.names = FALSE),
        check.names = FALSE
      )
    })
  )

  pivot_longer(
    Dat,
    cols = -c(Iteration, Chain),
    names_to = "Parameter",
    values_to = "Value"
  )
}


#' Trace Plot
MCMC.TracePlot <- function(chain) {
  df <- MCMC.ToLong(chain)
  ggplot(df, aes(x = Iteration, y = Value, color = Chain)) +
    geom_line(alpha = 0.7, linewidth = 0.4) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 2) +
    labs(x = "Iteration", y = NULL, color = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}


#' Density Plot
MCMC.DensPlot <- function(chain) {
  df <- MCMC.ToLong(chain)
  ggplot(df, aes(x = Value, color = Chain, fill = Chain)) +
    geom_density(alpha = 0.25, linewidth = 0.6) +
    facet_wrap(~Parameter, scales = "free", ncol = 2) +
    labs(x = NULL, y = "Density", color = NULL, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")
}


#' Simulate epidemic curves
MCMC.Simulate <- function(
  ModelParm = list(),
  comp,
  n_virus = 4L,
  after = as.Date("2019-01-01"),
  ObsDat = NULL
) {
  virus_cols <- paste0("Virus_", 1:n_virus)
  OverrideParm <- do.call(Parameter.Create, c(list(comp = comp), ModelParm))
  SimData <- Model.RunSim.ode(
    Parm = OverrideParm,
    after = after
  )
  # Modify SimData
  SimPlotData <- copy(SimData$Data)
  SimPlotData <- SimPlotData[,
    lapply(.SD, sum, na.rm = TRUE),
    by = .(ISOweek),
    .SDcols = virus_cols
  ]
  SimPlotData[,
    (virus_cols) := lapply(.SD, function(x) x / 10000 * 100),
    .SDcols = virus_cols
  ]
  SimPlotData <- melt(
    SimPlotData,
    id.vars = "ISOweek",
    measure.vars = patterns("Virus_"),
    variable.name = "Virus",
    value.name = "PosRate"
  )
  SimPlotData[, Group := "Simulated"]

  # Modify ObsDat
  PlotData <- copy(ObsDat)
  PlotData[,
    `:=`(
      Virus_1 = fifelse(IV_Tested > 0, IAV / IV_Tested * 100, NA_real_),
      Virus_2 = fifelse(IV_Tested > 0, IBV / IV_Tested * 100, NA_real_),
      Virus_3 = fifelse(RSV_Tested > 0, RSV / RSV_Tested * 100, NA_real_),
      Virus_4 = fifelse(RV_Tested > 0, RV / RV_Tested * 100, NA_real_)
    )
  ]
  PlotData[,
    c("IV_Tested", "IAV", "IBV", "RSV_Tested", "RSV", "RV_Tested", "RV") := NULL
  ]
  PlotData <- melt(
    PlotData,
    id.vars = c("Location", "ISOweek", "Monday"),
    measure.vars = patterns("Virus_"),
    variable.name = "Virus",
    value.name = "PosRate"
  )
  PlotData[, Group := "Observed"]

  MergeData <- rbindlist(
    list(PlotData, SimPlotData),
    use.names = TRUE,
    fill = TRUE
  )

  ggplot(
    MergeData,
    aes(x = ISOweek, y = PosRate, color = Group, group = Group)
  ) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    facet_wrap(~Virus, ncol = 1, scales = "free_y") +
    labs(
      x = "Time",
      y = "Positive Rate",
      color = ""
    ) +
    theme_minimal(base_size = 13)
}


#' Post-process MCMC chains: burn-in, thinning, diagnostics, and summary.
MCMC.PostProcess <- function(
  dat,
  n_virus = 4L,
  burn_in = 5000L,
  thin = 10L,
  conf = 0.95,
  include_eta = FALSE,
  plot = TRUE,
  after = as.Date("2015-08-31"),
  ModelParm = list(),
  ObsDat = NULL,
  return_chain = TRUE
) {
  # Decode, burn-in, thinning per chain
  chain_list <- vector("list", length(dat))
  for (i in seq_along(dat)) {
    decoded <- MCMC.DecodeChain(dat[[i]], include_eta = include_eta)
    decoded <- as.matrix(decoded)

    mat <- decoded[, seq_len(n_virus), drop = FALSE]
    colnames(mat) <- paste0("Virus_", seq_len(n_virus))

    n_raw <- nrow(mat)
    keep_idx <- seq.int(from = burn_in + 1L, to = n_raw, by = thin)

    # Try preserve original iteration index if rownames are numeric
    iter_raw <- suppressWarnings(as.numeric(rownames(mat)))
    if (length(iter_raw) != n_raw || anyNA(iter_raw)) {
      iter_raw <- seq_len(n_raw)
    }
    start_iter <- iter_raw[keep_idx[1L]]

    post_mat <- mat[keep_idx, , drop = FALSE]
    chain_list[[i]] <- coda::mcmc(post_mat, start = start_iter, thin = thin)
  }
  Chain <- coda::mcmc.list(chain_list)
  ChainSim <- lapply(Chain, \(x) x[, 1:(n_virus - 1L), drop = FALSE])

  # Classical diagnostics
  Geltest <- if (length(ChainSim) >= 2L) {
    tryCatch(
      coda::gelman.diag(ChainSim, autoburnin = FALSE, multivariate = TRUE),
      error = function(e) e
    )
  }
  Heitest <- tryCatch(coda::heidel.diag(ChainSim), error = function(e) e)
  Geweke <- tryCatch(coda::geweke.diag(Chain), error = function(e) e)
  ESS <- tryCatch(coda::effectiveSize(ChainSim), error = function(e) e)

  # Posterior summary
  draws_mat <- as.matrix(Chain) # combined chains
  param_names <- colnames(draws_mat)

  Mean <- colMeans(draws_mat)
  Median <- apply(draws_mat, 2, stats::median)
  SD <- apply(draws_mat, 2, stats::sd)

  hpd_obj <- coda::HPDinterval(coda::as.mcmc(draws_mat), prob = conf)
  CI <- hpd_obj[, c("lower", "upper"), drop = FALSE]

  q_low <- (1 - conf) / 2
  q_high <- 1 - q_low
  ETI <- t(apply(
    draws_mat,
    2,
    stats::quantile,
    probs = c(q_low, q_high),
    names = FALSE
  ))
  colnames(ETI) <- c("ETI_lower", "ETI_upper")

  ess_vec <- rep(NA_real_, length(param_names))
  names(ess_vec) <- param_names
  if (!inherits(ESS, "error") && is.numeric(ESS)) {
    hit <- intersect(names(ESS), param_names)
    ess_vec[hit] <- ESS[hit]
  }

  MCSE_mean_classic <- SD / sqrt(ess_vec)

  Summary <- data.frame(
    Parameter = param_names,
    Mean = Mean[param_names],
    Median = Median[param_names],
    SD = SD[param_names],
    HPD_lower = CI[param_names, "lower"],
    HPD_upper = CI[param_names, "upper"],
    ETI_lower = ETI[param_names, "ETI_lower"],
    ETI_upper = ETI[param_names, "ETI_upper"],
    ESS_classic = ess_vec[param_names],
    MCSE_mean_classic = MCSE_mean_classic[param_names],
    check.names = FALSE,
    row.names = NULL
  )

  posteriori <- sprintf(
    "%.6f (%.6f, %.6f)",
    Summary$Median,
    Summary$HPD_lower,
    Summary$HPD_upper
  )
  names(posteriori) <- Summary$Parameter

  # Plots
  Traceplot <- NULL
  Densplot <- NULL
  Rankplot <- NULL
  ACFplot <- NULL
  if (isTRUE(plot)) {
    Traceplot <- MCMC.TracePlot(Chain)
    # print(Traceplot)
    Densplot <- MCMC.DensPlot(Chain)
    # print(Densplot)
    # Rankplot <- bayesplot::mcmc_rank_overlay(Chain)
    # print(Rankplot)
    # ACFplot <- bayesplot::mcmc_acf(Chain)
    # print(ACFplot)
    SimPlot <- MCMC.Simulate(
      ModelParm = ModelParm,
      comp = Median,
      n_virus = n_virus,
      after = after,
      ObsDat = ObsDat
    )
  }

  list(
    Chain = if (isTRUE(return_chain)) Chain else NULL,
    Traceplot = Traceplot,
    Densplot = Densplot,
    # Rankplot = Rankplot,
    # ACFplot = ACFplot,
    SimPlot = SimPlot,
    Geltest = Geltest,
    # Heitest = Heitest,
    Geweke = Geweke,
    ESS = ESS,
    Summary = Summary,
    CI = CI, # backward-compatible name (HPD CI)
    posteriori = posteriori
  )
}
