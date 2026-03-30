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
  epsilon = 1e-9,
  add_roll_poisson = TRUE,
  roll_n = 52L,
  roll_align = "center",
  poisson_weight = 1,
  lambda_floor = 0
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

  llh_bb <- sum(
    LLH.BetaBinomial.PMF(
      y = round(y_mat[valid]),
      n = N_mat[valid],
      p = p_sim_mat[valid],
      rho = rho,
      epsilon = epsilon
    ),
    na.rm = TRUE
  )

  if (
    !isTRUE(add_roll_poisson) ||
      !is.finite(poisson_weight) ||
      poisson_weight == 0
  ) {
    return(llh_bb)
  }

  llh_roll <- LLH.RollingMaxPoisson(
    y_mat = y_mat,
    sim_mat = p_sim_mat * 100000,
    roll_n = roll_n,
    roll_align = roll_align,
    lambda_floor = lambda_floor
  )

  llh_bb + poisson_weight * llh_roll
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

  RollMax <- function(x) {
    frollmax(as.numeric(x), n = roll_n, align = roll_align, na.rm = TRUE)
  }
  Sim_roll <- apply(Sim, 2, RollMax)
  Obs_roll <- apply(Obs, 2, RollMax)
  Lambda <- (Sim / Sim_roll) * Obs_roll
  if (is.finite(lambda_floor) && lambda_floor > 0) {
    Lambda <- pmax(Lambda, lambda_floor)
  }

  Obs_use <- round(Obs)
  ok <- is.finite(Obs_use) & is.finite(Lambda) & Obs_use >= 0 & Lambda >= 0
  if (!any(ok)) {
    return(0)
  }

  ll <- matrix(NA_real_, nrow(Obs), ncol(Obs))
  ll[ok] <- dpois(Obs_use[ok], Lambda[ok], log = TRUE)
  sum(ll, na.rm = TRUE)
}


#' Calculate the log-likelihood for the Dirichlet method.
#' @param Y matrix of observed proportions
#' @param Alpha matrix of concentration parameters
LLH.DirMult.PMF <- function(Y, Alpha) {
  Y <- as.matrix(Y)
  Alpha <- as.matrix(Alpha)

  n <- rowSums(Y)
  alpha0 <- rowSums(Alpha)

  lgamma(n + 1) -
    rowSums(lgamma(Y + 1)) +
    lgamma(alpha0) -
    lgamma(n + alpha0) +
    rowSums(lgamma(Y + Alpha) - lgamma(Alpha))
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
LLH.CondDirMult <- function(
  y_mat,
  sim_mat,
  kappa,
  N_vec = NULL,
  c = 1e-8,
  weight = NULL,
  add_roll_poisson = TRUE,
  roll_n = 52L,
  roll_align = "center",
  poisson_weight = 1,
  lambda_floor = 0
) {
  Obs <- as.matrix(y_mat)
  Sim <- as.matrix(sim_mat)
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    return(-Inf)
  }

  keep <- apply(is.finite(Obs) & Obs >= 0, 1, all) &
    apply(is.finite(Sim), 1, all)
  if (!any(keep)) {
    return(-Inf)
  }

  Obs_keep <- Obs[keep, , drop = FALSE]
  if (any(abs(Obs_keep - round(Obs_keep)) > 1e-8)) {
    stop(
      "Conditional Dirichlet-multinomial requires y_mat to be counts (integers), not proportions."
    )
  }

  Y <- round(Obs_keep)
  M <- pmax(Sim[keep, , drop = FALSE], 0)

  # Optional sanity check only
  if (!is.null(N_vec)) {
    N_use <- N_vec[keep]
    if (any(is.finite(N_use) & rowSums(Y) > N_use)) {
      warning(
        "Some row sums of y_mat exceed N_vec. ",
        "Conditional Dirichlet-multinomial assumes mutually exclusive categories; ",
        "please check for coinfections or duplicated counting."
      )
    }
  }

  # Mean composition implied by the simulator
  P_sim <- sweep(M + c, 1, rowSums(M + c), "/")

  # Total concentration = kappa
  Alpha <- pmax(kappa * P_sim, 1e-12)
  llh_dm <- sum(LLH.DirMult.PMF(Y, Alpha))

  if (
    !isTRUE(add_roll_poisson) ||
      !is.finite(poisson_weight) ||
      poisson_weight == 0
  ) {
    return(llh_dm)
  }

  llh_roll <- LLH.RollingMaxPoisson(
    y_mat = y_mat,
    sim_mat = sim_mat,
    roll_n = roll_n,
    roll_align = roll_align,
    lambda_floor = lambda_floor
  )

  llh_dm + poisson_weight * llh_roll
}
