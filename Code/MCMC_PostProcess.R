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

  if (length(spec$idx$beta_seasonal) == 1L) {
    parts <- c(
      parts,
      list(matrix(
        plogis(chain[, spec$idx$beta_seasonal]),
        ncol = 1,
        dimnames = list(NULL, "beta_seasonal")
      ))
    )
  }

  if (length(spec$idx$phi) == 1L) {
    parts <- c(
      parts,
      list(matrix(
        Utility.Wrap365(chain[, spec$idx$phi]),
        ncol = 1,
        dimnames = list(NULL, "phi")
      ))
    )
  }

  if (length(spec$idx$npises) > 0L) {
    z <- chain[, spec$idx$npises, drop = FALSE]
    npm <- if (ia$NPISes_shared) {
      matrix(exp(z[, 1L]), ncol = 1, dimnames = list(NULL, "NPISes_shared"))
    } else {
      m <- exp(z)
      colnames(m) <- paste0("NPISes_", seq_len(ncol(m)))
      m
    }
    parts <- c(parts, list(npm))
  }

  do.call(cbind, parts)
}


Utility.CircMean365 <- function(phi, period = 365) {
  ang <- 2 * pi * (Utility.Wrap365(phi, period = period) - 1) / period
  mu <- atan2(mean(sin(ang), na.rm = TRUE), mean(cos(ang), na.rm = TRUE))
  if (mu < 0) {
    mu <- mu + 2 * pi
  }
  1 + period * mu / (2 * pi)
}

MCMC.CircSummary365 <- function(phi, conf = 0.95, period = 365) {
  center <- Utility.CircMean365(phi, period = period)
  phi_u <- center + Utility.CircDist365(phi, center = center, period = period)

  hpd <- HPDinterval(
    as.mcmc(matrix(phi_u, ncol = 1)),
    prob = conf
  )[1, ]

  q <- stats::quantile(
    phi_u,
    probs = c((1 - conf) / 2, 1 - (1 - conf) / 2),
    names = FALSE,
    na.rm = TRUE
  )

  c(
    Mean = Utility.Wrap365(center, period = period),
    Median = Utility.Wrap365(
      median(phi_u, na.rm = TRUE),
      period = period
    ),
    SD = sd(phi_u, na.rm = TRUE),
    HPD_lower = Utility.Wrap365(hpd["lower"], period = period),
    HPD_upper = Utility.Wrap365(hpd["upper"], period = period),
    ETI_lower = Utility.Wrap365(q[1], period = period),
    ETI_upper = Utility.Wrap365(q[2], period = period)
  )
}


MCMC.ToLong <- function(chain) {
  Dat <- do.call(
    rbind,
    lapply(seq_along(chain), function(i) {
      ch <- chain[[i]]
      mat <- as.matrix(ch)
      n <- nrow(mat)

      mcpar <- attr(ch, "mcpar")
      iter <- if (!is.null(mcpar)) {
        seq(from = mcpar[1], to = mcpar[2], by = mcpar[3])
      } else {
        seq_len(n)
      }
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
# MCMC.Simulate <- function(
#   ModelParm = list(),
#   comp,
#   n_virus = 4L,
#   after = as.Date("2019-01-01"),
#   ObsDat = NULL
# ) {
#   virus_cols <- paste0("Virus_", 1:n_virus)
#   OverrideParm <- do.call(Parameter.Create, c(list(comp = comp), ModelParm))
#   SimData <- Model.RunSim.ode(
#     Parm = OverrideParm,
#     after = after
#   )
#   # Modify SimData
#   SimPlotData <- copy(SimData$Data)
#   SimPlotData <- SimPlotData[,
#     lapply(.SD, sum, na.rm = TRUE),
#     by = .(ISOweek),
#     .SDcols = virus_cols
#   ]
#   SimPlotData[,
#     (virus_cols) := lapply(.SD, function(x) x / 10000 * 100),
#     .SDcols = virus_cols
#   ]
#   SimPlotData <- melt(
#     SimPlotData,
#     id.vars = "ISOweek",
#     measure.vars = patterns("Virus_"),
#     variable.name = "Virus",
#     value.name = "PosRate"
#   )
#   SimPlotData[, Group := "Simulated"]

#   # Modify ObsDat
#   PlotData <- copy(ObsDat)
#   PlotData[,
#     `:=`(
#       Virus_1 = fifelse(IV_Tested > 0, IAV / IV_Tested * 100, NA_real_),
#       Virus_2 = fifelse(IV_Tested > 0, IBV / IV_Tested * 100, NA_real_),
#       Virus_3 = fifelse(RSV_Tested > 0, RSV / RSV_Tested * 100, NA_real_),
#       Virus_4 = fifelse(RV_Tested > 0, RV / RV_Tested * 100, NA_real_)
#     )
#   ]
#   PlotData[,
#     c("IV_Tested", "IAV", "IBV", "RSV_Tested", "RSV", "RV_Tested", "RV") := NULL
#   ]
#   PlotData <- melt(
#     PlotData,
#     id.vars = c("Location", "ISOweek", "Monday"),
#     measure.vars = patterns("Virus_"),
#     variable.name = "Virus",
#     value.name = "PosRate"
#   )
#   PlotData[, Group := "Observed"]

#   MergeData <- rbindlist(
#     list(PlotData, SimPlotData),
#     use.names = TRUE,
#     fill = TRUE
#   )

#   MergeData[, Monday := as.Date(ISOweek::ISOweek2date(paste0(ISOweek, "-1")))]

#   ggplot(
#     MergeData,
#     aes(x = Monday, y = PosRate, color = Group, group = Group)
#   ) +
#     geom_line(linewidth = 0.9, na.rm = TRUE) +
#     facet_wrap(~Virus, ncol = 1, scales = "free_y") +
#     scale_x_date(date_labels = "%Y-%m", date_breaks = "4 months") +
#     labs(
#       x = "Time",
#       y = "Positive Rate",
#       color = ""
#     ) +
#     theme_minimal(base_size = 13) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     theme(legend.position = "none")
# }

MCMC.Simulate <- function(
  ModelParm = list(),
  comp = NULL,
  n_virus = 4L,
  after = as.Date("2019-01-01"),
  ObsDat = NULL
) {
  virus_cols <- paste0("Virus_", seq_len(n_virus))
  comp_nm <- paste0("comp_", seq_len(n_virus))

  if (
    !is.null(comp) && !is.null(names(comp)) && all(comp_nm %in% names(comp))
  ) {
    if (length(ModelParm) == 0L) {
      extra_nm <- setdiff(names(comp), comp_nm)
      if (length(extra_nm) > 0L) {
        ModelParm <- as.list(comp[extra_nm])
      }
    }
    comp <- as.numeric(comp[comp_nm])
  }

  if (is.null(comp)) {
    comp <- as.numeric(unlist(ModelParm[comp_nm], use.names = FALSE))
    ModelParm[comp_nm] <- NULL
  } else {
    comp <- as.numeric(comp)
    ModelParm[intersect(names(ModelParm), comp_nm)] <- NULL
  }

  ModelParm[c("phi_raw", "z_beta_seasonal", "kappa")] <- NULL

  OverrideParm <- do.call(
    Parameter.Create,
    c(list(comp = comp), ModelParm)
  )

  SimData <- Model.RunSim.ode(
    Parm = OverrideParm,
    after = after
  )

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
    measure.vars = virus_cols,
    variable.name = "Virus",
    value.name = "PosRate"
  )
  SimPlotData[, `:=`(
    Group = "Simulated",
    Location = NA_character_,
    Monday = as.Date(NA)
  )]

  PlotData <- copy(ObsDat)
  PlotData[,
    `:=`(
      Virus_1 = fifelse(IV_Tested > 0, IAV / IV_Tested * 100, NA_real_),
      Virus_2 = fifelse(IV_Tested > 0, IBV / IV_Tested * 100, NA_real_),
      Virus_3 = fifelse(RSV_Tested > 0, RSV / RSV_Tested * 100, NA_real_),
      Virus_4 = fifelse(RV_Tested > 0, RV / RV_Tested * 100, NA_real_)
    )
  ]
  drop_cols <- intersect(
    names(PlotData),
    c("IV_Tested", "IAV", "IBV", "RSV_Tested", "RSV", "RV_Tested", "RV")
  )
  if (length(drop_cols) > 0L) {
    PlotData[, (drop_cols) := NULL]
  }
  PlotData <- melt(
    PlotData,
    id.vars = intersect(c("Location", "ISOweek", "Monday"), names(PlotData)),
    measure.vars = virus_cols,
    variable.name = "Virus",
    value.name = "PosRate"
  )
  PlotData[, Group := "Observed"]

  MergeData <- rbindlist(
    list(PlotData, SimPlotData),
    use.names = TRUE,
    fill = TRUE
  )

  MergeData[
    is.na(Monday),
    Monday := as.Date(ISOweek::ISOweek2date(paste0(ISOweek, "-1")))
  ]
  MergeData[, Virus := factor(Virus, levels = virus_cols)]

  ggplot2::ggplot(
    MergeData,
    ggplot2::aes(x = Monday, y = PosRate, color = Group, group = Group)
  ) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::facet_wrap(~Virus, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_date(date_labels = "%Y-%m", date_breaks = "4 months") +
    ggplot2::labs(
      x = "Time",
      y = "Positive Rate",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}
#' Post-process MCMC chains: burn-in, thinning, diagnostics, and summary.
# MCMC.PostProcess <- function(
#   dat,
#   n_virus = 4L,
#   burn_in = 5000L,
#   thin = 10L,
#   conf = 0.95,
#   include_eta = FALSE,
#   plot = TRUE,
#   after = as.Date("2015-08-31"),
#   ModelParm = list(),
#   ObsDat = NULL,
#   return_chain = TRUE
# ) {
#   # Decode, burn-in, thinning per chain
#   chain_list <- vector("list", length(dat))
#   for (i in seq_along(dat)) {
#     decoded <- MCMC.DecodeChain(dat[[i]], include_eta = include_eta)
#     decoded <- as.matrix(decoded)

#     mat <- decoded[, seq_len(n_virus), drop = FALSE]
#     colnames(mat) <- paste0("Virus_", seq_len(n_virus))

#     n_raw <- nrow(mat)
#     keep_idx <- seq.int(from = burn_in + 1L, to = n_raw, by = thin)

#     # Try preserve original iteration index if rownames are numeric
#     iter_raw <- suppressWarnings(as.numeric(rownames(mat)))
#     if (length(iter_raw) != n_raw || anyNA(iter_raw)) {
#       iter_raw <- seq_len(n_raw)
#     }
#     start_iter <- iter_raw[keep_idx[1L]]

#     post_mat <- mat[keep_idx, , drop = FALSE]
#     chain_list[[i]] <- coda::mcmc(post_mat, start = start_iter, thin = thin)
#   }
#   Chain <- coda::mcmc.list(chain_list)
#   ChainSim <- lapply(Chain, \(x) x[, 1:(n_virus - 1L), drop = FALSE])

#   # Classical diagnostics
#   Geltest <- if (length(ChainSim) >= 2L) {
#     tryCatch(
#       coda::gelman.diag(ChainSim, autoburnin = FALSE, multivariate = TRUE),
#       error = function(e) e
#     )
#   }
#   Heitest <- tryCatch(coda::heidel.diag(ChainSim), error = function(e) e)
#   Geweke <- tryCatch(coda::geweke.diag(Chain), error = function(e) e)
#   ESS <- tryCatch(coda::effectiveSize(ChainSim), error = function(e) e)

#   # Posterior summary
#   draws_mat <- as.matrix(Chain) # combined chains
#   param_names <- colnames(draws_mat)

#   Mean <- colMeans(draws_mat)
#   Median <- apply(draws_mat, 2, stats::median)
#   SD <- apply(draws_mat, 2, stats::sd)

#   hpd_obj <- coda::HPDinterval(coda::as.mcmc(draws_mat), prob = conf)
#   CI <- hpd_obj[, c("lower", "upper"), drop = FALSE]

#   q_low <- (1 - conf) / 2
#   q_high <- 1 - q_low
#   ETI <- t(apply(
#     draws_mat,
#     2,
#     stats::quantile,
#     probs = c(q_low, q_high),
#     names = FALSE
#   ))
#   colnames(ETI) <- c("ETI_lower", "ETI_upper")

#   ess_vec <- rep(NA_real_, length(param_names))
#   names(ess_vec) <- param_names
#   if (!inherits(ESS, "error") && is.numeric(ESS)) {
#     hit <- intersect(names(ESS), param_names)
#     ess_vec[hit] <- ESS[hit]
#   }

#   MCSE_mean_classic <- SD / sqrt(ess_vec)

#   Summary <- data.frame(
#     Parameter = param_names,
#     Mean = Mean[param_names],
#     Median = Median[param_names],
#     SD = SD[param_names],
#     HPD_lower = CI[param_names, "lower"],
#     HPD_upper = CI[param_names, "upper"],
#     ETI_lower = ETI[param_names, "ETI_lower"],
#     ETI_upper = ETI[param_names, "ETI_upper"],
#     ESS_classic = ess_vec[param_names],
#     MCSE_mean_classic = MCSE_mean_classic[param_names],
#     check.names = FALSE,
#     row.names = NULL
#   )

#   posteriori <- sprintf(
#     "%.6f (%.6f, %.6f)",
#     Summary$Median,
#     Summary$HPD_lower,
#     Summary$HPD_upper
#   )
#   names(posteriori) <- Summary$Parameter

#   # Plots
#   Traceplot <- NULL
#   Densplot <- NULL
#   Rankplot <- NULL
#   ACFplot <- NULL
#   if (isTRUE(plot)) {
#     Traceplot <- MCMC.TracePlot(Chain)
#     # print(Traceplot)
#     Densplot <- MCMC.DensPlot(Chain)
#     # print(Densplot)
#     # Rankplot <- bayesplot::mcmc_rank_overlay(Chain)
#     # print(Rankplot)
#     # ACFplot <- bayesplot::mcmc_acf(Chain)
#     # print(ACFplot)
#     SimPlot <- MCMC.Simulate(
#       ModelParm = ModelParm,
#       comp = Median,
#       n_virus = n_virus,
#       after = after,
#       ObsDat = ObsDat
#     )
#   }

#   list(
#     Chain = if (isTRUE(return_chain)) Chain else NULL,
#     Traceplot = Traceplot,
#     Densplot = Densplot,
#     # Rankplot = Rankplot,
#     # ACFplot = ACFplot,
#     SimPlot = SimPlot,
#     Geltest = Geltest,
#     # Heitest = Heitest,
#     Geweke = Geweke,
#     ESS = ESS,
#     Summary = Summary,
#     CI = CI, # backward-compatible name (HPD CI)
#     posteriori = posteriori
#   )
# }

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
  decoded_chain_list <- vector("list", length(dat))
  diag_chain_list <- vector("list", length(dat))

  for (i in seq_along(dat)) {
    raw_mat <- as.matrix(dat[[i]])
    decoded_mat <- as.matrix(
      MCMC.DecodeChain(dat[[i]], include_eta = include_eta)
    )

    n_raw <- nrow(raw_mat)
    keep_idx <- seq.int(from = burn_in + 1L, to = n_raw, by = thin)

    mcpar <- attr(dat[[i]], "mcpar")
    iter_raw <- if (!is.null(mcpar)) {
      seq(from = mcpar[1], to = mcpar[2], by = mcpar[3])
    } else {
      seq_len(n_raw)
    }
    start_iter <- iter_raw[keep_idx[1L]]

    decoded_post <- decoded_mat[keep_idx, , drop = FALSE]
    raw_post <- raw_mat[keep_idx, , drop = FALSE]

    decoded_chain_list[[i]] <- mcmc(
      decoded_post,
      start = start_iter,
      thin = thin
    )
    diag_chain_list[[i]] <- mcmc(
      raw_post,
      start = start_iter,
      thin = thin
    )
  }

  Chain <- mcmc.list(decoded_chain_list)
  ChainDiag <- mcmc.list(diag_chain_list)

  # Classical diagnostics
  Geltest <- if (length(ChainDiag) >= 2L) {
    tryCatch(
      coda::gelman.diag(ChainDiag, autoburnin = FALSE, multivariate = TRUE),
      error = function(e) e
    )
  }
  Heitest <- tryCatch(coda::heidel.diag(ChainDiag), error = function(e) e)
  Geweke <- tryCatch(coda::geweke.diag(ChainDiag), error = function(e) e)
  ESS <- tryCatch(coda::effectiveSize(ChainDiag), error = function(e) e)

  # Posterior summary
  draws_mat <- as.matrix(Chain)
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
    ess_map <- stats::setNames(param_names, param_names)

    spec0 <- attr(dat[[1]], "Spec")
    if (!is.null(spec0)) {
      ia0 <- spec0$Infer
      raw_method <- if (length(spec0$idx$method_aux) > 0L) {
        spec0$aux_names[spec0$idx$method_aux]
      } else {
        character()
      }

      # if (
      #   spec0$Method == "BetaBinomial" &&
      #     "rho" %in% param_names &&
      #     length(raw_method) >= 1L
      # ) {
      #   ess_map["rho"] <- raw_method[1L]
      # }
      # if (
      #   spec0$Method == "Dirichlet" &&
      #     "kappa" %in% param_names &&
      #     length(raw_method) >= 1L
      # ) {
      #   ess_map["kappa"] <- raw_method[1L]
      # }
      # if (
      #   spec0$Method == "RatioAbs" &&
      #     "sigma_ratio" %in% param_names &&
      #     length(raw_method) >= 2L
      # ) {
      #   ess_map["sigma_ratio"] <- raw_method[1L]
      #   ess_map["sigma_abs"] <- raw_method[2L]
      # }
      # if (
      #   spec0$Method == "LogDiff" &&
      #     "sigma" %in% param_names &&
      #     length(raw_method) >= 1L
      # ) {
      #   ess_map["sigma"] <- raw_method[1L]
      # }

      raw_npises <- if (length(spec0$idx$npises) > 0L) {
        spec0$aux_names[spec0$idx$npises]
      } else {
        character()
      }
      if (length(raw_npises) > 0L) {
        if (isTRUE(ia0$NPISes_shared)) {
          if ("NPISes_shared" %in% param_names) {
            ess_map["NPISes_shared"] <- raw_npises[1L]
          }
        } else {
          for (j in seq_along(raw_npises)) {
            nm <- paste0("NPISes_", j)
            if (nm %in% param_names) {
              ess_map[nm] <- raw_npises[j]
            }
          }
        }
      }
    }

    if ("beta_seasonal" %in% param_names && "z_beta_seasonal" %in% names(ESS)) {
      ess_map["beta_seasonal"] <- "z_beta_seasonal"
    }
    if ("phi" %in% param_names && "phi_raw" %in% names(ESS)) {
      ess_map["phi"] <- "phi_raw"
    }

    hit <- names(ess_map)[ess_map %in% names(ESS)]
    ess_vec[hit] <- ESS[unname(ess_map[hit])]
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

  if ("phi" %in% Summary$Parameter) {
    circ <- MCMC.CircSummary365(draws_mat[, "phi"], conf = conf)
    idx_phi <- Summary$Parameter == "phi"
    Summary[idx_phi, names(circ)] <- as.list(circ)
    Summary$MCSE_mean_classic[idx_phi] <- NA_real_
  }

  posterior_median <- stats::setNames(Summary$Median, Summary$Parameter)

  posteriori <- sprintf(
    "%.6f (%.6f, %.6f)",
    Summary$Median,
    Summary$HPD_lower,
    Summary$HPD_upper
  )
  names(posteriori) <- Summary$Parameter

  if (isTRUE(plot)) {
    Traceplot <- MCMC.TracePlot(Chain)
    Densplot <- MCMC.DensPlot(Chain)
    # Rankplot <- bayesplot::mcmc_rank_overlay(Chain)
    # print(Rankplot)
    # ACFplot <- bayesplot::mcmc_acf(Chain)
    # print(ACFplot)
    sim_parm <- utils::modifyList(
      if (is.list(ModelParm)) ModelParm else as.list(ModelParm),
      as.list(posterior_median)
    )

    SimPlot <- MCMC.Simulate(
      ModelParm = sim_parm,
      n_virus = n_virus,
      after = after,
      ObsDat = ObsDat
    )
  }

  list(
    Chain = if (isTRUE(return_chain)) Chain else NULL,
    ChainDiag = if (isTRUE(return_chain)) ChainDiag else NULL,
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
