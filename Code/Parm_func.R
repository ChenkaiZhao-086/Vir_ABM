#' Create Parmeters for model
#'
#' @return a list of parameters
Parameter.Create <- function(
  # Model base
  num_of_agent = 100000,
  dt = 1,
  years = 16,
  year_start = "1999-01-01",
  year_end = "2024-03-31",
  initial_seeds = 10,
  added_cases = 2, # 每个病毒随机每天新引入的I，这里是每dt2个，相当于一周4个
  penalty = 0.5, # 感染多病毒的惩罚系数

  # R0
  R0 = c(1.41, 1.07, 1.7, 1.88),

  # Seasonal force
  beta_seasonal = 0.2,
  phi = 340,
  beta_amplify = 1,

  # Duration of infectious
  gamma = c(1 / 6, 1 / 4, 1 / 7.4, 1 / 10.9),

  # Duration of immunity of each virus
  omega = c(1 / (365 * 2), 1 / 424.1, 1 / 358.9, 1 / 36.5),

  # Virus competition
  comp = c(1, 1, 1, 1),

  # NPI susceptibility
  NPISes = c(1, 1, 1, 1),

  # Base immunity
  base_immune = 0.4,

  # NPI
  NPI = TRUE,
  NPI_value = c(0.8),
  NPI_start = c("2020-03-10"), # not use this param any more
  NPI_end = c("2020-04-10") # not use this param any more
  # NPI_value = c(0), # range[0,1] 1 means totally no transmission(Full NPI)
  # decay_coef = c(0.01) # range[0,Inf] 0 means no decay
) {
  beta0 <- (R0 * gamma) / (1 + beta_seasonal)

  NPI_start_num <- as.numeric(as.Date(NPI_start))
  NPI_end_num <- as.numeric(as.Date(NPI_end))
  NPI_value <- as.numeric(NPI_value)

  Parmeters <-
    list(
      # Model base
      num_of_agent = num_of_agent,
      dt = dt,
      years = years,
      year_start = year_start,
      year_end = year_end,
      initial_seeds = initial_seeds,
      added_cases = added_cases,
      Penalty = penalty,

      # Seasonal force
      beta_seasonal = beta_seasonal,
      phi = phi,
      beta_amplify = beta_amplify,

      # Transmission rate
      beta0 = beta0,

      # Duration of infectious
      gamma = gamma,

      # Duration of immunity of each virus
      omega = omega,

      # Virus competition
      comp = comp,

      # NPI susceptibility
      NPISes = NPISes,

      # NPI
      NPI = NPI,
      NPI_start = NPI_start_num,
      NPI_end = NPI_end_num,
      NPI_value = NPI_value,
      # decay_coef = decay_coef,

      # Base immunity
      base_immune = base_immune
    )
  return(Parmeters)
}
