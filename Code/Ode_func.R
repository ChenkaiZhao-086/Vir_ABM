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
  dt <- Sim |> as.data.table()
  dt[, time := as.Date(time, origin = "1970-01-01")]
  inc_cols <- grep("^Inc_\\d+$", names(dt), value = TRUE)

  out <- dt[time > after, c("time", inc_cols), with = FALSE]

  if (!is.null(virus_names)) {
    stopifnot(length(virus_names) == length(inc_cols))
    data.table::setnames(out, inc_cols, virus_names)
  }
  return(out)
}

#' Run the ODE simulation with the given parameters and return the results.
#'
#' @param Parm: a list of parameters for the ODE model
#'
#' @return a list containing the raw simulation output and the processed incidence data
Model.RunSim.ode <- function(Parm) {
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
  SimResult_Inc <- Model.GetI(SimResult)

  return(list(
    DataRaw = SimResult,
    Data = SimResult_Inc
  ))
}
