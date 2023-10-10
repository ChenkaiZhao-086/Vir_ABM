File.CreateMainFolder <- function(path) {
  CurrentDate <- format(Sys.time(), "%Y%m%d")
  if (!file.exists(paste0(path, CurrentDate))) {
    dir.create(paste0(path, CurrentDate))
  }
  FolderPath <- paste0(path, CurrentDate, "/")
  cli_print("Create folder:", FolderPath)

  return(FolderPath)
}

File.CreateSubFolder <- function(path, folder) {
  if (!file.exists(paste0(path, folder))) {
    dir.create(paste0(path, folder))
  }
  return(paste0(path, folder, "/"))
}

Create.parameter <- function() {

}

#' @title Print message to command line interface
#'
#' @description Command line interface: print message
#'
#' @param ... (parts of the) message to print
#' @param WARNING boolean, to print the message in red
#'
#' @keywords internal
cli_print <- function(..., WARNING = FALSE, FORCED = FALSE) {
  # get function arguments
  function_arguments <- as.list(match.call(expand.dots = FALSE))$...

  # get function-call environment (to retrieve variable from that environment)
  pf <- parent.frame()

  # parse list => make character vector
  f_out <- " "
  for (i in 1:length(function_arguments)) {
    f_out <- cbind(f_out, eval(unlist(function_arguments[[i]]), envir = pf))
  }

  # add a space to each function arguments
  function_arguments <- paste(f_out, collapse = " ")

  # set text color: black (default) or red (warning)
  web_color_black <- "\033[0;30m"
  web_color_red <- "\033[0;31m"
  text_color <- ifelse(WARNING, web_color_red, web_color_black)

  # print time + arguments (without spaces)
  cli_out <- paste0(c(
    'echo "', text_color, "[", format(Sys.time(), "%H:%M:%S"), "]",
    function_arguments, web_color_black, '"'
  ), collapse = "")

  # print if function is called by master-node or first slave
  if (!exists("par_nodes_info") ||
    Sys.getpid() == par_nodes_info$pid_master ||
    FORCED) {
    system(cli_out)
  }


  # add to R warnings
  if (WARNING) {
    cli_warning <- paste0(c(
      text_color, "[", format(Sys.time(), "%H:%M:%S"), "]",
      function_arguments, web_color_black
    ), collapse = "")
    warning(cli_warning,
      call. = FALSE, immediate. = FALSE
    )
  }
}

#########################################
## PROGRESS BAR
##########################################
cli_progress <- function(i_current, i_total, time_stamp_loop) {
  # print if function is called by first node
  if (exists("par_nodes_info") &&
    Sys.getpid() == par_nodes_info$pid_slave1) {
    # calculate progress
    progress_scen <- floor(i_current / i_total * 100)
    progress_time <- round(difftime(Sys.time(), time_stamp_loop, units = "min"), digits = 1)

    # estimate remaining time (after 15%)
    time_label <- ""
    if (progress_scen > 15 & progress_scen < 99) {
      estim_time <- round(progress_time / progress_scen * (100 - progress_scen), digits = 1)
      if (estim_time < 1) {
        estim_time <- "<1"
      }
      time_label <- paste0("[", estim_time, " min remaining]")
    }

    cli_print("RUNNING...", i_current, "/", i_total, time_label, FORCED = TRUE)
  }
}

################################################################
## START CLUSTER WITH PARALLEL WORKERS
################################################################
Parallel.Regist <- function() {
  cli_print("START PARALLEL WORKERS")

  ## SETUP PARALLEL NODES
  # note: they will be removed after 600 seconds inactivity
  num_proc <- detectCores()
  par_cluster <- makeCluster(num_proc, cores = num_proc, timeout = 600)
  registerDoParallel(par_cluster)

  # store the process id (pid) of the first slave
  pid_slave1 <- clusterEvalQ(par_cluster, {
    Sys.getpid()
  })[[1]]

  # CREATE GLOBAL VARIABLE
  par_nodes_info <<- list(
    par_cluster = par_cluster,
    pid_master = Sys.getpid(),
    pid_slave1 = pid_slave1
  )
}

################################################################
## STOP CLUSTER WITH PARALLEL WORKERS
################################################################
Parallel.Stop <- function() {
  ## CLOSE NODES AND NODE INFO
  if (exists("par_nodes_info")) {
    cli_print("STOP PARALLEL WORKERS")

    stopCluster(par_nodes_info$par_cluster)
    rm(par_nodes_info, envir = .GlobalEnv) # REMOVE GLOBAL VARIABLE
  }
}

################################################################
## CHECK IF CLUSTER EXISTS AND START ONE IF NOT
################################################################
Parallel.Check <- function() {
  if (!exists("par_nodes_info")) {
    start_parallel_workers()
  } else if (!any(grepl(par_nodes_info$pid_slave1, system("ps -A", intern = TRUE)))) {
    start_parallel_workers()
  }
}

################################################################
## RESET CLUSTER
################################################################
# reset parallel workers
Parallel.Reset <- function() {
  stop_parallel_workers()
  gc()
  start_parallel_workers()
}


Parameter.Create <- function() {
  parameters <- list(
    # Model base
    num_of_agent = 100000,
    dt = 1 / 7, # Length of each time step in weeks [1/7]
    years = 10,
    initial_seeds = 10,
    daily_new_infection = 10, # 每个病毒随机每天新引入的I

    # Transmission rate
    beta_IFV = 1,
    beta_IFVA = 1,
    beta_IFVB = 1,
    beta_RSV = 1,
    beta_HPIV = 1,
    beta_HMPV = 1,
    beta_HCoV = 1,
    beta_HRV = 1,
    beta_HAdV = 1,

    # Seasonal force
    beta_seasonal = 1,

    # Duration of infectious
    gamma = 1 / 7,

    # Duration of immunity of each virus
    omega_IFV = 1 / 365,
    omega_IFVA = 1 / 150,
    omega_IFVB = 1 / 150,
    omega_RSV = 1 / 150,
    omega_HPIV = 1 / 150,
    omega_HMPV = 1 / 150,
    omega_HCoV = 1 / 150,
    omega_HRV = 1 / 150,
    omega_HAdV = 1 / 150,

    # Virus competition
    comp_IFV = 1 / 150,
    comp_IFVA = 1 / 150,
    comp_IFVB = 1 / 150,
    comp_RSV = 1 / 150,
    comp_HPIV = 1 / 150,
    comp_HMPV = 1 / 150,
    comp_HCoV = 1 / 150,
    comp_HRV = 1 / 150,
    comp_HAdV = 1 / 150
  )

  return(parameters)
}

Model.Sim <- function(Parm) {
  # Time steps to run simulation
  steps <- (52 * Parm$years) / Parm$dt
  initial_seeds <- Parm$initial_seeds
  N <- Parm$num_of_agent
  dt <- Parm$dt

  # Create storage for all agent
  agent_status <- data.frame(
    IFV = rep("S", Parm$num_of_agent),
    IFVA = rep("S", Parm$num_of_agent),
    IFVB = rep("S", Parm$num_of_agent),
    RSV = rep("S", Parm$num_of_agent),
    HPIV = rep("S", Parm$num_of_agent),
    HMPV = rep("S", Parm$num_of_agent),
    HCoV = rep("S", Parm$num_of_agent),
    HRV = rep("S", Parm$num_of_agent),
    HAdV = rep("S", Parm$num_of_agent)
  )

  # Seed the simulation with each virus
  # agent_status$IFV[1:initial_seeds]  <- "I"
  # agent_status$IFVA[1:initial_seeds] <- "I"
  # agent_status$IFVB[1:initial_seeds] <- "I"
  # agent_status$RSV[1:initial_seeds]  <- "I"
  # agent_status$HPIV[1:initial_seeds] <- "I"
  # agent_status$HMPV[1:initial_seeds] <- "I"
  # agent_status$HCoV[1:initial_seeds] <- "I"
  # agent_status$HRV[1:initial_seeds]  <- "I"
  # agent_status$HAdV[1:initial_seeds] <- "I"
  agent_status[1:initial_seeds, c("IFV", "IFVA", "IFVB", "RSV", "HPIV", "HMPV", "HCoV", "HRV", "HAdV")] <- "I"

  # Create storage for simulation results
  results <- data.frame(
    time = (0:steps) * dt,
    pop_IFV = 0,
    pop_IFVA = 0,
    pop_IFVB = 0,
    pop_RSV = 0,
    pop_HPIV = 0,
    pop_HMPV = 0,
    pop_HCoV = 0,
    pop_HRV = 0,
    pop_HAdV = 0,
    agent1 = 0
  )

  # Store the initial conditions
  results$pop_IFV[1] <- sum(agent_status[, "IFV"] == "I", na.rm = TRUE)
  results$pop_IFVA[1] <- sum(agent_status[, "IFVA"] == "I", na.rm = TRUE)
  results$pop_IFVB[1] <- sum(agent_status[, "IFVB"] == "I", na.rm = TRUE)
  results$pop_RSV[1] <- sum(agent_status[, "RSV"] == "I", na.rm = TRUE)
  results$pop_HPIV[1] <- sum(agent_status[, "HPIV"] == "I", na.rm = TRUE)
  results$pop_HMPV[1] <- sum(agent_status[, "HMPV"] == "I", na.rm = TRUE)
  results$pop_HCoV[1] <- sum(agent_status[, "HCoV"] == "I", na.rm = TRUE)
  results$pop_HRV[1] <- sum(agent_status[, "HRV"] == "I", na.rm = TRUE)
  results$pop_HAdV[1] <- sum(agent_status[, "HAdV"] == "I", na.rm = TRUE)
  # results$agent1[1] <- agent_status$IFV[1]

  bar <- txtProgressBar(min = 1, max = steps, style = 3)

  # Run the simulation
  for (i in 1:steps) {
    # Save population status for this states
    results$pop_IFV[1 + i] <- sum(agent_status[, "IFV"] == "I", na.rm = TRUE)
    results$pop_IFVA[1 + i] <- sum(agent_status[, "IFVA"] == "I", na.rm = TRUE)
    results$pop_IFVB[1 + i] <- sum(agent_status[, "IFVB"] == "I", na.rm = TRUE)
    results$pop_RSV[1 + i] <- sum(agent_status[, "RSV"] == "I", na.rm = TRUE)
    results$pop_HPIV[1 + i] <- sum(agent_status[, "HPIV"] == "I", na.rm = TRUE)
    results$pop_HMPV[1 + i] <- sum(agent_status[, "HMPV"] == "I", na.rm = TRUE)
    results$pop_HCoV[1 + i] <- sum(agent_status[, "HCoV"] == "I", na.rm = TRUE)
    results$pop_HRV[1 + i] <- sum(agent_status[, "HRV"] == "I", na.rm = TRUE)
    results$pop_HAdV[1 + i] <- sum(agent_status[, "HAdV"] == "I", na.rm = TRUE)

    # Calculate the total infection for each virus
    I_IFV <- sum(agent_status[, "IFV"] == "I", na.rm = TRUE)
    I_IFVA <- sum(agent_status[, "IFVA"] == "I", na.rm = TRUE)
    I_IFVB <- sum(agent_status[, "IFVB"] == "I", na.rm = TRUE)
    I_RSV <- sum(agent_status[, "RSV"] == "I", na.rm = TRUE)
    I_HPIV <- sum(agent_status[, "HPIV"] == "I", na.rm = TRUE)
    I_HMPV <- sum(agent_status[, "HMPV"] == "I", na.rm = TRUE)
    I_HCoV <- sum(agent_status[, "HCoV"] == "I", na.rm = TRUE)
    I_HRV <- sum(agent_status[, "HRV"] == "I", na.rm = TRUE)
    I_HAdV <- sum(agent_status[, "HAdV"] == "I", na.rm = TRUE)

    # Calculate force of infection of each virus
    lambda_IFV <- beta * delta * I_IFV / N
    lambda_IFVA <- beta * delta * I_IFVA / N
    lambda_IFVB <- beta * delta * I_IFVB / N
    lambda_RSV <- beta * delta * I_RSV / N
    lambda_HPIV <- beta * delta * I_HPIV / N
    lambda_HMPV <- beta * delta * I_HMPV / N
    lambda_HCoV <- beta * delta * I_HCoV / N
    lambda_HRV <- beta * delta * I_HRV / N
    lambda_HAdV <- beta * delta * I_HAdV / N

    # Get gamma
    gamma <- Parm$gamma

    # Get omega   # 如果omega要变成一个分布，可以考虑在这里进行抽样
    omega_IFV <- Parm$omega_IFV
    omega_IFVA <- Parm$omega_IFVA
    omega_IFVB <- Parm$omega_IFVB
    omega_RSV <- Parm$omega_RSV
    omega_HPIV <- Parm$omega_HPIV
    omega_HMPV <- Parm$omega_HMPV
    omega_HCoV <- Parm$omega_HCoV
    omega_HRV <- Parm$omega_HRV
    omega_HAdV <- Parm$omega_HAdV

    # p_gamma <- 1 - exp(-gamma * dt)
    # p_omega_IFV <- 1 - exp(-omega_IFV * dt)

    # Select number of event
    ## IFV
    IFV_S2I <- rbinom(n = 1, size = sum(agent_status[, "IFV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_IFV * dt))
    IFV_I2R <- rbinom(n = 1, size = sum(agent_status[, "IFV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    IFV_R2S <- rbinom(n = 1, size = sum(agent_status[, "IFV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_IFV * dt))
    ## IFVA
    IFVA_S2I <- rbinom(n = 1, size = sum(agent_status[, "IFVA"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_IFVA * dt))
    IFVA_I2R <- rbinom(n = 1, size = sum(agent_status[, "IFVA"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    IFVA_R2S <- rbinom(n = 1, size = sum(agent_status[, "IFVA"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_IFVA * dt))
    ## IFVB
    IFVB_S2I <- rbinom(n = 1, size = sum(agent_status[, "IFVB"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_IFVB * dt))
    IFVB_I2R <- rbinom(n = 1, size = sum(agent_status[, "IFVB"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    IFVB_R2S <- rbinom(n = 1, size = sum(agent_status[, "IFVB"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_IFVB * dt))
    ## RSV
    RSV_S2I <- rbinom(n = 1, size = sum(agent_status[, "RSV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_RSV * dt))
    RSV_I2R <- rbinom(n = 1, size = sum(agent_status[, "RSV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    RSV_R2S <- rbinom(n = 1, size = sum(agent_status[, "RSV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_RSV * dt))
    ## HPIV
    HPIV_S2I <- rbinom(n = 1, size = sum(agent_status[, "HPIV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_HPIV * dt))
    HPIV_I2R <- rbinom(n = 1, size = sum(agent_status[, "HPIV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    HPIV_R2S <- rbinom(n = 1, size = sum(agent_status[, "HPIV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_HPIV * dt))
    ## HMPV
    HMPV_S2I <- rbinom(n = 1, size = sum(agent_status[, "HMPV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_HMPV * dt))
    HMPV_I2R <- rbinom(n = 1, size = sum(agent_status[, "HMPV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    HMPV_R2S <- rbinom(n = 1, size = sum(agent_status[, "HMPV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_HMPV * dt))
    ## HCoV
    HCoV_S2I <- rbinom(n = 1, size = sum(agent_status[, "HCoV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_HCoV * dt))
    HCoV_I2R <- rbinom(n = 1, size = sum(agent_status[, "HCoV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    HCoV_R2S <- rbinom(n = 1, size = sum(agent_status[, "HCoV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_HCoV * dt))
    ## HRV
    HRV_S2I <- rbinom(n = 1, size = sum(agent_status[, "HRV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_HRV * dt))
    HRV_I2R <- rbinom(n = 1, size = sum(agent_status[, "HRV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    HRV_R2S <- rbinom(n = 1, size = sum(agent_status[, "HRV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_HRV * dt))
    ## HAdV
    HAdV_S2I <- rbinom(n = 1, size = sum(agent_status[, "HAdV"] == "S", na.rm = TRUE), prob = 1 - exp(-lambda_HAdV * dt))
    HAdV_I2R <- rbinom(n = 1, size = sum(agent_status[, "HAdV"] == "I", na.rm = TRUE), prob = 1 - exp(-gamma * dt))
    HAdV_R2S <- rbinom(n = 1, size = sum(agent_status[, "HAdV"] == "R", na.rm = TRUE), prob = 1 - exp(-omega_HAdV * dt))

    # Choose agent
    ## IFV
    agent_status[sample(which(agent_status[, "IFV"] == "S"), IFV_S2I), "IFV"] <- "I"
    agent_status[sample(which(agent_status[, "IFV"] == "I"), IFV_I2R), "IFV"] <- "R"
    agent_status[sample(which(agent_status[, "IFV"] == "R"), IFV_R2S), "IFV"] <- "S"
    ## IFVA
    agent_status[sample(which(agent_status[, "IFVA"] == "S"), IFVA_S2I), "IFVA"] <- "I"
    agent_status[sample(which(agent_status[, "IFVA"] == "I"), IFVA_I2R), "IFVA"] <- "R"
    agent_status[sample(which(agent_status[, "IFVA"] == "R"), IFVA_R2S), "IFVA"] <- "S"
    ## IFVB
    agent_status[sample(which(agent_status[, "IFVB"] == "S"), IFVB_S2I), "IFVB"] <- "I"
    agent_status[sample(which(agent_status[, "IFVB"] == "I"), IFVB_I2R), "IFVB"] <- "R"
    agent_status[sample(which(agent_status[, "IFVB"] == "R"), IFVB_R2S), "IFVB"] <- "S"
    ## RSV
    agent_status[sample(which(agent_status[, "RSV"] == "S"), RSV_S2I), "RSV"] <- "I"
    agent_status[sample(which(agent_status[, "RSV"] == "I"), RSV_I2R), "RSV"] <- "R"
    agent_status[sample(which(agent_status[, "RSV"] == "R"), RSV_R2S), "RSV"] <- "S"
    ## HPIV
    agent_status[sample(which(agent_status[, "HPIV"] == "S"), HPIV_S2I), "HPIV"] <- "I"
    agent_status[sample(which(agent_status[, "HPIV"] == "I"), HPIV_I2R), "HPIV"] <- "R"
    agent_status[sample(which(agent_status[, "HPIV"] == "R"), HPIV_R2S), "HPIV"] <- "S"
    ## HMPV
    agent_status[sample(which(agent_status[, "HMPV"] == "S"), HMPV_S2I), "HMPV"] <- "I"
    agent_status[sample(which(agent_status[, "HMPV"] == "I"), HMPV_I2R), "HMPV"] <- "R"
    agent_status[sample(which(agent_status[, "HMPV"] == "R"), HMPV_R2S), "HMPV"] <- "S"
    ## HCoV
    agent_status[sample(which(agent_status[, "HCoV"] == "S"), HCoV_S2I), "HCoV"] <- "I"
    agent_status[sample(which(agent_status[, "HCoV"] == "I"), HCoV_I2R), "HCoV"] <- "R"
    agent_status[sample(which(agent_status[, "HCoV"] == "R"), HCoV_R2S), "HCoV"] <- "S"
    ## HRV
    agent_status[sample(which(agent_status[, "HRV"] == "S"), HRV_S2I), "HRV"] <- "I"
    agent_status[sample(which(agent_status[, "HRV"] == "I"), HRV_I2R), "HRV"] <- "R"
    agent_status[sample(which(agent_status[, "HRV"] == "R"), HRV_R2S), "HRV"] <- "S"
    ## HAdV
    agent_status[sample(which(agent_status[, "HAdV"] == "S"), HAdV_S2I), "HAdV"] <- "I"
    agent_status[sample(which(agent_status[, "HAdV"] == "I"), HAdV_I2R), "HAdV"] <- "R"
    agent_status[sample(which(agent_status[, "HAdV"] == "R"), HAdV_R2S), "HAdV"] <- "S"


    setTxtProgressBar(bar, i)
    if (i == steps) {
      close(bar)
    }
  }

  return(results)
}
