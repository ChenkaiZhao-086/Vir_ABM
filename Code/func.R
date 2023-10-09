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
    omega_IFV = 1 / 150,
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

Model.Sim <- function(parameters) {
  si <- rbinom(n = , size = , prob = )
}
