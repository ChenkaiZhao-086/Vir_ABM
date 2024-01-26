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

Create.Parmeter <- function() {

}

#' @title Print message to command line interface
#'
#' @description Command line interface: print message
#'
#' @Parm ... (parts of the) message to print
#' @Parm WARNING boolean, to print the message in red
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


#' @title Create Parmeters for model
#'
#' @description Create Parmeters for model
Parameter.Create <- function(
    # Model base
    num_of_agent = 100000,
    dt = 1 / 7 * 2, # Length of each time step in weeks [1/7]
    years = 16,
    initial_seeds = 10,
    added_cases = 2, # 每个病毒随机每天新引入的I，这里是每dt2个，相当于一周14个
    Penal = 0.5, # 感染多病毒的惩罚系数
    # R0
    R0_IFVA = 1.25,
    R0_IFVB = 1.07,
    R0_RSV = 4.5,
    R0_PIV = 2.5,
    R0_MPV = 1.2,
    R0_sCoV = 4.18,
    R0_RV = 1.88,
    R0_AdV = 5.1,
    # Transmission rate
    beta_IFVA = 1,
    beta_IFVB = 1,
    beta_RSV = 1,
    beta_PIV = 1,
    beta_MPV = 1,
    beta_sCoV = 1,
    beta_RV = 1,
    beta_AdV = 1,
    # Seasonal force
    beta_seasonal = 0.2,
    phi = 12,
    # Duration of infectious
    # gamma = 1 / 7, #
    gamma_IFVA = 1 / 6 * 2, # 1 / 7, #
    gamma_IFVB = 1 / 4 * 2, # 1 / 7, #
    gamma_RSV = 1 / 6.4 * 2, # 1 / 7, #
    gamma_PIV = 1 / 9.8 * 2, # 1 / 7, #
    gamma_MPV = 1 / 5 * 2, # 1 / 7, #
    gamma_sCoV = 1 / 15.2 * 2, # 1 / 7, #
    gamma_RV = 1 / 10.9 * 2, # 1 / 7, #
    gamma_AdV = 1 / 7 * 2,
    # Duration of immunity of each virus
    omega_IFVA = 1 / (365 * 2) * 2, # 150, #
    omega_IFVB = 1 / 424.1 * 2, # 150, #
    omega_RSV = 1 / 358.9 * 2, # 150, #
    omega_PIV = 1 / 303 * 2, # 150, #
    omega_MPV = 1 / 657 * 2, # 150, #
    omega_sCoV = 1 / 259 * 2, # 150, #
    omega_RV = 1 / 36.5 * 2, # 150, #
    omega_AdV = 1 / 365 * 2, # 150, #
    # Virus competition
    comp_IFVA = 1,
    comp_IFVB = 1,
    comp_RSV = 1,
    comp_PIV = 1,
    comp_MPV = 1,
    comp_sCoV = 1,
    comp_RV = 1,
    comp_AdV = 1,
    # NPI
    NPI_start = 624, # 第12年, 428第八年的三月
    NPI_end = 676, # 第13年
    NPI_value = 0.8, # range[0,1] 1 means totally no transmission(Full NPI), 0 means no NPI(free transimiison)
    decay_coef = 0.01, # range[0,Inf] 0 means no decay
    # Base immunity
    base_immune_IFVA = 40000,
    base_immune_IFVB = 40000,
    base_immune_RSV = 40000,
    base_immune_PIV = 40000,
    base_immune_MPV = 40000,
    base_immune_sCoV = 40000,
    base_immune_RV = 40000,
    base_immune_AdV = 40000) {
  beta_IFVA <- (R0_IFVA * gamma_IFVA) / (1 + beta_seasonal) #  R0_IFVA * gamma_IFVA #
  beta_IFVB <- (R0_IFVB * gamma_IFVB) / (1 + beta_seasonal) # R0_IFVB * gamma_IFVB # 
  beta_RSV <- (R0_RSV * gamma_RSV) / (1 + beta_seasonal) # R0_RSV * gamma_RSV # 
  beta_PIV <- (R0_PIV * gamma_PIV) / (1 + beta_seasonal) # R0_PIV * gamma_PIV # 
  beta_MPV <- (R0_MPV * gamma_MPV) / (1 + beta_seasonal) # R0_MPV * gamma_MPV # 
  beta_sCoV <- (R0_sCoV * gamma_sCoV) / (1 + beta_seasonal) # R0_sCoV * gamma_sCoV # 
  beta_RV <- (R0_RV * gamma_RV) / (1 + beta_seasonal) # R0_RV * gamma_RV # 
  beta_AdV <- (R0_AdV * gamma_AdV) / (1 + beta_seasonal) # R0_AdV * gamma_AdV # 
  Parmeters <-
    list(
      # Model base
      num_of_agent = num_of_agent,
      dt = dt,
      years = years,
      initial_seeds = initial_seeds,
      added_cases = added_cases,
      Penal = Penal,

      # Seasonal force
      beta_seasonal = beta_seasonal,
      phi = phi,

      # Duration of infectious
      gamma = gamma,

      # Transmission rate
      beta_IFVA = beta_IFVA,
      beta_IFVB = beta_IFVB,
      beta_RSV = beta_RSV,
      beta_PIV = beta_PIV,
      beta_MPV = beta_MPV,
      beta_sCoV = beta_sCoV,
      beta_RV = beta_RV,
      beta_AdV = beta_AdV,

      # Duration of infectious
      gamma_IFVA = gamma_IFVA,
      gamma_IFVB = gamma_IFVB,
      gamma_RSV = gamma_RSV,
      gamma_PIV = gamma_PIV,
      gamma_MPV = gamma_MPV,
      gamma_sCoV = gamma_sCoV,
      gamma_RV = gamma_RV,
      gamma_AdV = gamma_AdV,

      # Duration of immunity of each virus
      omega_IFVA = omega_IFVA,
      omega_IFVB = omega_IFVB,
      omega_RSV = omega_RSV,
      omega_PIV = omega_PIV,
      omega_MPV = omega_MPV,
      omega_sCoV = omega_sCoV,
      omega_RV = omega_RV,
      omega_AdV = omega_AdV,

      # Virus competition
      comp_IFVA = comp_IFVA,
      comp_IFVB = comp_IFVB,
      comp_RSV = comp_RSV,
      comp_PIV = comp_PIV,
      comp_MPV = comp_MPV,
      comp_sCoV = comp_sCoV,
      comp_RV = comp_RV,
      comp_AdV = comp_AdV,

      # NPI
      NPI_start = NPI_start,
      NPI_end = NPI_end,
      NPI_value = NPI_value,
      decay_coef = decay_coef,

      # Base immunity
      base_immune_IFVA = base_immune_IFVA,
      base_immune_IFVB = base_immune_IFVB,
      base_immune_RSV = base_immune_RSV,
      base_immune_PIV = base_immune_PIV,
      base_immune_MPV = base_immune_MPV,
      base_immune_sCoV = base_immune_sCoV,
      base_immune_RV = base_immune_RV,
      base_immune_AdV = base_immune_AdV
    )

  return(Parmeters)
}

#' @title Run simulation of model and plot
#'
#' @description Run simulation of model and plot, the ModelSimCpp function is written in C++
#' @Parm Parm Parmeters of model
#' @Parm ncores number of cores to run simulation
#' @Parm NPI whether to run NPI
#' @Parm BaseImmu whether to run base immunity
#' @Parm seeds random seeds
#' @return Return a list of simulation result and plots
Model.RunSim <- function(Parm, ncores = 6, NPI = FALSE, BaseImmu = FALSE, seeds = NULL) {
  if (!is.null(seeds)) set.seed(seeds)
  SimResult <- ModelSimCpp(Parm = Parm, ncores = ncores, NPI = NPI, BaseImmu = BaseImmu)
  colnames(SimResult) <- c(
    "time", "IFVA", "IFVB", "RSV", "PIV", "MPV", "sCoV", "RV", "AdV"
    # , "S_IFVA", "S_IFVB", "S_RSV", "S_PIV", "S_MPV", "S_sCoV", "S_RV", "S_AdV"
  )

  fig1 <- SimResult %>%
    as.data.frame() %>%
    dplyr::select(1:9) %>%
    pivot_longer(cols = !time, names_to = "virus", values_to = "cases") %>%
    mutate(virus = fct_relevel(virus, c("IFVA", "IFVB", "RSV", "PIV", "MPV", "sCoV", "RV", "AdV"))) %>%
    filter(time > 53*4) %>%
    ggplot(., aes(x = time, y = cases)) +
    geom_line() +
    scale_x_continuous(breaks = seq(1, length(SimResult[, 1]), by = 52)) +
    # scale_y_continuous(limits = c(0, 10000)) +
    theme_bw() +
    facet_wrap(vars(virus), nrow = 3, scales = "free_y")

  # fig3 <- SimResult %>%
  #   as.data.frame() %>%
  #   dplyr::select(1, 10:17) %>%
  #   pivot_longer(cols = !time, names_to = "virus", values_to = "cases") %>%
  #   mutate(virus = fct_relevel(virus, c("S_IFVA", "S_IFVB", "S_RSV", "S_PIV", "S_MPV", "S_sCoV", "S_RV", "S_AdV"))) %>%
  #   # filter(time > 100) %>%
  #   ggplot(., aes(x = time, y = cases)) +
  #   geom_line() +
  #   scale_x_continuous(breaks = seq(1, length(SimResult[, 1]), by = 52)) +
  #   # scale_y_continuous(limits = c(0, 10000)) +
  #   theme_bw() +
  #   facet_wrap(vars(virus), nrow = 3, scales = "free_y")

  fig2 <- SimResult %>%
    as.data.frame() %>%
    dplyr::select(1:9) %>%
    pivot_longer(cols = !time, names_to = "virus", values_to = "cases") %>%
    mutate(virus = fct_relevel(virus, c("IFVA", "IFVB", "RSV", "PIV", "MPV", "sCoV", "RV", "AdV"))) %>%
    filter(time > 53*4) %>%
    ggplot(., aes(x = time, y = cases, colour = virus)) +
    geom_line(alpha = 0.7) +
    scale_x_continuous(breaks = seq(1, length(SimResult[, 1]), by = 52)) +
    # scale_y_continuous(limits = c(0, 10000)) +
    theme_bw()

  # fig4 <- SimResult %>%
  #   as.data.frame() %>%
  #   dplyr::select(1, 10:17) %>%
  #   pivot_longer(cols = !time, names_to = "virus", values_to = "cases") %>%
  #   mutate(virus = fct_relevel(virus, c("S_IFVA", "S_IFVB", "S_RSV", "S_PIV", "S_MPV", "S_sCoV", "S_RV", "S_AdV"))) %>%
  #   # filter(time > 100) %>%
  #   ggplot(., aes(x = time, y = cases, colour = virus)) +
  #   geom_line(alpha = 0.7) +
  #   scale_x_continuous(breaks = seq(1, length(SimResult[, 1]), by = 52)) +
  #   # scale_y_continuous(limits = c(0, 10000)) +
  #   theme_bw()

  # plot(fig1)
  # plot(fig3)
  # plot(fig2)
  # plot(fig4)
  return(list(
    Data = SimResult,
    fig1 = fig1,
    # fig3 = fig3,
    fig2 = fig2 # ,
    # fig4 = fig4
  ))
}

#' @title Find peak of each virus
#'
#' @description Find peak of each virus
#' @Parm RawDat the result of ModelSimCpp
#' @Parm StartTime the start time of virus onset, usually 290
#' @Parm NPI_start the start time of NPI, extracted from Parm, should same as the one used in ModelSimCpp
#' @Parm NPI_end the end time of NPI, extracted from Parm, should same as the one used in ModelSimCpp
#' @Parm Threshold the threshold of virus onset, default is 0.2
#' @Parm span the span of findpeaks function, default is 24
#' @Parm Offset the offset of findpeaks function, avoid the impact of the last wave before NPI (which is lower than the
#'               normal wave), default is 25
#' @Parm Method the method to calculate likelihood, Loss function or just Peak
#' @Parm TargetDat the target data to compare with, get form fread("data/target.csv")(which is the result from review)
#' @return Return a list of peak and time of each virus, or the RMSE of Loss function, or the likelihood of Likelihood function
FindPeak <- function(RawDat, StartTime = 290, NPI_start, NPI_end, Threshold = 0.2, span = 25, Offset = 25,
                     Method = c("Peak", "Loss", "Likelihood"), TargetDat) {
  RawDat <- as.data.table(RawDat)
  RawDat$time <- floor(RawDat$time)
  dat <- RawDat[, lapply(.SD, sum), by = time, .SDcols = IFVA:AdV]

  BeforeNPI <- dat[time > StartTime & time < NPI_start - Offset, ]
  # AfterNPI <- dat[time > NPI_end - Offset, ]
  AfterNPI <- dat[time > NPI_start , ] 
  # 在V5模型中，修改了新的竞争模式，此时一些病毒可能在NPI未完全消失的时候就开始流行，因此放宽了之前在NPI结束时间提前25周寻找峰值的限制，
  # 改为在NPI开始时间后寻找峰值，这样可以避免找不到NPI期间的峰值的问题

  PeaksBeforeNPI <- lapply(BeforeNPI[, 2:9], \(x) x[splus2R::peaks(x, span = span)]) # find peaks before NPI
  PeakThreshold <- lapply(PeaksBeforeNPI, \(x) mean(x) * Threshold) # set threshold for virus onset
  PeaksAfterNPI <- lapply(AfterNPI[, 2:9], \(x) x[splus2R::peaks(x, span = span)]) # find peaks after NPI
  AfterTime <- AfterNPI[, 1]

  CheckOnset <- mapply(function(x, y) x > y, PeaksAfterNPI, PeakThreshold) # CheckOnset is a list of TRUE/FALSE
  FindPeak_after_Identify <- lapply(AfterNPI[, 2:9], splus2R::peaks, span = span) # find peaks after NPI

  Identifier <- as.list(colnames(dat)[-1])
  PeakAndTime <- mapply(function(Onset, PeakDat, PeakDatIdentify, AfterTime, Identifier) {
    if (sum(Onset) == 0) { # Sum of Onset is 0 means no peak found after NPI
      warning(sprintf("No peak found after NPI for %s", Identifier))
      PeakIdentify <- NA
      PeakTime <- tail(AfterTime, 1)
    } else {
      for (i in seq_len(length(Onset))) {
        if (Onset[i] == TRUE) {
          PeakIdentify <- PeakDat[i] # locate peak
          PeakTime <- AfterTime[PeakDatIdentify][i] # locate time
          break
        }
      }
    }
    return(list(PeakIdentify, PeakTime))
  }, CheckOnset, PeaksAfterNPI, FindPeak_after_Identify, AfterTime, Identifier)

  if (Method == "Peak") {
    return(PeakAndTime)
  } else if (Method == "Loss") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[, c("peak", "time") := lapply(.SD, as.numeric), .SDcols = c("peak", "time")][, ":="(PredInterval = (time - NPI_start) * 7, predict = (time - NPI_start) * 7)][, error := (mean - predict)^2]
    Result <- sqrt(sum(CombineTable$error) / 8) # RMSE
    return(Result)
  } else if (Method == "Likelihood") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[, c("peak", "time") := lapply(.SD, as.numeric), .SDcols = c("peak", "time")
                                 ][, ":="(PredInterval = (time - NPI_start) * 7, lambda = 1 / ((time - NPI_start) * 7))
                                   ][, density := dexp(mean, lambda, log = TRUE)]
    Result <- sum(CombineTable$density)
    return(Result)
  }
}

#' @title Run simulation and calculate likelihood
#'
#' @description Run simulation and calculate likelihood
#' @Parm ... same as Model.RunSim and FindPeak
#' @return Return the likelihood of the simulation
Model.RunSim.LLH <- function(Parm, ncores = 6, NPI = FALSE, BaseImmu = FALSE, seeds = NULL, StartTime = 290,
                             Threshold = 0.2, span = 25, Offset = 25, Method = "Likelihood", TargetDat) {
  if (!is.null(seeds)) set.seed(seeds)

  SimResult <- ModelSimCpp(Parm = Parm, ncores = ncores, NPI = NPI, BaseImmu = BaseImmu)
  colnames(SimResult) <- c("time", "IFVA", "IFVB", "RSV", "PIV", "MPV", "sCoV", "RV", "AdV")

  Likelihood <- FindPeak(
    RawDat = SimResult, StartTime = StartTime, NPI_start = Parm$NPI_start, NPI_end = Parm$NPI_end,
    Threshold = Threshold, span = span, Offset = Offset, Method = "Likelihood", TargetDat = Target
  )
  return(Likelihood)
}


Plot.SimResult <- function(data) {
  SummResult <- lapply(1:length(data), function(id) {
    dt <- as.data.table(data[[id]][1])
    setnames(dt, c("time", "IFVA", "IFVB", "RSV", "PIV", "MPV", "sCoV", "RV", "AdV"))
    # 注意data.table的setnames函数名称全部为小写，setNames是用来处理dataframe的
    dt[, MatrixID := id]
    return(dt)
  })
  SummResult <- rbindlist(SummResult)
  # Calculate the median of each time point for each column
  MedianSumm <- SummResult[, lapply(.SD, median), by = time]

  ggplot(SummResult, aes(x = time, y = IFVA, group = MatrixID)) +
    geom_line(alpha = 0.1) +
    geom_line(aes(x = time, y = IFVA), data = MedianSumm, colour = "#650404", linewidth = 1.5) +
    theme_bw()
}

MCMC.Proposal <- function(Parm, step = 5) { # mean = 0, sd = 10
  # ParmReal <- log((1 - Parm) / Parm) / -0.1
  # ParmUpdate <- ParmReal + runif(8, -step, step) #  rnorm(8, mean, sd)
  # return(1 / (1 + exp(-0.1 * ParmUpdate)))
  ParmReal <- asin((Parm - 0.5) * 2) / 0.05
  ParmUpdate <- ParmReal + runif(8, -step, step) #  rnorm(8, mean, sd)
  return(0.5 * sin(0.05 * ParmUpdate) + 0.5)
}

MCMC.MH <- function(Prior, n_iterations, ncores = 4, step = 10, TargetDat = TargetDat) { # mean = 0, sd = 0.5,
  chain <- matrix(NA, nrow = n_iterations, ncol = 8)
  chain[1, ] <- Prior

  current_log_likelihood <- Model.RunSim.LLH(
    Parm = Parameter.Create(
      comp_IFVA = chain[1, 1], comp_IFVB = chain[1, 2], comp_RSV = chain[1, 3], comp_PIV = chain[1, 4],
      comp_MPV = chain[1, 5], comp_sCoV = chain[1, 6], comp_RV = chain[1, 7], comp_AdV = chain[1, 8]
    ), ncores = ncores, NPI = TRUE, StartTime = 290, Threshold = 0.4,
    span = 25, Offset = 25, TargetDat = TargetDat
  )
  pb <- progress_bar$new(total = n_iterations, clear = TRUE, format = "  [:bar] :percent :etas")
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal(Parm = chain[i - 1, ], step = step) # mean = mean, sd = sd
    proposal_log_likelihood <- Model.RunSim.LLH(
      Parm = Parameter.Create(
        comp_IFVA = proposal[1], comp_IFVB = proposal[2], comp_RSV = proposal[3], comp_PIV = proposal[4],
        comp_MPV = proposal[5], comp_sCoV = proposal[6], comp_RV = proposal[7], comp_AdV = proposal[8]
      ), ncores = ncores, NPI = TRUE, StartTime = 290, Threshold = 0.4,
      span = 25, Offset = 25, TargetDat = TargetDat
    )

    print(sprintf("n_iteration is: %d Current LLH is: %f Proposal LLH is: %f", i, current_log_likelihood, proposal_log_likelihood))

    acceptance_ratio <- exp(proposal_log_likelihood - current_log_likelihood)
    if (runif(1) < acceptance_ratio) {
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
    } else {
      chain[i, ] <- chain[i - 1, ]
    }

    print(chain[i, ])
    pb$tick()
  }
  return(chain)
}
