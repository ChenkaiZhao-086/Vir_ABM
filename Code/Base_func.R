#' @title Print message to command line interface
#'
#' @description Command line interface: print message
#'
#' @param ... (parts of the) message to print
#' @param WARNING boolean, to print the message in red
#'
#' @keywords internal
CLI.Print <- function(..., WARNING = FALSE, FORCED = FALSE) {
    # get function arguments
    function_arguments <- as.list(match.call(expand.dots = FALSE))$...
    # return: CLI.Print(... = pairlist("Create folder:", FolderPath))

    # get function-call environment (to retrieve variable from that environment)
    pf <- parent.frame()
    # 每当在R中调用一个函数时，都会创建一个新的环境来存放该函数的局部变量。这个新环境是当前环境的子环境。
    # parent.frame() 返回的是这个子环境的父环境，即调用当前函数所在的环境。在这个函数中是调用CLI.Print创造的环境

    # parse list => make character vector
    f_out <- " "
    for (i in 1:length(function_arguments)) {
        f_out <- cbind(f_out, eval(unlist(function_arguments[[i]]), envir = pf))
    }

    # add a space to each function arguments
    function_arguments <- paste(f_out, collapse = " ")

    # set text color: black (default) or red (warning)
    WordColorNormal <- "\033[0;36m"
    WordColorWarning <- "\033[0;31m"
    text_color <- ifelse(WARNING, WordColorWarning, WordColorNormal)

    # print time + arguments (without spaces)
    cli_out <- paste0(
        c(
            'echo "',
            text_color,
            "[",
            format(Sys.time(), "%H:%M:%S"),
            "]",
            function_arguments,
            WordColorNormal,
            '"'
        ),
        collapse = ""
    )
    # 'echo "'得到的结果是"echo \""
    # 输出的结构: echo \"[11:21:27]  Create folder: Output/20240129/\"
    # echo后要有反斜线，输出内容用引号扩住，最后要加正斜线

    # print if function is called by master-node or first slave
    if (
        !exists("ParallelNodesInfo") ||
            Sys.getpid() == ParallelNodesInfo$pid_master ||
            FORCED
    ) {
        system(cli_out)
    }

    # add to R warnings
    if (WARNING) {
        cli_warning <- paste0(
            c(
                text_color,
                "[",
                format(Sys.time(), "%H:%M:%S"),
                "]",
                function_arguments,
                WordColorNormal
            ),
            collapse = ""
        )
        warning(cli_warning, call. = FALSE, immediate. = FALSE)
    }
}


cli_progress <- function(i_current, i_total, time_stamp_loop) {
    # print if function is called by first node
    if (
        exists("ParallelNodesInfo") &&
            Sys.getpid() == ParallelNodesInfo$pid_slave1
    ) {
        # calculate progress
        progress_scen <- floor(i_current / i_total * 100)
        progress_time <- round(
            difftime(Sys.time(), time_stamp_loop, units = "min"),
            digits = 1
        )

        # estimate remaining time (after 15%)
        time_label <- ""
        if (progress_scen > 15 & progress_scen < 99) {
            estim_time <- round(
                progress_time / progress_scen * (100 - progress_scen),
                digits = 1
            )
            if (estim_time < 1) {
                estim_time <- "<1"
            }
            time_label <- paste0("[", estim_time, " min remaining]")
        }

        CLI.Print(
            "RUNNING...",
            i_current,
            "/",
            i_total,
            time_label,
            FORCED = TRUE
        )
    }
}


#' @title Create one Main folder for Result output
#'
#' @description Create one Main folder for output. The folder name is current date or specific name
File.CreateMainFolder <- function(path, FolderName = NULL, Date = NULL) {
    if (!is.null(FolderName)) {
        if (!file.exists(paste0(path, FolderName))) {
            dir.create(paste0(path, FolderName))
        }
        FolderPath <- paste0(path, FolderName, "/")
    }

    if (!is.null(Date)) {
        CurrentDate <- format(Sys.time(), "%Y%m%d")
        if (!file.exists(paste0(path, CurrentDate))) {
            dir.create(paste0(path, CurrentDate))
        }
        FolderPath <- paste0(path, CurrentDate, "/")
    }

    CLI.Print("Create folder:", FolderPath)

    return(FolderPath)
}

#' @title Create one sub folder for Result output
#'
#' @description Create one sub folder for Result output.
File.CreateSubFolder <- function(path, SubFolderName) {
    if (!file.exists(paste0(path, SubFolderName))) {
        dir.create(paste0(path, SubFolderName))
    }
    return(paste0(path, SubFolderName, "/"))
}


Parallel.Regist <- function(ncores = NULL, seed = NULL) {
    ## SETUP PARALLEL NODES
    # note: they will be removed after 600 seconds inactivity
    if (is.null(ncores)) {
        ncores <- detectCores() - 2
    }
    par_cluster <- makeCluster(ncores, cores = ncores) # , timeout = 600
    registerDoParallel(par_cluster)

    # set seed for each node
    if (!is.null(seed)) {
        clusterSetRNGStream(par_cluster, seed)
    }

    # store the process id (pid) of the first slave
    pid_slave1 <- clusterEvalQ(par_cluster, {
        # 用于在集群的每个节点上执行相同的表达式
        Sys.getpid() # 这是一个基础 R 函数，用于获取当前 R 会话的进程 ID
    })[[1]]

    # CREATE GLOBAL VARIABLE
    ParallelNodesInfo <<- list(
        # 节点信息
        par_cluster = par_cluster,
        pid_master = Sys.getpid(),
        pid_slave1 = pid_slave1
    )

    CLI.Print("START PARALLEL WORKERS")
}


#' @param Func_Vir A list, The name of functions and virables that need to be exported
Parallel.Import <- function(Func_Vir) {
    clusterExport(ParallelNodesInfo[[1]], Func_Vir)

    clusterEvalQ(ParallelNodesInfo[[1]], library("Rcpp"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("RcppArmadillo"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("data.table"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("tidyverse"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("extraDistr"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("progress"))
    clusterEvalQ(ParallelNodesInfo[[1]], library("truncnorm"))

    # clusterEvalQ(ParallelNodesInfo[[1]], library("BaoJiaoZa"))
    clusterEvalQ(ParallelNodesInfo[[1]], {
        library("deSolve")
        sourceCpp("Code/ParmInferenceCpp_V4.cpp")
    }) # 在使用c++的时候，要在传递的时候直接把代码编译进去
}


Parallel.Stop <- function() {
    ## CLOSE NODES AND NODE INFO
    if (exists("ParallelNodesInfo")) {
        CLI.Print("STOP PARALLEL WORKERS")

        stopCluster(ParallelNodesInfo$par_cluster)
        rm(ParallelNodesInfo, envir = .GlobalEnv) # REMOVE GLOBAL VARIABLE
    }
}


Parallel.Check <- function() {
    if (!exists("ParallelNodesInfo")) {
        Parallel.Regist()
    } else if (
        !any(grepl(
            ParallelNodesInfo$pid_slave1,
            system("ps -A", intern = TRUE)
        ))
    ) {
        Parallel.Regist()
    }
}


#' @title Create Parmeters for model
#'
#' @description Reset cluster and start a new cluster
# reset parallel workers
Parallel.Reset <- function() {
    Parallel.Stop()
    gc()
    Parallel.Regist()
}


#' @title Create Parmeters for model
#'
#' @description Create Parmeters for model
Parameter.Create <- function(
    # Model base
    num_of_agent = 100000,
    dt = 1,
    years = 16,
    year_start = "1999-01-01",
    year_end = "2024-03-31",
    initial_seeds = 10,
    added_cases = 2, # 每个病毒随机每天新引入的I，这里是每dt2个，相当于一周4个
    Penal = 0.5, # 感染多病毒的惩罚系数
    # R0
    R0_IFVA = 1.41,
    R0_IFVB = 1.07,
    R0_RSV = 1.7, # 3, #
    R0_RV = 1.88,
    # Transmission rate
    # beta_IFVA = 1,
    # beta_IFVB = 1,
    # beta_RSV = 1,
    # beta_RV = 1,
    # Seasonal force
    beta_seasonal = 0.2,
    phi = 340,
    beta_amplify = 1,
    # Duration of infectious
    gamma_IFVA = 1 / 6,
    gamma_IFVB = 1 / 4,
    gamma_RSV = 1 / 7.4,
    gamma_RV = 1 / 10.9,
    # Duration of immunity of each virus
    omega_IFVA = 1 / (365 * 2),
    omega_IFVB = 1 / 424.1,
    omega_RSV = 1 / 358.9,
    omega_RV = 1 / 36.5,
    # Virus competition
    comp_IFVA = 1,
    comp_IFVB = 1,
    comp_RSV = 1,
    comp_RV = 1,
    # NPI susceptibility
    NPISes_IFVA = 1,
    NPISes_IFVB = 1,
    NPISes_RSV = 1,
    NPISes_RV = 1,
    # Base immunity
    base_immune = 0.4,
    # NPI
    NPI = TRUE,
    NPI_value = c(0.8)
    # NPI_start = c("2020-03-10"), # not use this param any more
    # NPI_end = c("2020-04-10"), # not use this param any more
    # NPI_value = c(0), # range[0,1] 1 means totally no transmission(Full NPI)
    # decay_coef = c(0.01) # range[0,Inf] 0 means no decay
) {
    beta_IFVA <- (R0_IFVA * gamma_IFVA) / (1 + beta_seasonal)
    beta_IFVB <- (R0_IFVB * gamma_IFVB) / (1 + beta_seasonal)
    beta_RSV <- (R0_RSV * gamma_RSV) / (1 + beta_seasonal)
    beta_RV <- (R0_RV * gamma_RV) / (1 + beta_seasonal)
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
            Penal = Penal,

            # Seasonal force
            beta_seasonal = beta_seasonal,
            phi = phi,
            beta_amplify = beta_amplify,

            # Transmission rate
            # beta_IFVA = beta_IFVA,
            # beta_IFVB = beta_IFVB,
            # beta_RSV = beta_RSV,
            # beta_RV = beta_RV,
            beta0 = c(beta_IFVA, beta_IFVB, beta_RSV, beta_RV),

            # Duration of infectious
            # gamma_IFVA = gamma_IFVA,
            # gamma_IFVB = gamma_IFVB,
            # gamma_RSV = gamma_RSV,
            # gamma_RV = gamma_RV,
            gamma = c(gamma_IFVA, gamma_IFVB, gamma_RSV, gamma_RV),

            # Duration of immunity of each virus
            # omega_IFVA = omega_IFVA,
            # omega_IFVB = omega_IFVB,
            # omega_RSV = omega_RSV,
            # omega_RV = omega_RV,
            omega = c(omega_IFVA, omega_IFVB, omega_RSV, omega_RV),

            # Virus competition
            # comp_IFVA = comp_IFVA,
            # comp_IFVB = comp_IFVB,
            # comp_RSV = comp_RSV,
            # comp_RV = comp_RV,
            comp = c(comp_IFVA, comp_IFVB, comp_RSV, comp_RV),

            # NPI susceptibility
            NPISes_IFVA = NPISes_IFVA,
            NPISes_IFVB = NPISes_IFVB,
            NPISes_RSV = NPISes_RSV,
            NPISes_RV = NPISes_RV,

            # NPI
            NPI = NPI,
            # NPI_start = as.numeric(as.Date(NPI_start)),
            # NPI_end = as.numeric(as.Date(NPI_end)),
            NPI_value = NPI_value,
            # decay_coef = decay_coef,

            # Base immunity
            base_immune = base_immune,
            R0_IFVA = R0_IFVA,
            R0_IFVB = R0_IFVB,
            R0_RSV = R0_RSV,
            R0_RV = R0_RV
        )
    return(Parmeters)
}

# Parameter.Create <- function(
#     # Model base
#     num_of_agent = 100000,
#     dt = 1,
#     years = 16,
#     year_start = "1999-01-01",
#     year_end = "2024-03-31",
#     initial_seeds = 10,
#     added_cases = 2, # 每个病毒随机每天新引入的I，这里是每dt2个，相当于一周4个
#     Penal = 0.5, # 感染多病毒的惩罚系数
#     # R0
#     R0_IFVA = 1.41,
#     R0_IFVB = 1.07,
#     R0_RSV = 3.0,
#     R0_RV = 1.88,
#     R0_v1 = 1.41,
#     R0_v2 = 1.07,
#     R0_v3 = 3.0,
#     R0_v4 = 1.88,
#     # Transmission rate
#     # beta_IFVA = 1,
#     # beta_IFVB = 1,
#     # beta_RSV = 1,
#     # beta_RV = 1,
#     # Seasonal force
#     beta_seasonal = 0.2,
#     phi = 330,
#     beta_amplify = 1,
#     # Duration of infectious
#     gamma_IFVA = 1 / 6,
#     gamma_IFVB = 1 / 4,
#     gamma_RSV = 1 / 6.4,
#     gamma_RV = 1 / 10.9,
#     gamma_v1 = 1 / 6,
#     gamma_v2 = 1 / 4,
#     gamma_v3 = 1 / 6.4,
#     gamma_v4 = 1 / 10.9,
#     # Duration of immunity of each virus
#     omega_IFVA = 1 / (365 * 2),
#     omega_IFVB = 1 / 424.1,
#     omega_RSV = 1 / 358.9,
#     omega_RV = 1 / 36.5,
#     omega_v1 = 1 / (365 * 2),
#     omega_v2 = 1 / 424.1,
#     omega_v3 = 1 / 358.9,
#     omega_v4 = 1 / 36.5,
#     # Virus competition
#     comp_IFVA = 1,
#     comp_IFVB = 1,
#     comp_RSV = 1,
#     comp_RV = 1,
#     comp_v1 = 1,
#     comp_v2 = 1,
#     comp_v3 = 1,
#     comp_v4 = 1,
#     # NPI susceptibility
#     NPISes_IFVA = 1,
#     NPISes_IFVB = 1,
#     NPISes_RSV = 1,
#     NPISes_RV = 1,
#     NPISes_v1 = 1,
#     NPISes_v2 = 1,
#     NPISes_v3 = 1,
#     NPISes_v4 = 1,
#     # Base immunity
#     base_immune = 0.4,
#     # NPI
#     NPI = TRUE,
#     NPI_value = c(0.8)
#     # NPI_start = c("2020-03-10"), # not use this param any more
#     # NPI_end = c("2020-04-10"), # not use this param any more
#     # NPI_value = c(0), # range[0,1] 1 means totally no transmission(Full NPI)
#     # decay_coef = c(0.01) # range[0,Inf] 0 means no decay
#     ) {
#     beta_IFVA <- (R0_IFVA * gamma_IFVA) / (1 + beta_seasonal)
#     beta_IFVB <- (R0_IFVB * gamma_IFVB) / (1 + beta_seasonal)
#     beta_RSV <- (R0_RSV * gamma_RSV) / (1 + beta_seasonal)
#     beta_RV <- (R0_RV * gamma_RV) / (1 + beta_seasonal)
#     beta_v1 <- (R0_v1 * gamma_v1) / (1 + beta_seasonal)
#     beta_v2 <- (R0_v2 * gamma_v2) / (1 + beta_seasonal)
#     beta_v3 <- (R0_v3 * gamma_v3) / (1 + beta_seasonal)
#     beta_v4 <- (R0_v4 * gamma_v4) / (1 + beta_seasonal)
#     Parmeters <-
#         list(
#             # Model base
#             num_of_agent = num_of_agent,
#             dt = dt,
#             years = years,
#             year_start = year_start,
#             year_end = year_end,
#             initial_seeds = initial_seeds,
#             added_cases = added_cases,
#             Penal = Penal,

#             # Seasonal force
#             beta_seasonal = beta_seasonal,
#             phi = phi,
#             beta_amplify = beta_amplify,

#             # Transmission rate
#             beta_IFVA = beta_IFVA,
#             beta_IFVB = beta_IFVB,
#             beta_RSV = beta_RSV,
#             beta_RV = beta_RV,
#             beta_v1 = beta_v1,
#             beta_v2 = beta_v2,
#             beta_v3 = beta_v3,
#             beta_v4 = beta_v4,

#             # Duration of infectious
#             gamma_IFVA = gamma_IFVA,
#             gamma_IFVB = gamma_IFVB,
#             gamma_RSV = gamma_RSV,
#             gamma_RV = gamma_RV,
#             gamma_v1 = gamma_v1,
#             gamma_v2 = gamma_v2,
#             gamma_v3 = gamma_v3,
#             gamma_v4 = gamma_v4,

#             # Duration of immunity of each virus
#             omega_IFVA = omega_IFVA,
#             omega_IFVB = omega_IFVB,
#             omega_RSV = omega_RSV,
#             omega_RV = omega_RV,
#             omega_v1 = omega_v1,
#             omega_v2 = omega_v2,
#             omega_v3 = omega_v3,
#             omega_v4 = omega_v4,

#             # Virus competition
#             comp_IFVA = comp_IFVA,
#             comp_IFVB = comp_IFVB,
#             comp_RSV = comp_RSV,
#             comp_RV = comp_RV,
#             comp_v1 = comp_v1,
#             comp_v2 = comp_v2,
#             comp_v3 = comp_v3,
#             comp_v4 = comp_v4,

#             # NPI susceptibility
#             NPISes_IFVA = NPISes_IFVA,
#             NPISes_IFVB = NPISes_IFVB,
#             NPISes_RSV = NPISes_RSV,
#             NPISes_RV = NPISes_RV,
#             NPISes_v1 = NPISes_v1,
#             NPISes_v2 = NPISes_v2,
#             NPISes_v3 = NPISes_v3,
#             NPISes_v4 = NPISes_v4,

#             # NPI
#             NPI = NPI,
#             # NPI_start = as.numeric(as.Date(NPI_start)),
#             # NPI_end = as.numeric(as.Date(NPI_end)),
#             NPI_value = NPI_value,
#             # decay_coef = decay_coef,

#             # Base immunity
#             base_immune = base_immune
#         )

#     return(Parmeters)
# }
