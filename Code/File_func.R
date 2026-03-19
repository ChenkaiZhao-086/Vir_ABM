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
