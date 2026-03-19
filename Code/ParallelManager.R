.tail_lines <- function(path, n = 80) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(character(0))
  }
  x <- readLines(path, warn = FALSE)
  tail(x, n)
}

.format_calls <- function(calls = sys.calls()) {
  if (length(calls) == 0) {
    return("")
  }
  paste(
    sprintf(
      "[%02d] %s",
      seq_along(calls),
      vapply(
        calls,
        function(cl) paste(deparse(cl), collapse = ""),
        character(1)
      )
    ),
    collapse = "\n"
  )
}

.mcmc_mgr_state <- function(mgr) attr(mgr, ".state_env", exact = TRUE)


print.mcmc_future_manager <- function(x, ...) {
  st <- .mcmc_mgr_state(x)
  if (is.null(st) || !is.environment(st)) {
    cat("<invalid mcmc_future_manager>\n")
    return(invisible(x))
  }

  cat("<mcmc_future_manager>\n")
  cat(" label      :", st$label, "\n")
  cat(" backend    :", st$backend, "\n")
  cat(" workers    :", st$workers, "\n")
  cat(" debug_mode :", st$debug_mode, "\n")
  cat(" closed     :", st$closed, "\n")
  cat(" jobs       :", length(st$jobs), "\n")
  if (length(st$jobs) > 0) {
    cat(" job_ids    :", paste(names(st$jobs), collapse = ", "), "\n")
  }
  invisible(x)
}


MCMC.Future.ManagerCreate <- function(
  workers = max(1L, future::availableCores() - 1L),
  backend = c("multisession", "multicore", "sequential"),
  debug_mode = FALSE,
  future_debug = debug_mode,
  fallback_to_multisession = TRUE,
  label = "mcmc_manager"
) {
  backend <- match.arg(backend)
  workers <- as.integer(max(1L, workers))

  # Debug mode: force sequential + workers = 1
  if (isTRUE(debug_mode)) {
    backend <- "sequential"
    workers <- 1L
  }

  # Multicore support check
  if (identical(backend, "multicore") && !isTRUE(future::supportsMulticore())) {
    msg <- paste0(
      "backend='multicore' is not supported in this session/OS/IDE. ",
      "Use 'multisession' instead."
    )
    if (isTRUE(fallback_to_multisession)) {
      warning(msg, " Falling back to 'multisession'.", call. = FALSE)
      backend <- "multisession"
    } else {
      stop(msg, call. = FALSE)
    }
  }

  if (identical(backend, "sequential")) {
    workers <- 1L
  }

  st <- new.env(parent = emptyenv())
  st$label <- label
  st$workers <- workers
  st$backend <- backend
  st$debug_mode <- isTRUE(debug_mode)
  st$future_debug <- isTRUE(future_debug)
  st$created_at <- Sys.time()
  st$jobs <- list()
  st$closed <- FALSE

  st$old_plan <- future::plan()
  st$old_future_debug <- getOption("future.debug")

  options(future.debug = st$future_debug)

  switch(
    backend,
    sequential = plan(sequential),
    multisession = plan(multisession, workers = workers),
    multicore = plan(multicore, workers = workers)
  )

  #' Generate a unique job ID based on timestamp and existing jobs
  MakeJobID <- function() {
    idx <- length(st$jobs) + 1L
    id <- sprintf("job_%s_%03d", format(Sys.time(), "%Y%m%d_%H%M%S"), idx)
    while (!is.null(st$jobs[[id]])) {
      idx <- idx + 1L
      id <- sprintf("job_%s_%03d", format(Sys.time(), "%Y%m%d_%H%M%S"), idx)
    }
    return(id)
  }

  safe_value <- function(f) {
    tryCatch(
      future::value(f),
      error = function(e) {
        list(
          ok = FALSE,
          value = NULL,
          error = paste0("Future-level error: ", conditionMessage(e)),
          error_class = paste(class(e), collapse = "/"),
          call_stack = NA_character_,
          warnings = character(0),
          pid = NA_integer_,
          elapsed_sec = NA_real_,
          started_at = as.POSIXct(NA),
          ended_at = as.POSIXct(NA),
          seed = NA_integer_
        )
      }
    )
  }

  #' Get the status of each chain in a job as a data frame
  #' @param job A job object containing futures
  #' @param timeout Time in seconds to wait for each future to resolve when checking status
  #' @return A data frame with columns: chain, resolved, state, ok, pid, elapsed_sec, error
  ChainStatusDataFrame <- function(job, timeout = 0) {
    fs <- job$futures
    nm <- names(fs) %||% paste0("chain_", seq_along(fs))

    rows <- lapply(seq_along(fs), function(i) {
      f <- fs[[i]]
      resolved <- future::resolved(f, timeout = timeout)

      state <- if (resolved) "resolved" else "running"
      ok <- NA
      pid <- NA_integer_
      elapsed_sec <- NA_real_
      err <- NA_character_
      err_cls <- NA_character_
      warn_n <- NA_integer_

      if (resolved) {
        v <- safe_value(f)
        ok <- isTRUE(v$ok)
        state <- if (ok) "success" else "failed"
        pid <- as.integer(v$pid %||% NA_integer_)
        elapsed_sec <- as.numeric(v$elapsed_sec %||% NA_real_)
        err <- as.character(v$error %||% NA_character_)
        err_cls <- as.character(v$error_class %||% NA_character_)
        warn_n <- length(v$warnings %||% character(0))
      }

      data.frame(
        chain = nm[i],
        resolved = resolved,
        state = state,
        ok = ok,
        pid = pid,
        elapsed_sec = elapsed_sec,
        warnings_n = warn_n,
        error_class = err_cls,
        error = err,
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  #' Submit a new MCMC job with multiple chains
  #' @param job_id Optional job ID (auto-generated if NULL)
  #' @param browser_at_start If TRUE, call browser() at the start of each chain (only works in sequential backend)
  #' @param woeker_sleep Number of seconds to sleep before starting the MCMC run in each worker (can be used to stagger start times)
  Submit <- function(
    n_chain = 4L,
    InitialList = NULL,
    n_iterations,
    step = 0.1,
    TargetDat = NULL,
    Prep = NULL,
    Method = c("BetaBinomial", "RatioAbs", "LogDiff", "Dirichlet"),
    BaseParm = Parameter.Create(),
    after = as.Date("2015-08-31"),
    LLHArg = list(),
    PriorArg = list(),
    Adapt = list(),
    InferArg = list(),
    seed = NULL,
    source_r = NULL,
    source_cpp = NULL,
    job_id = NULL,
    browser_at_start = FALSE,
    worker_sleep = 0
  ) {
    if (isTRUE(st$closed)) {
      stop("Manager is shutdown. Create a new manager.")
    }

    Method <- match.arg(Method)
    n_chain <- as.integer(n_chain)
    worker_sleep <- as.numeric(worker_sleep)

    if (is.null(seed)) {
      seed <- 42L
    } else {
      seed <- as.integer(seed)
    }

    if (is.null(job_id)) {
      job_id <- MakeJobID()
    } else {
      job_id <- as.character(job_id)
      if (!is.null(st$jobs[[job_id]])) {
        stop("`job_id` already exists: ", job_id)
      }
    }

    n_virus <- length(BaseParm[["beta0"]])

    if (is.null(Prep)) {
      if (is.null(TargetDat)) {
        stop("Provide either `Prep` or `TargetDat`.")
      }
      Prep <- Inference.Setup(
        TargetDat = TargetDat,
        BaseParm = BaseParm,
        after = after,
        viruses = LLHArg$viruses %||% paste0("Virus_", seq_len(n_virus))
      )
    }

    if (is.null(InitialList)) {
      InitialList <- MCMC.MakeInitialList(
        n_chain = n_chain,
        Initial = NULL,
        Method = Method,
        n_virus = n_virus,
        jitter = 0.3,
        seed = seed
      )
    }
    stopifnot(length(InitialList) == n_chain)

    backend_local <- st$backend

    fs <- lapply(seq_len(n_chain), function(k) {
      chain_seed <- as.integer(seed + k)

      future::future(
        {
          pid <- Sys.getpid()
          t0 <- Sys.time()
          warnings_buf <- character(0)

          if (is.finite(worker_sleep) && worker_sleep > 0) {
            Sys.sleep(worker_sleep)
          }

          if (
            isTRUE(browser_at_start) && identical(backend_local, "sequential")
          ) {
            browser()
          }

          set.seed(chain_seed)

          out <- tryCatch(
            withCallingHandlers(
              {
                if (!is.null(source_r)) {
                  source(source_r, local = TRUE)
                }

                if (!is.null(source_cpp)) {
                  if (!requireNamespace("Rcpp", quietly = TRUE)) {
                    stop(
                      "`source_cpp` is set but package 'Rcpp' is not installed."
                    )
                  }
                  cache_dir <- file.path(
                    tempdir(),
                    sprintf("rcpp_cache_%s_pid%s_chain%s", job_id, pid, k)
                  )
                  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
                  Rcpp::sourceCpp(
                    source_cpp,
                    cacheDir = cache_dir,
                    rebuild = FALSE
                  )
                }

                ans <- MCMC.MH.Adaptive(
                  Initial = InitialList[[k]],
                  n_iterations = n_iterations,
                  step = step,
                  Prep = Prep,
                  Method = Method,
                  BaseParm = BaseParm,
                  after = after,
                  LLHArg = LLHArg,
                  PriorArg = PriorArg,
                  InferArg = InferArg,
                  Adapt = Adapt,
                  verbose = FALSE
                )

                list(
                  ok = TRUE,
                  value = ans,
                  error = NULL,
                  error_class = NA_character_,
                  call_stack = NA_character_
                )
              },
              warning = function(w) {
                warnings_buf <<- c(warnings_buf, conditionMessage(w))
                invokeRestart("muffleWarning")
              }
            ),
            error = function(e) {
              list(
                ok = FALSE,
                value = NULL,
                error = conditionMessage(e),
                error_class = paste(class(e), collapse = "/"),
                call_stack = .format_calls()
              )
            }
          )

          t1 <- Sys.time()
          out$warnings <- warnings_buf
          out$pid <- pid
          out$seed <- chain_seed
          out$started_at <- t0
          out$ended_at <- t1
          out$elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
          out
        },
        seed = chain_seed
      )
    })

    names(fs) <- paste0("chain_", seq_len(n_chain))

    st$jobs[[job_id]] <- list(
      job_id = job_id,
      submitted_at = Sys.time(),
      n_chain = n_chain,
      method = Method,
      seed = seed,
      backend = st$backend,
      debug_mode = st$debug_mode,
      futures = fs,
      collected = FALSE,
      collected_at = as.POSIXct(NA)
    )

    invisible(job_id)
  }

  #' Get the status of one or more jobs as a data frame
  Status <- function(job_id = NULL, timeout = 0) {
    ids <- if (is.null(job_id)) names(st$jobs) else as.character(job_id)

    if (length(ids) == 0) {
      return(data.frame(
        job_id = character(0),
        state = character(0),
        backend = character(0),
        n_chain = integer(0),
        resolved = integer(0),
        success = integer(0),
        failed = integer(0),
        collected = logical(0),
        submitted_at = as.POSIXct(character(0)),
        stringsAsFactors = FALSE
      ))
    }

    rows <- lapply(ids, function(id) {
      job <- st$jobs[[id]]
      if (is.null(job)) {
        return(data.frame(
          job_id = id,
          state = "not_found",
          backend = NA_character_,
          n_chain = NA_integer_,
          resolved = NA_integer_,
          success = NA_integer_,
          failed = NA_integer_,
          collected = NA,
          submitted_at = as.POSIXct(NA),
          stringsAsFactors = FALSE
        ))
      }

      cdf <- chain_status_df(job, timeout = timeout)
      n_chain <- nrow(cdf)
      n_resolved <- sum(cdf$resolved)
      n_success <- sum(cdf$state == "success")
      n_failed <- sum(cdf$state == "failed")

      state <- if (n_resolved < n_chain) {
        if (n_failed > 0) "running_with_error" else "running"
      } else {
        if (n_failed > 0) "finished_with_error" else "finished"
      }

      data.frame(
        job_id = id,
        state = state,
        backend = job$backend %||% st$backend,
        n_chain = n_chain,
        resolved = n_resolved,
        success = n_success,
        failed = n_failed,
        collected = isTRUE(job$collected),
        submitted_at = job$submitted_at,
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  #' Get the status of each chain in a job as a data frame
  ChainStatus <- function(job_id, timeout = 0) {
    job_id <- as.character(job_id)
    job <- st$jobs[[job_id]]
    if (is.null(job)) {
      stop("Job not found: ", job_id)
    }
    ChainStatus(job, timeout = timeout)
  }

  #' Check if all chains in a job are ready (resolved)
  Ready <- function(job_id, timeout = 0) {
    job_id <- as.character(job_id)
    job <- st$jobs[[job_id]]
    if (is.null(job)) {
      stop("Job not found: ", job_id)
    }
    return(all(vapply(
      job$futures,
      future::resolved,
      logical(1),
      timeout = timeout
    )))
  }

  #' Collect results from a completed job, optionally stopping on error and/or removing the job
  collect <- function(job_id, stop_on_error = TRUE, remove = FALSE) {
    job_id <- as.character(job_id)
    job <- st$jobs[[job_id]]
    if (is.null(job)) {
      stop("Job not found: ", job_id)
    }

    vals <- lapply(job$futures, safe_value)
    chain_names <- names(job$futures) %||%
      paste0("chain_", seq_along(job$futures))

    warnings_n <- vapply(
      vals,
      function(x) length(x$warnings %||% character(0)),
      integer(1)
    )
    warnings_txt <- vapply(
      vals,
      function(x) {
        w <- x$warnings %||% character(0)
        if (length(w) == 0) NA_character_ else paste(w, collapse = " | ")
      },
      character(1)
    )

    report <- data.frame(
      chain = chain_names,
      ok = vapply(vals, function(x) isTRUE(x$ok), logical(1)),
      pid = vapply(
        vals,
        function(x) as.integer(x$pid %||% NA_integer_),
        integer(1)
      ),
      elapsed_sec = vapply(
        vals,
        function(x) as.numeric(x$elapsed_sec %||% NA_real_),
        numeric(1)
      ),
      warnings_n = warnings_n,
      warnings = warnings_txt,
      error_class = vapply(
        vals,
        function(x) as.character(x$error_class %||% NA_character_),
        character(1)
      ),
      error = vapply(
        vals,
        function(x) as.character(x$error %||% NA_character_),
        character(1)
      ),
      call_stack = vapply(
        vals,
        function(x) as.character(x$call_stack %||% NA_character_),
        character(1)
      ),
      stringsAsFactors = FALSE
    )

    st$jobs[[job_id]]$collected <- TRUE
    st$jobs[[job_id]]$collected_at <- Sys.time()

    if (stop_on_error && any(!report$ok)) {
      bad <- report[!report$ok, , drop = FALSE]
      lines <- paste0(" - ", bad$chain, ": ", bad$error)
      stop(
        paste0(
          "Some chains failed in ",
          job_id,
          ":\n",
          paste(lines, collapse = "\n")
        ),
        call. = FALSE
      )
    }

    chains <- lapply(vals, `[[`, "value")
    names(chains) <- chain_names
    attr(chains, "report") <- report
    attr(chains, "job_id") <- job_id
    attr(chains, "submitted_at") <- job$submitted_at
    attr(chains, "collected_at") <- Sys.time()

    if (isTRUE(remove)) {
      st$jobs[[job_id]] <- NULL
    }

    chains
  }

  #' Cancel a running job, optionally removing it from the manager
  cancel <- function(job_id, remove = FALSE) {
    job_id <- as.character(job_id)
    job <- st$jobs[[job_id]]
    if (is.null(job)) {
      stop("Job not found: ", job_id)
    }

    has_cancel <- "cancel" %in% getNamespaceExports("future")
    fs <- job$futures
    nm <- names(fs) %||% paste0("chain_", seq_along(fs))

    rows <- lapply(seq_along(fs), function(i) {
      f <- fs[[i]]
      resolved <- future::resolved(f, timeout = 0)
      attempted <- FALSE
      canceled <- FALSE
      msg <- NA_character_

      if (!resolved) {
        attempted <- TRUE
        if (has_cancel) {
          canceled <- tryCatch(
            {
              future::cancel(f)
              TRUE
            },
            error = function(e) {
              msg <<- conditionMessage(e)
              FALSE
            }
          )
          if (isTRUE(canceled)) msg <- "cancel requested"
        } else {
          msg <- "future::cancel() not available in current future version."
        }
      } else {
        msg <- "already resolved"
      }

      data.frame(
        chain = nm[i],
        attempted = attempted,
        canceled = canceled,
        message = msg,
        stringsAsFactors = FALSE
      )
    })

    out <- do.call(rbind, rows)
    rownames(out) <- NULL

    if (isTRUE(remove)) {
      st$jobs[[job_id]] <- NULL
    }
    out
  }

  remove <- function(job_id, force = FALSE) {
    job_id <- as.character(job_id)
    job <- st$jobs[[job_id]]
    if (is.null(job)) {
      stop("Job not found: ", job_id)
    }

    resolved_all <- all(vapply(
      job$futures,
      future::resolved,
      logical(1),
      timeout = 0
    ))
    if (!resolved_all && !force) {
      stop("Job is still running. Use `force = TRUE` to cancel then remove.")
    }

    if (!resolved_all && force) {
      cancel(job_id, remove = FALSE)
    }

    st$jobs[[job_id]] <- NULL
    invisible(TRUE)
  }

  jobs <- function() names(st$jobs)

  config <- function() {
    list(
      label = st$label,
      backend = st$backend,
      workers = st$workers,
      debug_mode = st$debug_mode,
      future_debug = st$future_debug,
      created_at = st$created_at,
      closed = st$closed
    )
  }

  #' Shutdown the manager, optionally canceling running jobs and restoring previous future plan and options
  shutdown <- function(
    cancel_running = FALSE,
    restore_plan = TRUE,
    restore_options = TRUE
  ) {
    if (isTRUE(st$closed)) {
      return(invisible(TRUE))
    }

    ids <- names(st$jobs)
    if (length(ids) > 0) {
      unresolved <- vapply(
        ids,
        function(id) {
          job <- st$jobs[[id]]
          !all(vapply(job$futures, future::resolved, logical(1), timeout = 0))
        },
        logical(1)
      )
      running_ids <- ids[unresolved]

      if (length(running_ids) > 0 && !cancel_running) {
        stop(
          "There are running jobs: ",
          paste(running_ids, collapse = ", "),
          ". Use `cancel_running = TRUE`."
        )
      }

      if (length(running_ids) > 0 && cancel_running) {
        for (id in running_ids) {
          try(cancel(id, remove = FALSE), silent = TRUE)
        }
      }

      st$jobs <- list()
    }

    if (isTRUE(restore_plan) && !is.null(st$old_plan)) {
      future::plan(st$old_plan)
    }

    if (isTRUE(restore_options)) {
      if (is.null(st$old_future_debug)) {
        options(future.debug = NULL)
      } else {
        options(future.debug = st$old_future_debug)
      }
    }

    st$closed <- TRUE
    invisible(TRUE)
  }

  mgr <- list(
    Submit = Submit,
    Status = Status,
    ChainStatus = ChainStatus,
    Ready = Ready,
    collect = collect,
    cancel = cancel,
    remove = remove,
    jobs = jobs,
    config = config,
    shutdown = shutdown
  )

  class(mgr) <- "mcmc_future_manager"
  attr(mgr, ".state_env") <- st
  return(mgr)
}
