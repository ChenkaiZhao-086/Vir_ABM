#' @title Find peak of each virus
#'
#' @description Find peak of each virus
#' @param RawDat the result of ModelSimCpp
#' @param StartTime the start time of virus onset, usually 290
#' @param NPI_start the start time of NPI, extracted from Parm, should same as the one used in ModelSimCpp
#' @param NPI_end the end time of NPI, extracted from Parm, should same as the one used in ModelSimCpp
#' @param Threshold the threshold of virus onset, default is 0.2
#' @param span the span of findpeaks function, default is 24
#' @param Offset the offset of findpeaks function, avoid the impact of the last wave before NPI (which is lower than the
#'               normal wave), default is 25
#' @param Method the method to calculate likelihood, Loss function or just Peak
#' @param TargetDat the target data to compare with, get form fread("data/target.csv")(which is the result from review)
#' @return Return a list of peak and time of each virus, or the RMSE of Loss function, or the likelihood of Likelihood function
FindPeak <- function(
  RawDat,
  StartTime = 290,
  NPI_start,
  NPI_end,
  Threshold = 0.2,
  span = 25,
  Offset = 25,
  Method = c("Peak", "Loss", "Likelihood"),
  TargetDat
) {
  RawDat <- as.data.table(RawDat[-1, ]) # 去掉第一行初始值
  RawDat$time <- floor(RawDat$time)
  dat <- RawDat[, lapply(.SD, sum), by = time, .SDcols = IFVA:AdV]

  BeforeNPI <- dat[time > StartTime & time < NPI_start - Offset, ]
  # AfterNPI <- dat[time > NPI_end - Offset, ]
  AfterNPI <- dat[time > NPI_start, ]
  # 在V5模型中，修改了新的竞争模式，此时一些病毒可能在NPI未完全消失的时候就开始流行，因此放宽了之前在NPI结束时间提前25周寻找峰值的限制，
  # 改为在NPI开始时间后寻找峰值，这样可以避免找不到NPI期间的峰值的问题

  PeaksBeforeNPI <- lapply(BeforeNPI[, 2:9], \(x) {
    x[splus2R::peaks(x, span = span)]
  }) # find peaks before NPI
  PeakThreshold <- lapply(PeaksBeforeNPI, \(x) mean(x) * Threshold) # set threshold for virus onset
  PeaksAfterNPI <- lapply(AfterNPI[, 2:9], \(x) {
    x[splus2R::peaks(x, span = span)]
  }) # find peaks after NPI
  AfterTime <- AfterNPI[, 1]

  CheckOnset <- mapply(
    function(x, y) x > y,
    PeaksAfterNPI,
    PeakThreshold,
    SIMPLIFY = F
  ) # CheckOnset is a list of TRUE/FALSE
  FindPeak_after_Identify <- lapply(
    AfterNPI[, 2:9],
    splus2R::peaks,
    span = span
  ) # find peaks after NPI

  Identifier <- as.list(colnames(dat)[-1])
  PeakAndTime <- mapply(
    function(Onset, PeakDat, PeakDatIdentify, AfterTime, Identifier) {
      if (sum(Onset) == 0) {
        # Sum of Onset is 0 means no peak found after NPI
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
    },
    CheckOnset,
    PeaksAfterNPI,
    FindPeak_after_Identify,
    AfterTime,
    Identifier
  )

  if (Method == "Peak") {
    return(PeakAndTime)
  } else if (Method == "Loss") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[,
      c("peak", "time") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time")
    ][,
      ":="(
        PredInterval = (time - NPI_start) * 7,
        predict = (time - NPI_start) * 7
      )
    ][,
      error := (PredInterval - predict)^2
    ]
    Result <- sqrt(sum(CombineTable$error) / 8) # RMSE
    return(Result)
  } else if (Method == "Likelihood") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[,
      c("peak", "time") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time")
    ][,
      ":="(
        PredInterval = (time - NPI_start) * 7,
        lambda = 1 / ((time - NPI_start) * 7)
      )
    ][,
      density := dexp(mean, lambda, log = TRUE)
    ]
    Result <- sum(CombineTable$density)
    return(Result)
  }
}


FindPeak.V2 <- function(
  RawDat,
  span = 25,
  Method = c("Peak", "Loss", "Likelihood"),
  TargetDat,
  Threshold = 0.2,
  OnsetThreshold = 1000
) {
  RawDat <- as.data.table(RawDat[-1, ])
  RawDat$time <- floor(RawDat$time)
  dat <- RawDat[, lapply(.SD, sum), by = time, .SDcols = IFVA:AdV]

  PeaksAfterNPI <- lapply(dat[, 2:9], \(x) {
    x[splus2R::peaks(x, span = span, endbehavior = 1)]
  })
  # 增加endbehavior=1可以识别在半个span之内的峰值，这样可以识别出某些病毒在一开始就马上出峰以及结束时候的峰值
  PeakThreshold <- lapply(PeaksAfterNPI, \(x) mean(x) * Threshold) # set threshold for virus onset
  AfterTime <- dat[, 1]

  CheckOnset <- mapply(
    function(x, y) x > y & x > OnsetThreshold,
    PeaksAfterNPI,
    PeakThreshold,
    SIMPLIFY = F
  ) # CheckOnset is a list of TRUE/FALSE
  # !!这里如果不添加simplify=F，当每一个PeaksAfterNPI中的元素长度都恰好相等的时候，得到的结果会是一个matrix
  FindPeak_after_Identify <- lapply(
    dat[, 2:9],
    splus2R::peaks,
    span = span,
    endbehavior = 1
  )

  Identifier <- as.list(colnames(dat)[-1])
  PeakAndTime <- mapply(
    function(Onset, PeakDat, PeakDatIdentify, AfterTime, Identifier) {
      if (sum(Onset) == 0) {
        # Sum of Onset is 0 means no peak found after NPI
        warning(sprintf("No peak found after NPI for %s", Identifier))
        PeakIdentify <- NA
        PeakTime <- tail(AfterTime, 1) + 1000 # 对于被压制到没有流行的病毒给一个惩罚
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
    },
    CheckOnset,
    PeaksAfterNPI,
    FindPeak_after_Identify,
    AfterTime,
    Identifier
  )

  if (Method == "Peak") {
    return(PeakAndTime)
  } else if (Method == "Loss") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[,
      c("peak", "time") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time")
    ][,
      ":="(
        PredInterval = (time - NPI_start) * 7,
        predict = (time - NPI_start) * 7
      )
    ][,
      error := (PredInterval - predict)^2
    ]
    Result <- sqrt(sum(CombineTable$error) / 8) # RMSE
    return(Result)
  } else if (Method == "Likelihood") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")

    VirName <- rownames(PeakAndTime)
    PeakAndTime <- as.data.frame(PeakAndTime)
    PeakAndTime$Virus_name <- VirName
    PeakAndTime$RefPeak <- PeakAndTime$time[PeakAndTime$Virus_name == "RSV"] # 参照RSV计算似然用

    CombineTable <- merge(TargetDat, PeakAndTime)

    CombineTable <- CombineTable[,
      c("peak", "time", "RefPeak") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time", "RefPeak")
    ][,
      mu := (time - RefPeak) * 7
    ][,
      density := dlaplace(RefInterval, mu = mu, sigma = 90, log = TRUE)
    ]

    # 根据NPI结束到恢复流行峰值计算Likelihood用
    # CombineTable <- CombineTable[, c("peak", "time") := lapply(.SD, as.numeric), .SDcols = c("peak", "time")
    # ][, ":="(PredInterval = time * 7, lambda = 1 / (time * 7))
    # ][, density := dexp(NPI_Interval, lambda, log = TRUE)]
    Result <- sum(CombineTable$density)
    return(Result)
  }
}


## 增加第二个波次的计算
FindPeak.V3 <- function(
  RawDat,
  span = 25,
  Method = c("Peak", "Loss", "Likelihood"),
  TargetDat,
  Threshold = 0.2,
  OnsetThreshold = 100
) {
  RawDat <- as.data.table(RawDat[-1, ])
  RawDat$time <- floor(RawDat$time)
  dat <- RawDat[, lapply(.SD, sum), by = time, .SDcols = names(RawDat)[-1]]

  PeaksAfterNPI <- lapply(dat[, 2:dim(dat)[2]], \(x) {
    x[splus2R::peaks(x, span = span, endbehavior = 1)]
  })
  # 增加endbehavior=1可以识别在半个span之内的峰值，这样可以识别出某些病毒在一开始就马上出峰以及结束时候的峰值
  PeakThreshold <- lapply(PeaksAfterNPI, \(x) mean(x) * Threshold) # set threshold for virus onset
  AfterTime <- dat[, 1]

  CheckOnset <- mapply(
    function(x, y) x > y & x > OnsetThreshold,
    PeaksAfterNPI,
    PeakThreshold,
    SIMPLIFY = F
  ) # CheckOnset is a list of TRUE/FALSE
  # !!这里如果不添加simplify=F，当每一个PeaksAfterNPI中的元素长度都恰好相等的时候，得到的结果会是一个matrix
  FindPeak_after_Identify <- lapply(
    dat[, 2:dim(dat)[2]],
    splus2R::peaks,
    span = span,
    endbehavior = 1
  )

  Identifier <- as.list(colnames(dat)[-1])
  PeakAndTime <- mapply(
    function(Onset, PeakDat, PeakDatIdentify, AfterTime, Identifier) {
      if (sum(Onset) == 0) {
        # Sum of Onset is 0 means no peak found after NPI
        warning(sprintf("No peak found after NPI for %s", Identifier))
        PeakIdentify <- NA
        PeakTime <- tail(AfterTime, 1) + 1000 # 对于被压制到没有流行的病毒给一个惩罚
      } else {
        for (i in seq_len(length(Onset))) {
          WaveMarker <- which(Onset)[1]
          PeakIdentify <- PeakDat[WaveMarker]
          TimeMarker <- which(PeakDatIdentify)[WaveMarker]
          PeakTime <- AfterTime[TimeMarker]
        }
      }
      return(list(PeakIdentify, PeakTime))
    },
    CheckOnset,
    PeaksAfterNPI,
    FindPeak_after_Identify,
    AfterTime,
    Identifier
  )

  PeakAndTimeSecWave <- mapply(
    function(Onset, PeakDat, PeakDatIdentify, AfterTime, Identifier) {
      if (sum(Onset) == 0) {
        # Sum of Onset is 0 means no peak found after NPI
        warning(sprintf("No peak found after NPI for %s", Identifier))
        PeakDatIdentifySec <- NA
        PeakTimeSec <- tail(AfterTime, 1) + 1000 # 对于被压制到没有流行的病毒给一个惩罚
      } else if (sum(Onset) < 2) {
        warning(sprintf("Only one peak found after NPI for %s", Identifier))
        PeakDatIdentifySec <- NA
        PeakTimeSec <- tail(AfterTime, 1) + 1000
      } else {
        SecWaveMarker <- which(Onset)[2]
        PeakDatIdentifySec <- PeakDat[SecWaveMarker]
        SecPeakMarker <- which(PeakDatIdentify)[SecWaveMarker]
        PeakTimeSec <- AfterTime[SecPeakMarker]
      }
      return(list(PeakDatIdentifySec, PeakTimeSec))
    },
    CheckOnset,
    PeaksAfterNPI,
    FindPeak_after_Identify,
    AfterTime,
    Identifier
  )

  if (Method == "Peak") {
    OutTable <- rbind(PeakAndTime, PeakAndTimeSecWave)
    rownames(OutTable) <- c(
      "Case_Wave1",
      "Time_Wave1",
      "Case_Wave2",
      "Time_Wave2"
    )
    return(OutTable)
  } else if (Method == "Loss") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    CombineTable <- cbind(TargetDat, PeakAndTime)
    CombineTable <- CombineTable[,
      c("peak", "time") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time")
    ][,
      ":="(
        PredInterval = (time - NPI_start) * 7,
        predict = (time - NPI_start) * 7
      )
    ][,
      error := (PredInterval - predict)^2
    ]
    Result <- sqrt(sum(CombineTable$error) / 8) # RMSE
    return(Result)
  } else if (Method == "Likelihood") {
    PeakAndTime <- t(PeakAndTime)
    colnames(PeakAndTime) <- c("peak", "time")
    VirName <- rownames(PeakAndTime)
    PeakAndTime <- as.data.frame(PeakAndTime)

    PeakAndTime$Virus_name <- VirName
    PeakAndTime$RefPeak <- PeakAndTime$time[PeakAndTime$Virus_name == "RSV"] # 参照RSV计算似然用

    CombineTable <- merge(TargetDat[Index_of_Wave == 1, ], PeakAndTime)

    CombineTable <- CombineTable[,
      c("peak", "time", "RefPeak") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time", "RefPeak")
    ][,
      mu := (time - RefPeak) * 7
    ][,
      density := dlaplace(RefInterval, mu = mu, sigma = 90, log = TRUE)
    ]

    ## SecWave
    PeakAndTimeSecWave <- t(PeakAndTimeSecWave)
    colnames(PeakAndTimeSecWave) <- c("peak", "time")
    VirName <- rownames(PeakAndTimeSecWave)
    PeakAndTimeSecWave <- as.data.frame(PeakAndTimeSecWave)
    PeakAndTimeSecWave$time <- as.numeric(PeakAndTimeSecWave$time) -
      as.numeric(PeakAndTime$time) +
      1

    PeakAndTimeSecWave$Virus_name <- VirName
    # PeakAndTimeSecWave$RefPeak <- PeakAndTime$time[PeakAndTime$Virus_name == "RSV"] # 这里也要参照第一次的RSV
    CombineTableSecWave <- merge(
      TargetDat[Index_of_Wave == 2, ],
      PeakAndTimeSecWave
    )

    CombineTableSecWave <- CombineTableSecWave[,
      c("peak", "time", "RefInterval") := lapply(.SD, as.numeric),
      .SDcols = c("peak", "time", "RefInterval")
    ][,
      lambda := 1 / (time * 7)
    ][,
      density := dexp(RefInterval, rate = lambda, log = TRUE)
    ]

    # 根据NPI结束到恢复流行峰值计算Likelihood用
    # CombineTable <- CombineTable[, c("peak", "time") := lapply(.SD, as.numeric), .SDcols = c("peak", "time")
    # ][, ":="(PredInterval = time * 7, lambda = 1 / (time * 7))
    # ][, density := dexp(NPI_Interval, lambda, log = TRUE)]
    Result <- sum(CombineTable$density) + sum(CombineTableSecWave$density)
    return(Result)
  }
}


MCMC.Proposal <- function(Parm, step = 5) {
  # mean = 0, sd = 10
  # ParmReal <- log((1 - Parm) / Parm) / -0.1
  # ParmUpdate <- ParmReal + runif(8, -step, step) #  rnorm(8, mean, sd)
  # return(1 / (1 + exp(-0.1 * ParmUpdate)))
  ParmReal <- asin((Parm - 0.5) * 2) / 0.05
  ParmUpdate <- ParmReal + runif(8, -step, step) #  rnorm(8, mean, sd)
  return(0.5 * sin(0.05 * ParmUpdate) + 0.5)
}

MCMC.Proposal.Generator <- function(n) {
  data <- rexp(n, 1)
  geom_mean <- exp(mean(log(data))) # 计算原始数据的几何均值
  adjusted_data <- data / geom_mean # 调整数据以确保其几何均值为1
  return(adjusted_data)
}


MCMC.Proposal.V2 <- function(Parm, step1 = 1, step2 = 5) {
  ParmReal_comp <- log(Parm[1:8]) # log((Parm[1:8]-1)*(exp(0.5)-1)/999+1)/0.5
  ParmUpdate_comp <- ParmReal_comp + runif(8, -step1, step1) # -0.005, 0.005
  # 这里的ParmMean是为了保证参数更新后，所有参数的几何均值和原来一致。
  # 因为后面的参数要进行指数转换，所以这里直接计算均值。
  # 这里要保证所有的几何均值与更新前一致，是因为，如果不这样做多个MCMC链可能会在不同的位置收敛，
  # 这样看起来的结果可能是不收敛。而经过这样的转换后所有的参数都围绕1进行波动，同时也保留了参数之间的差异
  ParmMean <- mean(ParmUpdate_comp)
  ParmUpdate_comp <- ParmUpdate_comp - ParmMean
  NewCompParm <- exp(ParmUpdate_comp) # 1+(999*(exp(0.5*ParmUpdate_comp)-1))/(exp(0.5)-1)

  # Base immunity range: 0-1
  # ParmReal_pop <- asin((Parm[9] - 0.5) / 0.5) / 0.05
  # ParmUpdate_pop <- ParmReal_pop + runif(1, -step2, step2)
  # NewPopParm <- 0.5 * sin(0.05 * ParmUpdate_pop) + 0.5
  # Base immunity range: 0.2-0.6
  ParmReal_pop <- asin((Parm[9] - 0.4) / 0.2) / 0.05
  ParmUpdate_pop <- ParmReal_pop + runif(1, -step2, step2)
  NewPopParm <- 0.2 * sin(0.05 * ParmUpdate_pop) + 0.4

  # ParmReal_Penal <- asin((Parm[10] - 0.5) * 2) / 0.05
  # ParmUpdate_Penal <- ParmReal_Penal + runif(1, -step, step)
  # NewPenalParm <- 0.5 * sin(0.05 * ParmUpdate_Penal) + 0.5
  NewParm <- c(NewCompParm, NewPopParm) # , NewPenalParm

  return(NewParm)
}


MCMC.MH <- function(
  Prior,
  n_iterations,
  ncores = 4,
  step = 10,
  TargetDat = TargetDat,
  Threshold = 0.4
) {
  # mean = 0, sd = 0.5,
  chain <- matrix(NA, nrow = n_iterations, ncol = 8)
  chain[1, ] <- Prior

  current_log_likelihood <- Model.RunSim.LLH(
    Parm = Parameter.Create(
      comp_IFVA = chain[1, 1],
      comp_IFVB = chain[1, 2],
      comp_RSV = chain[1, 3],
      comp_PIV = chain[1, 4],
      comp_MPV = chain[1, 5],
      comp_sCoV = chain[1, 6],
      comp_RV = chain[1, 7],
      comp_AdV = chain[1, 8]
    ),
    ncores = ncores,
    NPI = TRUE,
    StartTime = 290,
    Threshold = Threshold,
    span = 25,
    Offset = 25,
    TargetDat = TargetDat
  )
  pb <- progress_bar$new(
    total = n_iterations,
    clear = TRUE,
    format = "  [:bar] :percent :etas"
  )
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal(Parm = chain[i - 1, ], step = step) # mean = mean, sd = sd
    proposal_log_likelihood <- Model.RunSim.LLH(
      Parm = Parameter.Create(
        comp_IFVA = proposal[1],
        comp_IFVB = proposal[2],
        comp_RSV = proposal[3],
        comp_PIV = proposal[4],
        comp_MPV = proposal[5],
        comp_sCoV = proposal[6],
        comp_RV = proposal[7],
        comp_AdV = proposal[8]
      ),
      ncores = ncores,
      NPI = TRUE,
      StartTime = 290,
      Threshold = Threshold,
      span = 25,
      Offset = 25,
      TargetDat = TargetDat
    )

    Info <- sprintf(
      "n_iteration is: %d Current LLH is: %f Proposal LLH is: %f",
      i,
      current_log_likelihood,
      proposal_log_likelihood
    )
    CLI.Print(Info)

    acceptance_ratio <- min(
      1,
      exp(proposal_log_likelihood - current_log_likelihood)
    )
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


MCMC.MH.V2 <- function(
  Prior,
  n_iterations,
  ncores = 4,
  step1 = 0.5,
  step2 = 10,
  TargetDat = TargetDat,
  Threshold = 0.2,
  OnsetThreshold = 1000
) {
  chain <- matrix(NA, nrow = n_iterations, ncol = 9)
  chain[1, ] <- Prior

  current_log_likelihood <- Model.RunSim.LLH.V2(
    Parm = Parameter.Create(
      comp_IFVA = chain[1, 1],
      comp_IFVB = chain[1, 2],
      comp_RSV = chain[1, 3],
      comp_PIV = chain[1, 4],
      comp_MPV = chain[1, 5],
      comp_sCoV = chain[1, 6],
      comp_RV = chain[1, 7],
      comp_AdV = chain[1, 8],
      base_immune = chain[1, 9],
      years = 4
    ),
    ncores = ncores,
    Threshold = Threshold,
    OnsetThreshold = OnsetThreshold,
    span = 25,
    TargetDat = TargetDat
  )

  pb <- progress_bar$new(
    total = n_iterations,
    clear = TRUE,
    format = "  [:bar] :percent :etas"
  )
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal.V2(
      Parm = chain[i - 1, ],
      step1 = step1,
      step2 = step2
    ) # mean = mean, sd = sd
    proposal_log_likelihood <- Model.RunSim.LLH.V2(
      Parm = Parameter.Create(
        comp_IFVA = proposal[1],
        comp_IFVB = proposal[2],
        comp_RSV = proposal[3],
        comp_PIV = proposal[4],
        comp_MPV = proposal[5],
        comp_sCoV = proposal[6],
        comp_RV = proposal[7],
        comp_AdV = proposal[8],
        base_immune = proposal[9],
        years = 4
      ),
      ncores = ncores,
      Threshold = Threshold,
      OnsetThreshold = OnsetThreshold,
      span = 25,
      TargetDat = TargetDat
    )

    Info <- sprintf(
      "n_iteration is: %d Current LLH is: %f Proposal LLH is: %f",
      i,
      current_log_likelihood,
      proposal_log_likelihood
    )
    CLI.Print(Info)

    acceptance_ratio <- min(
      1,
      exp(proposal_log_likelihood - current_log_likelihood)
    )
    # MH算法的正统写法，直接使用后面的exp部分也正确，因为propose比当前的状态高，带有最小值的算法接受的概率是1，不带有的算法接受概率是1，不会影响实际结果
    if (runif(1) < acceptance_ratio) {
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
    print(proposal)
    print(chain[i, ])
    pb$tick()
  }
  return(chain)
}

#' @title Run simulation and calculate likelihood with adaptive step
#'
#' @description
#' @param
#' @return
MCMC.MH.V3 <- function(
  Prior,
  n_iterations,
  ncores = 4,
  step1 = 0.5,
  step2 = 10,
  TargetDat = TargetDat,
  Threshold = 0.2,
  OnsetThreshold = 1000,
  TargetAcceptance_low = 0.23,
  TargetAcceptance_high = 0.5,
  AdaptInterval = 100
) {
  chain <- matrix(NA, nrow = n_iterations, ncol = 9)
  chain[1, ] <- Prior
  accepted <- 0

  current_log_likelihood <- Model.RunSim.LLH.V2(
    Parm = Parameter.Create(
      comp_IFVA = chain[1, 1],
      comp_IFVB = chain[1, 2],
      comp_RSV = chain[1, 3],
      comp_PIV = chain[1, 4],
      comp_MPV = chain[1, 5],
      comp_sCoV = chain[1, 6],
      comp_RV = chain[1, 7],
      comp_AdV = chain[1, 8],
      base_immune = chain[1, 9],
      years = 4
    ),
    ncores = ncores,
    Threshold = Threshold,
    OnsetThreshold = OnsetThreshold,
    span = 25,
    TargetDat = TargetDat
  )

  pb <- progress_bar$new(
    total = n_iterations,
    clear = TRUE,
    format = "  [:bar] :percent :etas"
  )
  pb$tick()
  for (i in 2:n_iterations) {
    proposal <- MCMC.Proposal.V2(
      Parm = chain[i - 1, ],
      step1 = step1,
      step2 = step2
    ) # mean = mean, sd = sd
    proposal_log_likelihood <- Model.RunSim.LLH.V2(
      Parm = Parameter.Create(
        comp_IFVA = proposal[1],
        comp_IFVB = proposal[2],
        comp_RSV = proposal[3],
        comp_PIV = proposal[4],
        comp_MPV = proposal[5],
        comp_sCoV = proposal[6],
        comp_RV = proposal[7],
        comp_AdV = proposal[8],
        base_immune = proposal[9],
        years = 4
      ),
      ncores = ncores,
      Threshold = Threshold,
      OnsetThreshold = OnsetThreshold,
      span = 25,
      TargetDat = TargetDat
    )

    Info <- sprintf(
      "n_iteration is: %d Current LLH is: %f Proposal LLH is: %f",
      i,
      current_log_likelihood,
      proposal_log_likelihood
    )
    CLI.Print(Info)

    acceptance_ratio <- min(
      1,
      exp(proposal_log_likelihood - current_log_likelihood)
    )
    if (runif(1) < acceptance_ratio) {
      chain[i, ] <- proposal
      current_log_likelihood <- proposal_log_likelihood
      accepted <- accepted + 1
    } else {
      chain[i, ] <- chain[i - 1, ]
    }
    print(proposal)
    print(chain[i, ])
    pb$tick()

    # Adjust propose step every adapt interval iterations based on acceptance rate
    if (i %% AdaptInterval == 0) {
      acceptance_rate <- accepted / AdaptInterval
      if (acceptance_rate < TargetAcceptance_low) {
        step1 <- step1 * 0.9
        # step2 <- step2 * 0.9
      } else if (acceptance_rate > TargetAcceptance_high) {
        step1 <- step1 * 1.1
        # step2 <- step2 * 1.1
      }
      accepted <- 0
      print(acceptance_rate)
      print(step1)
    }
  }
  return(chain)
}
