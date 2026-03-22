library(readxl)
library(tidyverse)
library(data.table)
# library(foreach)
# library(doParallel)
library(Rcpp)
# library(RcppArmadillo)
library(deSolve)
library(progress)
library(extraDistr)
library(truncnorm)
library(coda)
library(bayesplot)
library(future)
library(future.apply)


# packages <- c(
#   "readxl",
#   "tidyverse",
#   "data.table",
#   "foreach",
#   "doParallel",
#   "Rcpp",
#   "RcppArmadillo",
#   "progress",
#   "extraDistr",
#   "truncnorm",
#   "coda"
# )
# for (i in packages) {
#   suppressPackageStartupMessages(library(
#     i,
#     character.only = TRUE,
#     quietly = TRUE
#   ))
# }

# sourceCpp("Code/ModelSimCpp_V9.cpp")
sourceCpp("Code/ParmInferenceCpp_V6.cpp")

source("Code/File_func.R")
source("Code/Parm_func.R")
source("Code/Ode_func.R")
source("Code/MCMC.R")
source("Code/ParallelManager.R")

source("Code/func.R")
# source("Code/ABM_func.R")
# source("Code/BM_PA_func.R")

### Prepaer observed data and Google mobility report data
# source("Code/CanadaData.R")
# source("Code/CanadaGMR.R")

###
load("Data/HQDat.RData")
load("Data/HQGMR.RData")

###
# R0可以查一下IV，RSV等的SAR（二代发病率），SAR在一定程度下可以
# re-infection的时间分布可以代表omega的时间，但是re-infection会低估一些omega

FilePath <- File.CreateMainFolder(
  path = "Output/",
  FolderName = "ParmInference",
  Date = NULL
)


#  SimData <-
Model.RunSim.ode(
  Parm = Parameter.Create(
    # Seasonal force
    # years = 4,
    # year_start = "1999-01-01",
    # year_end = "2022-12-31",
    penalty = 0,

    # R0 = c(1.41, 1.07, 1.7, 1.88),

    # Virus competition
    comp = c(0.8, 0.3, -0.97, -0.13),
    #c(6.7, 3.4, -0.3, -9.8),
    #c(1.68, 0.45, -1.19, 0.94),
    # c(56, 28, 23, 0), # c(2.6, 0.8, 0.8, 0.5),
    # c(21,16,4,0), # c(3.05, 2.82, 1.49, -2.74), # c(1, 1, 1, 1),

    # # NPI susceptibility
    # NPISes = c(1, 1, 1, 1),

    # # Base immunity
    # beta_seasonal = Test1[nrow(Test1), 9],
    # NPI_start = c("2020-03-10", "2020-10-10", "2021-01-10", "2021-05-10"), # not use this param any more
    # NPI_end = c("2020-06-10", "2020-12-10", "2021-01-10", "2021-05-10"), # not use this param any more
    # NPI_value = c(0.8, 0.8, 0.8, 0.8), # range[0,1] 1 means totally no transmission(Full NPI)
    # decay_coef = c(0.01, 0.01, 0.01, 0.01), # range[0,Inf] 0 means no decay
    NPI = FALSE
  ),
  after = as.Date("2015-08-31")
)

Model.RunSim.LLH(
  Parm = Parameter.Create(
    NPI = FALSE
  ),
  TargetDat = RealData_c_Long[Location == "CANADA"],
  Method = "BetaBinomial"
  # "Dirichlet" # "LogDiff" # "RatioAbs" #
)

system.time(MCMC.MH(
  Initial = c(0.2, 5, 1, 0.8),
  n_iterations = 50,
  step = 0.1,
  TargetDat = RealData_c_Long[Location == "CANADA"],
  Method = "BetaBinomial"
  # "Dirichlet" # "LogDiff" # "RatioAbs" #
))

source("Code/MCMC.R")
source("Code/Untitled-2.R")
# Beta binomial likelihood
fit_bb <- MCMC.MH(
  Initial = list(comp = rep(1, 4), rho = 0.02),
  n_iterations = 5000,
  TargetDat = RealData_c_Long[Location == "CANADA"],
  step = list(comp = 0.08, z_rho = 0.12),
  Method = "BetaBinomial",
  BaseParm = Parameter.Create(),
  LLHArg = list(proxy = "Inc_sum"),
  PriorArg = list(sd_comp = 1),
  verbose = TRUE
)
post_bb <- MCMC.DecodeChain(fit_bb)
MCMC.TracePlot(post_bb)
MCMC.PostProcess(post_bb, burn_in = 1000, thin = 10)

# Dirichlet likelihood
fit_dir <- MCMC.MH(
  Initial = list(comp = rep(1, 4), kappa = 100),
  n_iterations = 5000,
  step = list(comp = 0.08, log_kappa = 0.10),
  Prep = Prep,
  Method = "Dirichlet",
  BaseParm = Parameter.Create(),
  LLHArg = list(
    proxy = "Inc_sum",
    weight = "sum_y_capN",
    c = 1e-8
  ),
  PriorArg = list(sd_comp = 1),
  verbose = TRUE
)

post_dir <- MCMC.DecodeChain(fit_dir)
MCMC.TracePlot(post_dir)
MCMC.PostProcess(post_dir, burn_in = 1000, thin = 10)


# Ratio-absolute difference likelihood
fit_ra <- MCMC.MH(
  Initial = list(comp = rep(1, 4), sigma_ratio = 0.7, sigma_abs = 0.7),
  n_iterations = 5000,
  step = list(comp = 0.08, log_sigma_ratio = 0.08, log_sigma_abs = 0.08),
  Prep = Prep,
  Method = "RatioAbs",
  BaseParm = Parameter.Create(),
  LLHArg = list(proxy = "Inc_sum", pseudo = 0.5),
  verbose = TRUE
)
post_ra <- MCMC.DecodeChain(fit_ra)
MCMC.TracePlot(post_ra)
MCMC.PostProcess(post_ra, burn_in = 1000, thin = 10)


# Log-difference likelihood
fit_ld <- MCMC.MH(
  Initial = list(comp = rep(1, 4), sigma = 0.7),
  n_iterations = 5000,
  step = list(comp = 0.08, log_sigma = 0.08),
  Prep = Prep,
  Method = "LogDiff",
  BaseParm = Parameter.Create(),
  LLHArg = list(proxy = "Inc_sum", pseudo = 0.5, transform = "logit"),
  verbose = TRUE
)
post_ld <- MCMC.DecodeChain(fit_ld)
MCMC.TracePlot(post_ld)
MCMC.PostProcess(post_ld, burn_in = 1000, thin = 10)


bb <- MCMC.MH(
  Initial = list(comp = rep(0.5, 4), rho = 0.02),
  n_iterations = 5000,
  TargetDat = RealData_c_Long[Location == "CANADA"],
  step = list(comp = 0.08, z_rho = 0.12),
  Method = "RatioAbs",
  BaseParm = Parameter.Create(),
  LLHArg = list(proxy = "Inc_sum"),
  PriorArg = list(sd_comp = 1),
  verbose = TRUE
)
MCMC.DecodeChain(bb)

cc <- MCMC.MH(
  Initial = list(comp = rep(10, 4), rho = 0.02),
  n_iterations = 5000,
  TargetDat = RealData_c_Long[Location == "CANADA"],
  step = list(comp = 0.08, z_rho = 0.12),
  Method = "BetaBinomial",
  BaseParm = Parameter.Create(),
  LLHArg = list(proxy = "Inc_sum"),
  PriorArg = list(sd_comp = 1),
  verbose = TRUE
)
a <- MCMC.DecodeChain(cc)

Prep <- Inference.Setup(
  TargetDat = RealData_c_Long[Location == "Manitoba"],
  # "Atlantic" # "British Columbia" # "CANADA" # "Manitoba" # "Prairies" # "Province of Ontario"
  # "Province of Quebec" # "Saskatchewan"
  BaseParm = Parameter.Create(),
  after = as.Date("2015-08-31")
)
chain1 <- MCMC.MH.Adaptive(
  Initial = c(1, -1, -2, 2),
  n_iterations = 5,
  step = 0.15,
  # Prep = Prep,
  TargetDat = RealData_c_Long_Dir[Location == "CANADA"],
  Method = "Dirichlet",
  BaseParm = Parameter.Create(),
  LLHArg = list(
    proxy = "Inc_sum",
    weight = "sum_y_capN",
    c = 1e-8
  ),
  PriorArg = list(sd_comp = 1),
  Adapt = list(
    use = TRUE,
    start = 500,
    end = 5000,
    every = 25,
    target_accept = 0.234
  ),
  verbose = TRUE,
  InferArg = list(
    infer_penal = FALSE,
    infer_NPISes = FALSE
  )
)
dd <- MCMC.DecodeChain(chain1)
MCMC.TracePlot(dd)


dd <- MCMC.RunChains.Future(
  n_chain = 4,
  n_iterations = 5000,
  step = 0.15,
  Prep = Prep,
  Method = "BetaBinomial",
  BaseParm = Parameter.Create(),
  LLHArg = list(epsilon = 1e-9),
  PriorArg = list(sd_comp = 1),
  Adapt = list(start = 200, end = 5000, every = 25),
  workers = 4,
  seed = 20250312,
  source_r = normalizePath("Code/Untitled-2.R"),
  source_cpp = normalizePath("Code/ParmInferenceCpp_V6.cpp"),
  wait = FALSE
)
MCMC.Future.Ready(dd)
dd <- MCMC.Future.Collect(dd)


a <- Model.RunSim(
  Parm = Parameter.Create(
    dt = 3.5,
    years = 2,
    comp_IFVA = 0,
    comp_IFVB = 0,
    comp_RSV = 0,
    comp_RV = 0,
    beta_seasonal = 0.2,
    phi = 330,
    R0_IFVA = 1.5,
    R0_IFVB = 2,
    R0_RSV = 2.5,
    R0_RV = 5,
    Penal = 1,
    base_immune = 0.8,
  ),
  ncores = 10,
  NPI = FALSE,
  BaseImmu = TRUE
) # , seeds = 971889
a$fig1
a$fig2
FindPeak(
  a[[1]],
  StartTime = 290,
  NPI_start = 624,
  NPI_end = 676,
  Threshold = 0.4,
  span = 25,
  Offset = 25,
  Method = "Likelihood",
  TargetDat = Target
)


Model.RunSim.Multi(
  Parm = Parameter.Create(
    dt = 3.5,
    years = 2,
    comp_IFVA = 0,
    comp_IFVB = 0,
    comp_RSV = 0,
    comp_RV = 0,
    beta_seasonal = 0.2,
    phi = 330,
    R0_IFVA = 1.5,
    R0_IFVB = 2,
    R0_RSV = 2.5,
    R0_RV = 5,
    Penal = 1,
    base_immune = 0.8,
  ),
  n_simulations = 3,
  ncores = 10,
  NPI = FALSE,
  BaseImmu = TRUE
)


# 画多次模拟的结果
Plot.SimResult(results)
a1 <- c(1, 2, 3, 1, 4, 2, 1, 3)
a2 <- c(2, 5, 6, 5, 1, 2, 3, 5)


#########
SimDat <- copy(SimData$Data)
# SimDat <- SimDat[time >= as.Date("2016-09-01")]

# Find Max for each virus in 12 months, and calculate the proportion
SimDat[, Monday := Time - (as.integer(strftime(Time, "%u")) - 1L)]
SimDat <- SimDat[,
  lapply(.SD, sum, na.rm = TRUE),
  .SDcols = !c("Time", "ISOweek", "Monday"),
  by = c("Monday", "ISOweek")
]
SimDat[,
  paste0(c("Virus_1", "Virus_2", "Virus_3", "Virus_4"), "_roll") := lapply(
    .SD,
    function(x) frollmax(x, n = 52, align = "center", na.rm = TRUE)
  ),
  .SDcols = c("Virus_1", "Virus_2", "Virus_3", "Virus_4")
]


SimDat[,
  ":="(
    Virus_1 = Virus_1 / Virus_1_roll,
    Virus_2 = Virus_2 / Virus_2_roll,
    Virus_3 = Virus_3 / Virus_3_roll,
    Virus_4 = Virus_4 / Virus_4_roll
  )
]
SimDat[,
  c(
    "Virus_1_roll",
    "Virus_2_roll",
    "Virus_3_roll",
    "Virus_4_roll"
  ) := NULL
]
SimDat <- melt(
  SimDat,
  id.vars = c("Monday", "ISOweek"),
  measure.vars = patterns("Virus_"),
  variable.name = "Virus",
  value.name = "Prop"
)

# Calculate the adjusted observed cases based on the proportion
ObsDat <- copy(RealData_c_Long_Dir[Location == "CANADA"])
setorder(ObsDat, Virus, Monday)
ObsDat[,
  y_max_52 := frollmax(y, n = 52, align = "center", na.rm = TRUE),
  by = Virus
]

MergeDat <- merge(
  SimDat,
  ObsDat,
  by = c("ISOweek", "Monday", "Virus"),
  all.x = TRUE
)
MergeDat[, Virus_est := Prop * y_max_52][,
  llh := dpois(round(y), Virus_est, log = TRUE)
]
# Calsulate the log likelihood
LLH <- sum(MergeDat$llh, na.rm = TRUE)

MCMC.DecodeChain(resA[[1]], include_eta = FALSE)
a <- MCMC.PostProcess(
  resA,
  burn_in = 5000,
  thin = 10,
  ModelParm = list(penalty = 0, NPI = FALSE),
  ObsDat = RealData_c_Dir[Location == "CANADA"],
)
