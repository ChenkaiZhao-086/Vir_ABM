packages <- c("readxl", "tidyverse", "data.table","foreach", "doParallel")

for (i in packages) {
  suppressPackageStartupMessages(library(i, character.only = T, quietly = T))
}

source("Code/func.R")
FilePath <- CreateMainFolder(path = "Output/") # "Output/20230816_Final/"

set.seed(971889)

# cl <- makeCluster(12)
# registerDoParallel(cl)