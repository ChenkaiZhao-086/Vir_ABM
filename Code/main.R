# rm(list=ls())

packages <- c("readxl", "tidyverse", "data.table", "foreach", "doParallel")

for (i in packages) {
  suppressPackageStartupMessages(library(i, character.only = TRUE, quietly = TRUE))
}
# lapply(c("readxl", "tidyverse", "data.table","foreach", "doParallel"),
#        library,character.only=T)

source("Code/func.R")
FilePath <- File.CreateMainFolder(path = "Output/") # "Output/20230816_Final/"

set.seed(971889)

Parameter <- Parameter.Create()
