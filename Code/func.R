CreateMainFolder <- function(path){
  CurrentDate <- format(Sys.time(), "%Y%m%d")
  if (!file.exists(paste0(path, CurrentDate))) {
    dir.create(paste0(path, CurrentDate))
  }
  return(paste0(path, CurrentDate, "/"))
}

CreateSubFolder <- function(path, folder) {
  if (!file.exists(paste0(path, folder))) {
    dir.create(paste0(path, folder))
  }
  return(paste0(path, folder, "/"))
}
