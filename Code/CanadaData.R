### Canada data prepare ####
#################################
# Canadian data was split into two segments: pre-COVID surveillance (through November 19, 2022)
# and post-COVID surveillance (from November 20, 2022 onward), loaded separately.
#################################

### -----------------------------------------------------------------------
## 2015-2022
CanadaFile <- list.files(paste0(getwd(), "/Data/CanadaP1/"))
CanadaData <- rbindlist(
  lapply(CanadaFile, function(x) {
    dt <- fread(
      paste0(getwd(), "/Data/CanadaP1/", x),
      col.names = c(
        "Location",
        "IV_Tested",
        "A(H1)pdm09",
        "A(H3)",
        "A(UnS)",
        "IAV",
        "IBV",
        "RSV_Tested",
        "RSV",
        "PIV_Tested",
        "PIV_1",
        "PIV_2",
        "PIV_3",
        "PIV_4",
        "Other_PIV",
        "ADV_Tested",
        "ADV",
        "MPV_Tested",
        "MPV",
        "RV_Tested",
        "RV",
        "CoV_Tested",
        "CoV"
      )
    )
    date <- gsub("\\.csv$", "", basename(x))
    dt[, Date := as.Date(date)]
    return(dt)
  }),
  fill = TRUE
)

### -----------------------------------------------------------------------
## 2022-2024
CanadaFile2 <- list.files(paste0(getwd(), "/Data/CanadaP2/"))
CanadaData2 <- rbindlist(
  lapply(CanadaFile2, function(x) {
    dt <- fread(
      paste0(getwd(), "/Data/CanadaP2/", x),
      drop = c("SARS-CoV-2 Tested", "SARS-CoV-2 Positive"),
      col.names = c(
        "Location",
        "IV_Tested",
        "A(H1)pdm09",
        "A(H3)",
        "A(UnS)",
        "IAV",
        "IBV",
        "RSV_Tested",
        "RSV",
        "PIV_Tested",
        "PIV_1",
        "PIV_2",
        "PIV_3",
        "PIV_4",
        "Other_PIV",
        "ADV_Tested",
        "ADV",
        "MPV_Tested",
        "MPV",
        "RV_Tested",
        "RV",
        "CoV_Tested",
        "CoV"
      )
    )
    date <- gsub("\\.csv$", "", basename(x))
    dt[, Date := as.Date(date)]
    return(dt)
  })
)
CanadaData_All <- rbind(CanadaData, CanadaData2)

### tidy up the data
Dat <- copy(CanadaData_All)
Dat <- Dat[order(Location, Date)]
Dat <- setcolorder(Dat, c("Location", "Date"))
Dat[, Date := as.Date(Date)]

# fwrite(Dat, "Data/CanadaData.csv")

### -----------------------------------------------------------------------
## COVID-19 data
COVIDDat <- fread("Data/covid19-download.csv")
COVIDDat <- COVIDDat[
  date <= as.Date("2024/06/08") & prname != "Repatriated Travellers",
  c("prname", "date", "numtotal_last7", "ratecases_last7")
]
# numtotal_last7: Number of cases reported in the reporting week
# ratecases_last7: Case rate in the reporting week per 100,000 population
colnames(COVIDDat) <- c("Location", "Date", "COVID", "COVID_IR")

COVIDDat[,
  c("COVID", "COVID_IR") := lapply(.SD, function(x) {
    as.numeric(replace(x, x == "-", NA))
  }),
  .SDcols = c("COVID", "COVID_IR")
]

NameMapCOVID <- c(
  "Ontario" = "Province of Ontario",
  "Quebec" = "Province of Québec",
  "Newfoundland and Labrador" = "Newfoundland",
  "Canada" = "CANADA"
)

# 使用 fcoalesce 实现：如果在映射表中找到就用新名，找不到(NA)就用原名
COVIDDat[, Location := fcoalesce(NameMapCOVID[Location], Location)]

### -----------------------------------------------------------------------
## Merge all data
MergeDat <- merge(
  Dat,
  COVIDDat,
  by = c("Location", "Date"),
  all.x = TRUE
)


NameMapAll <- c(
  "Centre-du-Québec" = "Central Quebec",
  "Montréal-Laval" = "Montreal-Laval",
  "Montérégie" = "Monteregie",
  "Ouest du Québec" = "West of Quebec",
  "P.H.O.L. - Sault Ste. Marie" = "P.H.O.L. - Sault Sainte Marie",
  "Province of Québec" = "Province of Quebec",
  "Québec-Chaudière-Appalaches" = "Quebec-Chaudiere-Appalaches",
  "Région Nord-Est" = "North-East Region"
)

MergeDat[, Location := fcoalesce(NameMapAll[Location], Location)]
fwrite(MergeDat, "Data/CanadaData.csv")


### -----------------------------------------------------------------------
## Filter high quality data
HQDat <- copy(MergeDat)
HQDat <- HQDat[
  Location %in%
    c(
      "Atlantic",
      "Province of Quebec",
      "Province of Ontario",
      "Manitoba",
      "Saskatchewan",
      "Prairies",
      "British Columbia",
      "CANADA"
      # Alberta : missing data for RSV
    ),
  c(
    "Location",
    "Date",
    "IV_Tested",
    "IAV",
    "IBV",
    "RSV_Tested",
    "RSV",
    "RV_Tested",
    "RV",
    "COVID",
    "COVID_IR"
  )
]
HQDat[, IAVP := (IAV / IV_Tested)][, IBVP := (IBV / IV_Tested)][,
  RSVP := (RSV / RSV_Tested)
][, RVP := (RV / RV_Tested)]

LocaDateList <- expand.grid(
  Location = as.character(unique(HQDat$Location)),
  Date = seq(as.Date(min(HQDat$Date)), as.Date(max(HQDat$Date)), by = "week")
) |>
  setDT()

HQDat <- merge(HQDat, LocaDateList, by = c("Location", "Date"), all = TRUE)
HQDat[, Location := as.character(Location)]
fwrite(HQDat, "Data/HQData.csv")
save(HQDat, file = "Data/HQDat.RData")

########## Visualization
# a <- fread("Data/CanadaData.csv")
# a <- a[, Date := as.Date(Date)]
# unique(a$Location)

# aa <- table(a$Date, a$Location)
# View(aa)

# a %>%
#   mutate(reported = 1) %>%
#   complete(Date, Location, fill = list(reported = 0)) %>%
#   ggplot(aes(x = Date, y = Location, fill = factor(reported))) +
#   geom_tile()

# a %>%
#   ggplot(aes(x = Date, y = IV_Tested, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = IAV, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = IBV, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = RSV_Tested, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = RSV, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = RV_Tested, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   ggplot(aes(x = Date, y = RV, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   mutate(IAVP = IAV / IV_Tested * 100) %>%
#   ggplot(aes(x = Date, y = IAVP, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   mutate(IBVP = IBV / IV_Tested * 100) %>%
#   ggplot(aes(x = Date, y = IBVP, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   mutate(RSVP = RSV / RSV_Tested * 100) %>%
#   ggplot(aes(x = Date, y = RSVP, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")

# a %>%
#   mutate(RVP = RV / RV_Tested * 100) %>%
#   ggplot(aes(x = Date, y = RVP, )) +
#   geom_line() +
#   facet_wrap(~Location, scales = "free_y")
