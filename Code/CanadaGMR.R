### Google Mobility Report data
GMR <- fread("Data/Global_Mobility_Report.csv")
GMR <- GMR[country_region == "Canada" & sub_region_2 == "" & metro_area == "", ]
GMR[sub_region_1 == "", "sub_region_1"] <- "Canada"
GMR <- GMR[, c(
    "date",
    "sub_region_1",
    "retail_and_recreation_percent_change_from_baseline",
    "grocery_and_pharmacy_percent_change_from_baseline",
    "parks_percent_change_from_baseline",
    "transit_stations_percent_change_from_baseline",
    "workplaces_percent_change_from_baseline",
    "residential_percent_change_from_baseline"
)]
colnames(GMR) <- c(
    "Date",
    "Location",
    "retail",
    "grocery",
    "parks",
    "transit",
    "workplace",
    "residential"
)

GMR <- expand.grid(
    Location = unique(GMR$Location),
    Date = seq(as.Date("2020/02/15"), as.Date("2022/10/15"), by = "day")
) |>
    setDT() |>
    merge(GMR, by = c("Location", "Date"), all = TRUE)
fwrite(GMR, "Data/GMR_Canada.csv")


GMR <- GMR[, .(Location, Date, retail)][,
    retail := (as.numeric(retail) + 100) / 100
][,
    retail := frollmean(retail, 7, fill = NA),
    by = Location
]
GMR <- na.omit(GMR)
# save(GMR, file = "Data/GMR_Canada.RData")

HQGMR <- GMR[
    Location %in%
        c(
            "British Columbia",
            "Canada",
            "Manitoba",
            "New Brunswick",
            "Ontario",
            "Quebec",
            "Saskatchewan"
        ),
]
# 使用"New Brunswick"替代"Atlantic",

NameMapAll <- c(
    "Canada" = "CANADA",
    "New Brunswick" = "Atlantic",
    "Ontario" = "Province of Ontario",
    "Quebec" = "Province of Quebec"
)

HQGMR[, Location := fcoalesce(NameMapAll[Location], Location)]
save(HQGMR, file = "Data/HQGMR.RData")

# HQGMR |>
#     as.data.frame() |>
#     ggplot(aes(x = Date, y = retail, group = Location)) +
#     geom_line() +
#     facet_wrap(~Location)
