library(dplyr)
library(ggplot2)
library(sf)
# library(future)
# plan(multisession)

dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

cat12_only <- TRUE
# cat12_only <- FALSE

source("analysis/functions.R")

x <- sf::read_sf("data-raw/mpatt_survey_overlaps.gdb", type = 7)
table(x$Category_Simple)
table(x$Common_site_name_Site_Profile)
table(x$Desig_Tool_Detail)
table(x$SurveyOverlap)

f <- list.files("~/src/gfsynopsis-2021/report/data-cache-feb-2023/",
  full.names = TRUE
)
s <- readRDS(f[grep("arrow", f)])$survey_sets
s <- filter(s,survey_abbrev %in% c("SYN QCS", "SYN HS", "SYN WCHG"))

trawl_removed <- dplyr::filter(x, SurveyOverlap %in% c("Category 1", "Category 2", "Gwaii Haanas Site- Strict Protection"))

g <- ggplot(trawl_removed, aes(colour = SurveyOverlap)) + geom_sf() +
  theme_light() +
  theme(legend.position = "top")
  # geom_point(data = s, mapping = aes())
ggsave("figs/rough-map-trawl.png", width = 7, height = 7)

# filter(x, SurveyOverlap %in% "Gwaii Haanas Site- Strict Protection") |>
#   ggplot() + geom_sf(mapping = aes(colour = SurveyOverlap))
#
# filter(x, SurveyOverlap %in% "Category 1") |>
#   ggplot() + geom_sf(mapping = aes(colour = SurveyOverlap))
#
# filter(x, SurveyOverlap %in% "Category 2") |>
#   ggplot() + geom_sf(mapping = aes(colour = SurveyOverlap))
#
# g <- filter(x, SurveyOverlap %in% c("Category 1", "Category 2")) |>
#   ggplot() + geom_sf(mapping = aes(colour = SurveyOverlap))
# ggsave("~/Desktop/map.png", width = 8, height = 8)
#
# filter(x, SurveyOverlap %in% "Category 2") |>
#   ggplot() + geom_sf(mapping = aes(colour = SurveyOverlap))
#
#
# x <- sf::read_sf("data-raw/MPATT_P2_Nov25_limited_attributes.gdb/", type = 7)
# table(x$Category_Detailed)
# # Category 1
# # 171
# # Category 2
# # 12
# # Existing MPA/RCA - 'as-is, where-is'
# # 54
# # Existing MPA/RCA - 'as-is, where-is' *
# #   24
# # Existing MPA/RCA - 'as-is, where-is' **
# #   29
# # Existing MPA/RCA - 'as-is, where-is' ***
# #   67
# table(x$Prop_Desig_Tool)
# # BC WMA Fisheries Act + BC tools
# # 2                        7                        9
# # Fisheries Act tool                     MNWA                    NMCAR
# # 17                       51                       76
# # Oceans Act MPA                      TBD
# # 18                       12
# table(x$SUBREGION)
# # CC  HG  NC NVI
# # 83  84  68 122

# if (cat12_only) {
#   trawl_removed <- dplyr::filter(x, SurveyOverlap %in% c("Category 1", "Category 2", "Gwaii Haanas Site- Strict Protection", "Gwaii Haanas Site- Multi Use Gwaii"))
# } else {
#   stop("Not implemented")
#   # .cat <- c("Category 1", "Category 2", "Existing MPA/RCA - 'as-is, where-is'",
#   #   "Existing MPA/RCA - 'as-is, where-is' *", "Existing MPA/RCA - 'as-is, where-is' **",
#   #   "Existing MPA/RCA - 'as-is, where-is' ***")
#   # trawl_removed <- dplyr::filter(x, Category_Detailed %in% .cat)
# }

# x |>
#   dplyr::filter(Category_Detailed == "Category 1") |>
#   plot()


# x <- sf::read_sf("data-raw/Spatial_N1_May17_2021.gdb")
# # plot(x["hu_co_demersalfishing_bottomtrawling_d"])
# trawl_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("X"))
# plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

# https://gis.stackexchange.com/questions/389814/r-st-centroid-geos-error-unknown-wkb-type-12

# library(sf)
# library(gdalUtilities)
ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  sf::st_write(X, tmp1)
  gdalUtilities::ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- sf::st_read(tmp2)
  sf::st_sf(sf::st_drop_geometry(X), geom = sf::st_geometry(Y))
}
trawl_removed <- ensure_multipolygons(trawl_removed)

saveRDS(trawl_removed, file = "data-generated/Cat1_2_GH_Apr2023.rds")

assign_restricted_tows <- function(trawl_dat) {
  orig <- trawl_dat
  trawl_dat <- trawl_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(sf::st_crs(trawl_removed))
  intersected <- sf::st_intersects(trawl_dat, trawl_removed)
  remain <- which(lengths(intersected) == 0)
  lost <- which(lengths(intersected) > 0)
  trawl_dat$restricted <- lengths(intersected) > 0
  trawl_dat <- as.data.frame(trawl_dat) %>% select(-geometry)
  trawl_dat$longitude <- orig$longitude
  trawl_dat$latitude <- orig$latitude
  utm_zone9 <- 3156
  coords <- orig %>%
    sf::st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    sf::st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  trawl_dat$X <- coords$X / 1000
  trawl_dat$Y <- coords$Y / 1000
  trawl_dat %>% as_tibble()
}

f <- list.files("~/src/gfsynopsis-2021/report/data-cache-feb-2023/",
  full.names = TRUE
)
aleut <- readRDS(f[grep("arrow", f)])$survey_sets
aleut <- filter(aleut, year %in% c(2021, 2022))
# aleut <- filter(aleut, survey_abbrev == "SYN WCHG", year == 2015)
# aleut <- filter(aleut, survey_abbrev == "SYN WCHG", year == 2015)
aleut <- assign_restricted_tows(aleut)
ggplot(aleut, aes(longitude, latitude, size = density_kgpm2, colour = restricted)) + geom_point()

if (Sys.info()[["user"]] == "seananderson") {
  f <- list.files("~/src/gfsynopsis-2021/report/data-cache-feb-2023/",
    full.names = TRUE
  )
  f <- f[!grepl("cpue", f)]
  f <- f[!grepl("iphc", f)]
  synoptic_data <- purrr::map_dfr(seq_along(f), function(i) {
    cat(f[i], "\n")
    d <- readRDS(f[i])$survey_sets
    filter(d, survey_abbrev %in% c("SYN QCS", "SYN HS", "SYN WCHG")) %>%
      select(
        year, survey_abbrev, species_science_name, species_common_name, speed_mpm, duration_min,
        density_kgpm2, catch_weight, doorspread_m, tow_length_m, latitude, longitude, grouping_code, area_km2, depth_m
      )
  })
  saveRDS(synoptic_data, file = "data-raw/syn-survey-data.rds")
} else {
  synoptic_data <- readRDS("data-raw/syn-survey-data.rds")
}

dat <- assign_restricted_tows(synoptic_data)

frac_pos_df <- dat %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(
    frac_pos = sum(density_kgpm2 > 0) / length(density_kgpm2), .groups = "drop"
  )

tokeep <- frac_pos_df %>% filter(frac_pos > 0.05)
notkeep <- frac_pos_df %>% filter(frac_pos <= 0.05)

nrow(frac_pos_df)
nrow(tokeep)

dat_to_fit <- left_join(tokeep, dat, by = c("species_common_name", "survey_abbrev"))

dat_to_fit <- dat_to_fit |> filter(!(year == 2014 & survey_abbrev == "SYN WCHG"))

saveRDS(dat_to_fit, file = "data-generated/dat_to_fit.rds")

# HBLL -----------------

# x <- sf::read_sf("data-raw/Spatial_N1_May17_2021.gdb")
x <- sf::read_sf("data-raw/mpatt_survey_overlaps.gdb", type = 7, layer = "MPATT_Q1_full_march2023_Cat1_Cat2_GH_singlepart")
unique(x$SurveyOverlap)

g <- ggplot(x, aes(colour = SurveyOverlap)) + geom_sf() +
  theme_light() +
  theme(legend.position = "right")
# geom_point(data = s, mapping = aes())
ggsave("figs/rough-map-rca.png", width = 10, height = 7)

hbll_removed <- dplyr::filter(x, SurveyOverlap %in% c("Category 1", "Category 2", "Gwaii Haanas Site- Strict Protection - RCA overlap"))

g <- ggplot(hbll_removed, aes(colour = SurveyOverlap)) + geom_sf() +
  theme_light() +
  theme(legend.position = "top")
ggsave("figs/rough-map-hbll.png", width = 7, height = 7)

hbll_removed <- ensure_multipolygons(hbll_removed)
ll_removed <- hbll_removed

saveRDS(ll_removed, file = "data-generated/LL_Cat1_2_GH_Apr2023.rds")

# ll_removed <- trawl_removed
# plot(x["hu_co_demersalfishing_bottomlongline_d"])
# ll_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomlongline_d %in% c("X"))
# ll_removed <- sf::st_cast(ll_removed, "MULTIPOLYGON")
# saveRDS(ll_removed, file = "data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
# png("figs/ll-map.png", width = 5, height = 6, units = "in", res = 200)
# plot(ll_removed["hu_co_demersalfishing_bottomlongline_d"])
# dev.off()

assign_restricted_tows_hbll <- function(dat) {
  orig <- dat
  dat <- dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(sf::st_crs(ll_removed))
  intersected <- sf::st_intersects(dat, ll_removed)
  remain <- which(lengths(intersected) == 0)
  lost <- which(lengths(intersected) > 0)
  dat$restricted <- lengths(intersected) > 0
  dat <- as.data.frame(dat) %>% select(-geometry)
  dat$longitude <- orig$longitude
  dat$latitude <- orig$latitude
  utm_zone9 <- 3156
  coords <- orig %>%
    sf::st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    sf::st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  dat$X <- coords$X / 1000
  dat$Y <- coords$Y / 1000
  dat %>% as_tibble()
}

# aleut <- readRDS(f[grep("aleut", f)])$survey_sets
# aleut <- filter(aleut, survey_abbrev == "HBLL OUT N", year == 2015)
# aleut <- assign_restricted_tows_hbll(aleut)
# ggplot(aleut, aes(longitude, latitude, size = density_ppkm2, colour = restricted)) + geom_point()

if (Sys.info()[["user"]] == "seananderson") {
  f <- list.files("~/src/gfsynopsis-2021/report/data-cache-feb-2023/",
    full.names = TRUE
  )
  f <- f[!grepl("cpue", f)]
  f <- f[!grepl("iphc", f)]
  hbll_data <- purrr::map_dfr(seq_along(f), function(i) {
    cat(f[i], "\n")
    d <- readRDS(f[i])$survey_sets
    filter(d, survey_abbrev %in% c("HBLL OUT N")) %>%
      select(
        year, fishing_event_id, survey_abbrev, species_science_name, species_common_name,
        density_ppkm2, latitude, longitude, grouping_code, area_km2, depth_m,
        hook_count, catch_count
      )
  })
  saveRDS(hbll_data, file = "data-raw/hbll-survey-data.rds")
} else {
  if (Sys.info()[["user"]] == "dfomac") {
    f <- list.files("~/github/dfo/gfsynopsis/report/data-cache/",
      full.names = TRUE
    )
    f <- f[!grepl("cpue", f)]
    f <- f[!grepl("iphc", f)]
    hbll_data <- purrr::map_dfr(seq_along(f), function(i) {
      cat(f[i], "\n")
      d <- readRDS(f[i])$survey_sets
      filter(d, survey_abbrev %in% c("HBLL OUT N")) %>%
        select(
          year, fishing_event_id, survey_abbrev, species_science_name, species_common_name,
          density_ppkm2, latitude, longitude, grouping_code, area_km2, depth_m,
          hook_count, catch_count
        )
    })
    saveRDS(hbll_data, file = "data-raw/hbll-survey-data.rds")
  } else {
    hbll_data <- readRDS("data-raw/hbll-survey-data.rds")
  }
}

dat <- assign_restricted_tows_hbll(hbll_data)

frac_pos_df <- dat %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(frac_pos = sum(density_ppkm2 > 0) / length(density_ppkm2), .groups = "drop")

frac_pos_df %>%
  ungroup() %>%
  arrange(-frac_pos) %>%
  top_n(30) %>%
  as.data.frame()

tokeep <- frac_pos_df %>% filter(frac_pos > 0.05)
notkeep <- frac_pos_df %>% filter(frac_pos <= 0.05)

nrow(frac_pos_df)
nrow(tokeep)

dat_to_fit <- left_join(tokeep, dat, by = c("species_common_name", "survey_abbrev"))

saveRDS(dat_to_fit, file = "data-generated/dat_to_fit_hbll.rds")

dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
hbll_grid <- gfplot::hbll_n_grid$grid
utm_zone9 <- 3156
coords <- hbll_grid %>%
  sf::st_as_sf(crs = 4326, coords = c("X", "Y")) %>%
  sf::st_transform(utm_zone9)
coords_restr <- shrink_a_survey(coords, ll_removed, plot = TRUE)
coords <- coords %>%
  sf::st_coordinates() %>%
  as.data.frame()
coords$X <- coords$X / 1000
coords$Y <- coords$Y / 1000
coords$restricted <- coords_restr$restricted
coords$depth <- hbll_grid$depth
grid <- coords %>%
  expand_prediction_grid(years = sort(unique(dat_to_fit$year))) %>%
  as_tibble() %>%
  arrange(year, X, Y, depth)

saveRDS(grid, "data-generated/hbll-n-grid-w-restr.rds")

dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
syn_grid <- gfplot::synoptic_grid
syn_grid$X <- syn_grid$X * 1000
syn_grid$Y <- syn_grid$Y * 1000
utm_zone9 <- 3156
coords <- syn_grid %>%
  sf::st_as_sf(crs = utm_zone9, coords = c("X", "Y"))
coords_restr <- shrink_a_survey(coords, trawl_removed, plot = TRUE)
coords <- coords %>%
  sf::st_coordinates() %>%
  as.data.frame()
coords$X <- coords$X / 1000
coords$Y <- coords$Y / 1000
coords$restricted <- coords_restr$restricted
coords$survey_abbrev <- syn_grid$survey
grid <- coords %>%
  expand_prediction_grid(years = sort(unique(dat_to_fit$year))) %>%
  as_tibble() %>%
  arrange(survey_abbrev, year, X, Y)

saveRDS(grid, "data-generated/syn-grid-w-restr.rds")
