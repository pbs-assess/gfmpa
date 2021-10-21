library(dplyr)
library(sf)
library(future)
plan(multisession)

dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

source("analysis/functions.R")

x <- sf::read_sf("data-raw/Spatial_N1_May17_2021.gdb")
# plot(x["hu_co_demersalfishing_bottomtrawling_d"])
trawl_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("X"))
saveRDS(trawl_removed, file = "data-generated/demersalfishing_bottomtrawling_X.rds")
# plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

# aleut <- readRDS(f[grep("aleut", f)])$survey_sets
# aleut <- filter(aleut, survey_abbrev == "SYN HS", year == 2015)
# aleut <- assign_restricted_tows(aleut)
# ggplot(aleut, aes(longitude, latitude, size = density_kgpm2, colour = restricted)) + geom_point()

assign_restricted_tows <- function(trawl_dat) {
  orig <- trawl_dat
  trawl_dat <- trawl_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(sf::st_crs(trawl_removed))
  intersected <- sf::st_intersects(trawl_dat, trawl_removed)
  remain <- which(lengths(intersected) == 0)
  lost <- which(lengths(intersected) > 0)
  # remain_df <- pcod[remain,]
  # lost_df <- pcod[lost,]
  # plot(pcod$X, pcod$Y, col = "black")
  # points(lost_df$X, lost_df$Y, col = "red")
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

if (Sys.info()[['user']] == "seananderson") {
  f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
    full.names = TRUE)
  synoptic_data <- furrr::future_map_dfr(seq_along(f), function(i) {
    d <- readRDS(f[i])$survey_sets
    filter(d, survey_abbrev %in% c("SYN QCS", "SYN HS", "SYN WCHG")) %>%
      select(year, survey_abbrev, species_science_name, species_common_name,
        density_kgpm2, latitude, longitude, grouping_code, area_km2, depth_m)
  })
  saveRDS(synoptic_data, file = "data-raw/syn-survey-data.rds")
} else {
  synoptic_data <- readRDS("data-raw/syn-survey-data.rds")
}

dat <- assign_restricted_tows(synoptic_data)

frac_pos_df <- dat %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(frac_pos = sum(density_kgpm2 > 0) / length(density_kgpm2), .groups = "drop")

tokeep <- frac_pos_df %>% filter(frac_pos > 0.05)
notkeep <- frac_pos_df %>% filter(frac_pos <= 0.05)

nrow(frac_pos_df)
nrow(tokeep)

dat_to_fit <- left_join(tokeep, dat, by = c("species_common_name", "survey_abbrev"))

saveRDS(dat_to_fit, file = "data-generated/dat_to_fit.rds")

# HBLL -----------------

x <- sf::read_sf("data-raw/Spatial_N1_May17_2021.gdb")
# plot(x["hu_co_demersalfishing_bottomlongline_d"])
ll_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomlongline_d %in% c("X"))
ll_removed <- sf::st_cast(ll_removed, "MULTIPOLYGON")
saveRDS(ll_removed, file = "data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
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
  # remain_df <- pcod[remain,]
  # lost_df <- pcod[lost,]
  # plot(pcod$X, pcod$Y, col = "black")
  # points(lost_df$X, lost_df$Y, col = "red")
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

# library(ggplot2)
# aleut <- readRDS(f[grep("aleut", f)])$survey_sets
# aleut <- filter(aleut, survey_abbrev == "HBLL OUT N", year == 2015)
# aleut <- assign_restricted_tows_hbll(aleut)
# ggplot(aleut, aes(longitude, latitude, size = density_ppkm2, colour = restricted)) + geom_point()

if (Sys.info()[['user']] == "seananderson") {
  f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
    full.names = TRUE)
  hbll_data <- furrr::future_map_dfr(seq_along(f), function(i) {
    d <- readRDS(f[i])$survey_sets
    filter(d, survey_abbrev %in% c("HBLL OUT N")) %>%
      select(year, survey_abbrev, species_science_name, species_common_name,
        density_ppkm2, latitude, longitude, grouping_code, area_km2, depth_m)
  })
  saveRDS(hbll_data, file = "data-raw/hbll-survey-data.rds")
} else {
  hbll_data <- readRDS("data-raw/hbll-survey-data.rds")
}

dat <- assign_restricted_tows_hbll(hbll_data)

frac_pos_df <- dat %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(frac_pos = sum(density_ppkm2 > 0) / length(density_ppkm2), .groups = "drop")

frac_pos_df %>% ungroup() %>% arrange(-frac_pos) %>% top_n(30) %>% as.data.frame()

tokeep <- frac_pos_df %>% filter(frac_pos > 0.05)
notkeep <- frac_pos_df %>% filter(frac_pos <= 0.05)

nrow(frac_pos_df)
nrow(tokeep)

dat_to_fit <- left_join(tokeep, dat, by = c("species_common_name", "survey_abbrev"))

saveRDS(dat_to_fit, file = "data-generated/dat_to_fit_hbll.rds")

plan(sequential)

dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
hbll_grid <- gfplot::hbll_n_grid$grid
utm_zone9 <- 3156
coords <- hbll_grid %>%
  sf::st_as_sf(crs = 4326, coords = c("X", "Y")) %>%
  sf::st_transform(utm_zone9)
coords_restr <- shrink_a_survey(coords, ll_removed, plot = FALSE)
coords <- coords %>%
  sf::st_coordinates() %>%
  as.data.frame()
coords$X <- coords$X / 1000
coords$Y <- coords$Y / 1000
coords$restricted <- coords_restr$restricted
grid <- coords %>%
  expand_prediction_grid(years = sort(unique(dat_to_fit$year))) %>%
  as_tibble() %>%
  arrange(year, X, Y)

saveRDS(grid, "data-generated/hbll-n-grid-w-restr.rds")

dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
syn_grid <- gfplot::synoptic_grid
syn_grid$X <- syn_grid$X * 1000
syn_grid$Y <- syn_grid$Y * 1000
utm_zone9 <- 3156
coords <- syn_grid %>%
  sf::st_as_sf(crs = utm_zone9, coords = c("X", "Y"))
coords_restr <- shrink_a_survey(coords, trawl_removed, plot = FALSE)
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
