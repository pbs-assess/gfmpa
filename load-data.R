library(dplyr)
library(sf)
library(future)
plan(multisession)

dir.create("data-generated", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

source("functions.R")

x <- sf::read_sf("data-raw/Spatial_N1_May17_2021.gdb")
# plot(x["hu_co_demersalfishing_bottomtrawling_d"])
trawl_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("X"))
saveRDS(trawl_removed, file = "data-generated/demersalfishing_bottomtrawling_X.rds")
# plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
  full.names = TRUE)

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
  trawl_dat %>% as_tibble()
}

dat <- furrr::future_map_dfr(seq_along(f), function(i) {
  d <- readRDS(f[i])$survey_sets
  d <- filter(d, survey_abbrev %in% c("SYN QCS", "SYN HS", "SYN WCHG"))
  d2 <- assign_restricted_tows(d)
  d2
})

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

f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
  full.names = TRUE)

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
  dat %>% as_tibble()
}

library(ggplot2)
aleut <- readRDS(f[grep("aleut", f)])$survey_sets
aleut <- filter(aleut, survey_abbrev == "HBLL OUT N", year == 2015)
aleut <- assign_restricted_tows_hbll(aleut)
ggplot(aleut, aes(longitude, latitude, size = density_ppkm2, colour = restricted)) + geom_point()

dat <- furrr::future_map_dfr(seq_along(f), function(i) {
  d <- readRDS(f[i])$survey_sets
  d <- filter(d, survey_abbrev %in% c("HBLL OUT N"))
  d2 <- assign_restricted_tows_hbll(d)
  d2
})

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
