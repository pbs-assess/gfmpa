library(dplyr)
library(sf)
library(future)
plan(multisession)

dir.create("data-generated", showWarnings = FALSE)

x <- sf::read_sf("~/Downloads/MPAnetwork_Working_Draft_Scen2_VO1_IOAC_20210610.gdb/")
# plot(x["hu_co_demersalfishing_bottomtrawling_d"])
trawl_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("X"))
saveRDS(trawl_removed, file = "data-generated/demersalfishing_bottomtrawling_X.rds")

# plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

# trawl_empty <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c(""))
# plot(trawl_empty["hu_co_demersalfishing_bottomtrawling_d"])

trawl_empty_x <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("", "X"))
# plot(trawl_empty_x["hu_co_demersalfishing_bottomtrawling_d"])

f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
  full.names = TRUE)

# aleut <- readRDS(f[grep("aleut", f)])$survey_sets
# aleut <- filter(aleut, survey_abbrev == "SYN HS", year == 2015)
# aleut <- assign_restricted_tows(aleut)
# ggplot(aleut, aes(longitude, latitude, size = density_kgpm2, colour = restricted)) + geom_point()

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

