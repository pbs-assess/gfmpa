library(dplyr)
library(ggplot2)
library(sf)
library(future)
plan(multisession)
options(future.rng.onMisuse = "ignore")

dir.create("data-generated", showWarnings = FALSE)

x <- sf::read_sf("~/Downloads/MPAnetwork_Working_Draft_Scen2_VO1_IOAC_20210610.gdb/")
# plot(x["hu_co_demersalfishing_bottomtrawling_d"])
trawl_removed <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("X"))
saveRDS(trawl_removed, file = "data-generated/demersalfishing_bottomtrawling_X.rds")

# plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

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

#calculate design-based biomass estimate from output of get_survey_sets()
calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>% group_by(year, survey_id, area_km2, grouping_code) %>%
    summarise(density = mean(density_kgpm2*1e6), .groups = "drop_last") %>%
    group_by(year) %>%
    summarise(biomass = sum(density * area_km2), .groups = "drop_last") %>%
    pull(biomass)
}

boot_one_year <- function(x, reps) {
  b <- boot::boot(x, statistic = calc_bio, strata = x$grouping_code, R = reps)
  suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
  tibble::tibble(
    index = mean(b$t),
    median_boot = median(b$t),
    lwr = bci$percent[[4]],
    upr = bci$percent[[5]],
    cv = sd(b$t)/mean(b$t),
    biomass = calc_bio(x))
}

boot_wrapper <- function(dat, reps) {
  out <- dat %>% split(dat$year) %>%
    furrr::future_map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

f <- list.files("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/",
  full.names = TRUE)

out2 <- purrr::map_dfr(seq_along(f), function(i) {
  d <- readRDS(f[i])$survey_sets
  d <- filter(d, survey_abbrev %in% c("SYN QCS", "SYN HS", "SYN WCHG"))
  d2 <- assign_restricted_tows(d)

  cat(d$species_common_name[1], "\n")
  set.seed(1)
  out <- group_by(d2, survey_abbrev) %>%
    split(d2$survey_abbrev) %>%
    purrr::map_dfr(function(.x) boot_wrapper(.x, reps = 100L),
      .id = "survey_abbrev") %>%
    mutate(type = "Status quo")

  set.seed(1)
  out_restr <- d2 %>%
    filter(!restricted) %>%
    group_by(survey_abbrev) %>%
    split(.$survey_abbrev) %>%
    purrr::map_dfr(function(.x) boot_wrapper(.x, reps = 100L),
      .id = "survey_abbrev") %>%
    mutate(type = "Restricted")

  combined <- bind_rows(out, out_restr)
  combined$species_common_name <- d$species_common_name[1]
  combined
})

ggplot(out2, aes(year, biomass, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_grid(species_common_name~survey_abbrev) +
  geom_ribbon(alpha = 0.2, colour = NA)

#
# ggplot(out, aes(year, biomass, ymin = lwr, ymax = upr)) +
#   geom_line() +
#   geom_ribbon(alpha = 0.2)
#
# ii <- readRDS("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/pacific-cod.rds")$survey_index
#
# ii %>% filter(survey_abbrev %in% c("SYN QCS")) %>%
#   ggplot(aes(year, biomass, ymin = lowerci, ymax = upperci)) +
#   geom_line() +
#   geom_ribbon(alpha = 0.2)
