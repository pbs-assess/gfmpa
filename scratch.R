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

plot(trawl_removed["hu_co_demersalfishing_bottomtrawling_d"])

# trawl_empty <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c(""))
# plot(trawl_empty["hu_co_demersalfishing_bottomtrawling_d"])

trawl_empty_x <- dplyr::filter(x, hu_co_demersalfishing_bottomtrawling_d %in% c("", "X"))
plot(trawl_empty_x["hu_co_demersalfishing_bottomtrawling_d"])

assign_restricted_tows <- function(trawl_dat) {
  orig <- trawl_dat
  trawl_dat <- trawl_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(sf::st_crs(trawl_empty_x))
  intersected <- sf::st_intersects(trawl_dat, trawl_empty_x)
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
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

boot_wrapper_parallel <- function(dat, reps) {
  out <- dat %>% split(dat$year) %>%
    furrr::future_map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

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

index_syn <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    # cat(.x$species_common_name[1], "\n")
    set.seed(1)
    out <- split(.x, .x$survey_abbrev) %>%
      purrr::map_dfr(function(.y)
        boot_wrapper(.y, reps = 90L), .id = "survey_abbrev") %>%
      mutate(type = "Status quo")
    set.seed(1)
    out_restr <- .x %>%
      filter(!restricted) %>%
      split(.$survey_abbrev) %>%
      purrr::map_dfr(function(.y)
        boot_wrapper(.y, reps = 90L), .id = "survey_abbrev") %>%
      mutate(type = "Restricted")
    combined <- bind_rows(out, out_restr)
    combined$species_common_name <- .x$species_common_name[1]
    combined
  }, .progress = TRUE)

saveRDS(index_syn, file = "index-syn-x-empty.rds")

g <- ggplot(index_syn, aes(year, biomass, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
  geom_ribbon(alpha = 0.2, colour = NA)

ggsave("index-syn-restricted-x-empty.pdf", width = 8, height = 85, limitsize = FALSE)

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

no_na <- index_syn %>%
  group_by(survey_abbrev, species_common_name) %>%
  summarise(nas = sum(is.na(cv) > 0), zeros = sum(median_boot == 0)) %>%
  filter(nas == FALSE & zeros == 0)
nrow(no_na)

nrow(index_syn)
left_join(no_na, index_syn) %>% nrow()

x <- left_join(no_na, index_syn) %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  # summarise(cv = mean(cv, na.rm = TRUE), .groups = "drop_last") %>%
  summarise(cv_ratio = cv[type == "Restricted"] / cv[type == "Status quo"])

x %>%
  ggplot(aes(cv_ratio)) + facet_wrap(~survey_abbrev) + geom_histogram()

x %>% group_by(survey_abbrev) %>%
  summarise(mean_ratio = mean(cv_ratio))


# geostatistical??


expand_prediction_grid <- function(grid, years) {
  nd <- do.call("rbind",
    replicate(length(years), grid, simplify = FALSE))
  nd[["year"]] <- rep(years, each = nrow(grid))
  nd
}

library(sdmTMB)
geo_wrapper <- function(trawl_dat) {
  utm_zone9 <- 3156
  coords <- trawl_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  coords$X <- coords$X / 1000
  coords$Y <- coords$Y / 1000
  trawl_dat <- dplyr::bind_cols(trawl_dat, coords)
  trawl_dat$density <-   trawl_dat$density_kgpm2 * 1000
  mesh <- make_mesh(trawl_dat, c("X", "Y"), cutoff = 15)
  # mesh$mesh$n
  m <- sdmTMB(density ~ 0 + as.factor(year),
    data = trawl_dat,
    family = tweedie(),
    time = "year",
    spde = mesh
  )
  .grid <- gfplot::synoptic_grid %>%
    dplyr::filter(survey == unique(trawl_dat$survey_abbrev)) %>%
    expand_prediction_grid(years = unique(trawl_dat$year))
  set.seed(1)
  pred <- predict(m, newdata = .grid, xy_cols = c("X", "Y"), sims = 200L)
  ind <- get_index_sims(pred, area = rep(4, nrow(pred)))
  ind$species_common_name <- trawl_dat$species_common_name[1]
  ind$survey_abbrev <- trawl_dat$survey_abbrev[1]
  ind
}

index_syn_geo <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
  # purrr::map_dfr(function(.x) {
    cat(.x$species_common_name[1], .x$survey_abbrev[1], "\n")
    set.seed(1)
    out <- .x %>% geo_wrapper() %>%
      mutate(type = "Status quo")

    out_restr <- .x %>%
      geo_wrapper() %>%
      mutate(type = "Restricted")

    bind_rows(out, out_restr)
  }, .progress = TRUE)

saveRDS(index_syn, file = "index-syn-x-empty.rds")
