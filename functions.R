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

expand_prediction_grid <- function(grid, years) {
  nd <- do.call("rbind",
    replicate(length(years), grid, simplify = FALSE))
  nd[["year"]] <- rep(years, each = nrow(grid))
  nd
}
