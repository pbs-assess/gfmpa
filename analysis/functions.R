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


shrink_a_survey <- function(grid_dat, restriction_dat, plot = FALSE) {
  orig <- grid_dat
  grid_dat <- grid_dat %>% sf::st_transform(sf::st_crs(restriction_dat))
  intersected <- sf::st_intersects(grid_dat, restriction_dat)
  remain <- which(lengths(intersected) == 0)
  lost <- which(lengths(intersected) > 0)
  if (plot) {
    .d <- as.data.frame(sf::st_coordinates(orig))
    remain_df <- .d[remain,]
    lost_df <- .d[lost,]
    plot(remain_df$X, remain_df$Y, col = "black")
    points(lost_df$X, lost_df$Y, col = "red")
  }
  grid_dat$restricted <- lengths(intersected) > 0
  grid_dat <- as.data.frame(grid_dat) %>% select(-geometry)
  grid_dat$longitude <- orig$longitude
  grid_dat$latitude <- orig$latitude
  grid_dat %>% as_tibble()
}

do_sdmTMB_fit <- function(surv_dat, cutoff, pred_grid,
  family = sdmTMB::tweedie(link = "log"), return_model = FALSE, ...) {

  mesh <- make_mesh(surv_dat, c("X", "Y"), cutoff = cutoff)
  # plot(mesh)
  # mesh$mesh$n
  m <- try({
    sdmTMB(response ~ 0 + as.factor(year),
      data = surv_dat,
      family = family,
      time = "year",
      spde = mesh,
      ...
    )
  })
  if (class(m)[[1]] == "try-error") return(NULL)
  if (max(m$gradients) > 0.01) {
    m <- try({
      run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 0L)
    })
  }
  if (return_model) return(m)
  if (class(m)[[1]] == "try-error" || max(m$gradients) > 0.01) return(NULL)
  set.seed(1)
  pred <- try({
    predict(m, newdata = pred_grid, xy_cols = c("X", "Y"), sims = 300L)
  })
  if (class(pred)[[1]] == "try-error") return(NULL)
  pred
}

fit_geo_model <- function(surv_dat, pred_grid, shrink_survey = FALSE,
  survey = c("HBLL", "SYN"),
  family = c("tweedie", "binomial-gamma"),
  return_model = FALSE, ...) {

  survey <- match.arg(survey)
  family <- match.arg(family)

  utm_zone9 <- 3156
  coords <- surv_dat %>%
    sf::st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    sf::st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  coords$X <- coords$X / 1000
  coords$Y <- coords$Y / 1000
  surv_dat <- dplyr::bind_cols(surv_dat, coords)

  if (survey == "HBLL") surv_dat$density <- surv_dat$density_ppkm2
  if (survey == "SYN") surv_dat$density <- surv_dat$density_kgpm2 * 1000 # for computational reasons

  pred_grid <- filter(pred_grid, survey_abbrev == surv_dat$survey_abbrev[1])
  pred_grid <- filter(pred_grid, year %in% unique(surv_dat$year))
  cutoff <- if (survey == "HBLL") 10 else 15

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )

  if (family == "tweedie") {
    surv_dat$response <- surv_dat$density
    pred <- do_sdmTMB_fit(surv_dat, cutoff = cutoff, family = tweedie(),
      pred_grid = pred_grid, return_model = return_model, ...)
    if (return_model) return(pred)
    if (is.null(pred)) return(null_df)
    ind <- get_index_sims(pred, area = rep(4, nrow(pred))) # 2 x 2 km
  } else {
    surv_dat$present <- ifelse(surv_dat$density > 0, 1, 0)
    surv_dat$response <- surv_dat$present
    surv_dat_pos <- dplyr::filter(surv_dat, density > 0)
    surv_dat_pos$response <- surv_dat_pos$density
    pred_grid <- filter(pred_grid, year %in% unique(surv_dat_pos$year)) # in case 'pos' is missing some
    pred_bin <- do_sdmTMB_fit(surv_dat, cutoff = cutoff, family = binomial(),
      pred_grid = pred_grid, return_model = return_model, ...)
    if (is.null(pred_bin)) return(null_df)
    pred_pos <- do_sdmTMB_fit(surv_dat_pos, cutoff = cutoff,
      family = Gamma(link = "log"), pred_grid = pred_grid, return_model = return_model, ...)
    if (is.null(pred_pos)) return(null_df)
    if (return_model) list(bin = pred_bin, pos = pred_pos)

    pred_combined <- log(plogis(pred_bin) * exp(pred_pos))
    ind <- get_index_sims(pred_combined, area = rep(4, nrow(pred_combined)))
  }
  ind$species_common_name <- surv_dat$species_common_name[1]
  ind$survey_abbrev <- surv_dat$survey_abbrev[1]
  ind
}
