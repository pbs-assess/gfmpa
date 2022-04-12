# calculate design-based biomass estimate from output of get_survey_sets()
calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>%
    group_by(year, survey_id, area_km2, grouping_code) %>%
    summarise(density = mean(density_kgpm2 * 1e6), .groups = "drop_last") %>%
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
    cv = sd(b$t) / mean(b$t),
    biomass = calc_bio(x)
  )
}

boot_wrapper <- function(dat, reps) {
  out <- dat %>%
    split(dat$year) %>%
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

boot_wrapper_parallel <- function(dat, reps) {
  out <- dat %>%
    split(dat$year) %>%
    furrr::future_map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

expand_prediction_grid <- function(grid, years) {
  nd <- do.call(
    "rbind",
    replicate(length(years), grid, simplify = FALSE)
  )
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
    remain_df <- .d[remain, ]
    lost_df <- .d[lost, ]
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
                          MPA_trend,
                          # formula = response ~ 0 + as.factor(year) + offset,
                          family = sdmTMB::tweedie(link = "log"),
                          return_model = FALSE,
                          ...) {
  if ("hook_count" %in% names(surv_dat)) {
    survey_type <- "HBLL"
  } else {
    survey_type <- "SYN"
  }
  mesh <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  if (MPA_trend) {
    surv_dat <- surv_dat %>% mutate(
      MPA = as.integer(restricted),
      year_seq = year - min(surv_dat$year)
    )
    formula <- response ~ 1 + year_seq + MPA + year_seq:MPA
  } else {
    formula <- response ~ 0 + as.factor(year)
  }

  if (survey_type == "HBLL") {
    surv_dat$offset <- log(surv_dat$hook_count)
    # surv_dat$offset <- log(surv_dat$adjusted_hooks)
    surv_dat$response <- surv_dat$catch_count
  } else if (survey_type == "SYN") {
    surv_dat$offset <- log(surv_dat$tow_length_m * surv_dat$doorspread_m)
    surv_dat$response <- surv_dat$catch_weight
  } else {
    stop("Survey type not found", call. = FALSE)
  }

  m <- try({
    sdmTMB(
      formula = formula,
      data = surv_dat,
      family = family,
      time = "year",
      # silent = F,
      offset = surv_dat$offset,
      mesh = mesh,
      control = sdmTMBcontrol(newton_loops = 1L),
      ...
    )
  })

  if (class(m)[[1]] == "try-error") {
    return(NULL)
  }
  if (max(m$gradients) > 0.01) {
    m <- try({
      run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 1L)
    })
  }
  if (return_model) {
    return(m)
  }
  if (class(m)[[1]] == "try-error" || max(m$gradients) > 0.01) {
    return(NULL)
  }
  m
}

fit_geo_model <- function(surv_dat, pred_grid,
                          MPA_trend = FALSE,
                          mpa_dat_removed = FALSE,
                          survey = c("HBLL", "SYN"),
                          family = c(sdmTMB::tweedie(), sdmTMB::delta_gamma(), sdmTMB::nbinom2()),
                          return_model = FALSE, ...) {
  survey <- match.arg(survey)

  pred_grid <- dplyr::filter(pred_grid, survey_abbrev == surv_dat$survey_abbrev[1])
  pred_grid <- dplyr::filter(pred_grid, year %in% unique(surv_dat$year))

  cutoff <- if (survey == "HBLL") 10 else 10

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )

  cat("Fitting", survey, unique(surv_dat$species_common_name), "\n")

  .sp <- gsub(" ", "-", unique(surv_dat$species_science_name)[1])
  .sp <- gsub("/", "-", .sp)

  dir.create("data-generated/model-cache", showWarnings = FALSE)
  dir.create("data-generated/index-cache", showWarnings = FALSE)
  .file_model <- paste0(
    "data-generated/model-cache/model-",
    gsub(" ", "-", unique(surv_dat$survey_abbrev)), "-", .sp, "-",
    if (mpa_dat_removed) "-mpa-dat-removed-",
    paste(family$family, collapse = "-"), ".rds"
  )
  .file_ind <- gsub("model", "index", .file_model)

  if (!file.exists(.file_model)) {
    fit <- do_sdmTMB_fit(
      surv_dat,
      MPA_trend = MPA_trend,
      cutoff = cutoff,
      family = family,
      pred_grid = pred_grid,
      return_model = return_model,
      ...
    )
    saveRDS(fit, file = .file_model)
  } else {
    fit <- readRDS(.file_model)
    fit$tmb_obj$retape()
  }

  if (is.null(fit)) {
    return(null_df)
    message("NULL model; discarding.")
  }

  .diag <- sdmTMB:::get_convergence_diagnostics(fit$sd_report)
  if (max(.diag$final_grads) > 0.01 || isTRUE(.diag$bad_eig) || isFALSE(.diag$pdHess)) {
    message("Didn't converge; discarding.")
    return(null_df)
  }

  if (!file.exists(.file_ind)) {
    pred <- try({
      predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
    })

    if (return_model) {
      return(pred$data)
    }
    if (class(pred)[[1]] == "try-error") {
      return(NULL)
      message("Error on prediction; discarding.")
    }

    ind <- get_index(pred, area = 4, bias_correct = TRUE) # 2 x 2 km
    ind$region <- "all"
    # ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
    #   geom_line() + geom_ribbon(alpha = 0.2)

    if (length(unique(pred_grid$restricted)) > 1) {
      mpa_only <- pred_grid[pred_grid$restricted, , drop = FALSE]
      pred_mpa_only <- try({
        predict(fit, newdata = mpa_only, return_tmb_object = TRUE)
      })
      if (class(pred)[[1]] != "try-error") {
        ind2 <- get_index(pred_mpa_only, area = 4, bias_correct = TRUE)
        ind2$region <- "mpa"
        ind <- dplyr::bind_rows(ind, ind2)
      }
    }

    ind$pdHess <- .diag$pdHess
    ind$species_common_name <- surv_dat$species_common_name[1]
    ind$survey_abbrev <- surv_dat$survey_abbrev[1]
    saveRDS(ind, file = .file_ind)
  } else {
    ind <- readRDS(.file_ind)
  }
  ind
}


sim_mpa_surv <- function(surv_dat, grid,
                         survey = "HBLL",
                         # not sure if delta-gamma possible to sim effectively?
                         # so for now hard coded to tweedie
                         # family = sdmTMB::tweedie(link = "log"),
                         mpa_increase_per_year = log(1.05), ...) {
  SEED <- 2938
  set.seed(SEED)

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )

  surv_dat <- surv_dat %>% mutate(
    depth_mean = mean(surv_dat$depth),
    depth_sd = sd(surv_dat$depth),
    depth_scaled = (depth - depth_mean) / depth_sd,
    depth_scaled2 = depth_scaled^2
  )
  # if (survey == "HBLL") surv_dat$density <- surv_dat$density_ppkm2
  # if (survey == "SYN") surv_dat$density <- surv_dat$density_kgpm2 * 1000
  # surv_dat$response <- surv_dat$density
  # cutoff <- if (survey == "HBLL") 10 else 15
  # m <- do_sdmTMB_fit(surv_dat, cutoff,
  #                    MPA_trend = FALSE,
  #                    family = sdmTMB::tweedie(link = "log"),
  #                    return_model = TRUE, ...)

  m0 <- fit_geo_model(surv_dat,
    pred_grid = grid, survey = survey,
    family = "tweedie",
    MPA_trend = FALSE,
    return_model = TRUE, do_fit = FALSE
  )

  m <- try({
    sdmTMB::sdmTMB(
      # response ~ 1 + s(year, k = 5), #didn't fit as well for dogfish HBLL
      response ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
      data = m0$data,
      mesh = m0$spde,
      time = "year",
      family = sdmTMB::tweedie(link = "log"),
      # silent = FALSE,
      return_model = TRUE, ...
    )
  })
  # m

  if (class(m)[[1]] == "try-error") {
    return(null_df)
  }
  if (max(m$gradients) > 0.01) {
    m <- try({
      run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 0L)
    })
  }
  if (class(m)[[1]] == "try-error" || max(m$gradients) > 0.01) {
    return(null_df)
  }

  b1 <- tidy(m)
  b2 <- tidy(m, "ran_pars")
  b <- bind_rows(b1, b2)

  dat <- m$data
  dat$year_covariate <- dat$year - min(dat$year)
  pars <- m$tmb_obj$env$last.par.best
  omega_s <- pars[names(pars) == "omega_s"]
  epsilon_st <- matrix(pars[names(pars) == "epsilon_st"], ncol = length(unique(dat$year)))

  s <- sdmTMB::sdmTMB_simulate(
    formula = ~ 1 + depth_scaled + depth_scaled2 + restricted * year_covariate,
    data = dat,
    mesh = m$spde,
    family = sdmTMB::tweedie(link = "log"),
    time = "year",
    # rho = 2 * plogis(m$model$par[['ar1_phi']]) - 1, # TODO!?
    rho = NULL,
    sigma_O = b$estimate[b$term == "sigma_O"],
    sigma_E = b$estimate[b$term == "sigma_E"],
    phi = b$estimate[b$term == "phi"],
    tweedie_p = b$estimate[b$term == "tweedie_p"],
    range = b$estimate[b$term == "range"],
    fixed_re = list(omega_s = omega_s, epsilon_st = NULL, zeta_s = NULL),
    # (Intercept), depth_scaled, depth_scaled2, restrictedTRUE, year_covariate, restrictedTRUE:year_covariate
    B = c(
      mean(b[grep("year", b$term), "estimate"]),
      b[b$term == "depth_scaled", "estimate"],
      b[b$term == "depth_scaled2", "estimate"],
      0, 0, mpa_increase_per_year
    ),
    seed = SEED
  )
  s <- rename(s, restricted = restrictedTRUE)

  m2 <- try({
    sdmTMB::sdmTMB(
      observed ~ 1 + # depth_scaled + depth_scaled2 +
        restricted * year_covariate,
      data = s,
      mesh = m$spde,
      time = "year",
      family = sdmTMB::tweedie(link = "log"),
      # silent = FALSE,
      # ...
      spatial = "on",
      spatiotemporal = "IID"
      # spatial = "off",
      # spatiotemporal = "AR1"
    )
  })

  if (class(m2)[[1]] == "try-error") {
    return(null_df)
  }
  if (max(m2$gradients) > 0.01) {
    m2 <- try({
      run_extra_optimization(m2, newton_loops = 1L, nlminb_loops = 0L)
    })
  }
  if (class(m2)[[1]] == "try-error" || max(m2$gradients) > 0.05) {
    return(null_df)
  }

  pred_grid <- grid %>% mutate(
    depth_scaled = (depth - surv_dat$depth_mean[1]) / surv_dat$depth_sd[1],
    depth_scaled2 = depth_scaled^2
  )

  pred_grid <- filter(pred_grid, year %in% unique(dat$year))
  pred_grid$year_covariate <- pred_grid$year - min(pred_grid$year)

  p <- try({
    predict(m2, newdata = pred_grid, sims = 300L)
  })
  if (class(p)[[1]] == "try-error") {
    return(null_df)
  }

  p_mpa <- p[pred_grid$restricted, ]
  attr(p_mpa, "time") <- "year"

  p_out <- p[!pred_grid$restricted, ]
  attr(p_out, "time") <- "year"

  ind_mpa <- get_index_sims(p_mpa)
  ind_out <- get_index_sims(p_out)

  b3 <- tidy(m2, conf.int = T)

  # # could add 4. naive GLM BACI on raw observations?
  # m3 <- mgcv::gam(observed ~ restricted * year_covariate, data = m2$data, family = tw())
  # # above ignores spatial correlation and temporal autocorrelation and sampling location
  # summary(m3)

  eta <- s %>%
    group_by(restricted, year) %>%
    suppressMessages(summarise(mean_eta = mean(eta))) %>%
    ungroup()

  ind_mpa <- ind_mpa %>%
    suppressMessages(left_join(., filter(eta, restricted == 1))) %>%
    mutate(
      type = "Restricted",
      true_slope = mpa_increase_per_year,
      m_slope = b3[b3$term == "restricted:year_covariate", 2],
      m_slope_se = b3[b3$term == "restricted:year_covariate", 3],
      m_slope_lwr = b3[b3$term == "restricted:year_covariate", 4],
      m_slope_upr = b3[b3$term == "restricted:year_covariate", 5]
    )

  ind_out <- ind_out %>%
    suppressMessages(left_join(., filter(eta, restricted == 0))) %>%
    mutate(
      type = "Not restricted",
      true_slope = 0,
      m_slope = b3[b3$term == "year_covariate", 2],
      m_slope_se = b3[b3$term == "year_covariate", 3],
      m_slope_lwr = b3[b3$term == "year_covariate", 4],
      m_slope_upr = b3[b3$term == "year_covariate", 5]
    )
  ind <- bind_rows(ind_mpa, ind_out)

  # ind$true_slope <- mpa_increase_per_year
  # ind$m_slope <- b3[b3$term=="restricted:year_covariate",2]
  # ind$m_slope_se <- b3[b3$term=="restricted:year_covariate",3]
  ind$species_common_name <- surv_dat$species_common_name[1]
  ind$survey_abbrev <- surv_dat$survey_abbrev[1]
  ind
  # # not sure what to save?
  # list(
  #   mpa_increase_per_year = mpa_increase_per_year,
  #   m = m2, # model includes simulated data
  #   index = ind)
}
