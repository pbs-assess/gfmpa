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
                          MPA_trend, spatiotemporal = list("off", "iid"),
                          share_range = list(TRUE, TRUE),
                          # formula = response ~ 0 + as.factor(year) + offset,
                          family = sdmTMB::tweedie(link = "log"),
                          return_model = FALSE,
                          ...) {
  if ("hook_count" %in% names(surv_dat)) {
    survey_type <- "HBLL"
  } else {
    survey_type <- "SYN"
  }

  if (MPA_trend) {
    surv_dat <- surv_dat %>% mutate(
      MPA = as.integer(restricted),
      year_seq = year - min(surv_dat$year)
    )
    formula <- response ~ 1 + year_seq + MPA + year_seq:MPA
  } else {
    formula <- response ~ 0 + as.factor(year)
    # formula <- response ~ 1
  }

  if (survey_type == "HBLL") {
    surv_dat$offset <- log(surv_dat$hook_count)
    # surv_dat$offset <- log(surv_dat$adjusted_hooks)
    surv_dat$response <- surv_dat$catch_count
  } else if (survey_type == "SYN") {
    surv_dat$area_swept1 <- surv_dat$doorspread_m * surv_dat$tow_length_m
    surv_dat$area_swept2 <- surv_dat$doorspread_m * surv_dat$duration_min *
      surv_dat$speed_mpm
    surv_dat$area_swept <- ifelse(!is.na(surv_dat$tow_length_m),
      surv_dat$area_swept1, surv_dat$area_swept2)
    surv_dat <- dplyr::filter(surv_dat, !is.na(area_swept))
    surv_dat$offset <- log(surv_dat$area_swept * 0.00001)
    surv_dat$response <- surv_dat$catch_weight
  } else {
    stop("Survey type not found", call. = FALSE)
  }

  mesh <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  m <- try({
    sdmTMB(
      formula = formula,
      data = surv_dat,
      family = family,
      time = "year",
      spatiotemporal = spatiotemporal,
      offset = surv_dat$offset,
      mesh = mesh,
      anisotropy = FALSE,
      priors = sdmTMBpriors(
        matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
        matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
      ),
      silent = TRUE,
      share_range = share_range,
      control = sdmTMBcontrol(newton_loops = 1L),
      ...
    )
  })
  if (inherits(m, "try-error")) return(m)
  print(m)
  m
}

fit_geo_model <- function(surv_dat, pred_grid,
                          MPA_trend = FALSE,
                          mpa_dat_removed = FALSE, shrunk = FALSE,
                          survey = c("HBLL", "SYN"),
                          family = sdmTMB::delta_gamma(),
                          return_model = FALSE, ...) {
  survey <- match.arg(survey)

  pred_grid <- dplyr::filter(pred_grid, survey_abbrev == surv_dat$survey_abbrev[1])
  pred_grid <- dplyr::filter(pred_grid, year %in% unique(surv_dat$year))

  cutoff <- if (survey == "HBLL") 8 else 8

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )

  cat(
    "Fitting", unique(surv_dat$survey_abbrev),
    unique(surv_dat$species_common_name), "\n"
  )

  .sp <- gsub(" ", "-", unique(surv_dat$species_science_name)[1])
  .sp <- gsub("/", "-", .sp)

  .surv_name <- gsub(" ", "-", unique(surv_dat$survey_abbrev))

  dir.create("data-generated/model-cache", showWarnings = FALSE)
  dir.create("data-generated/index-cache", showWarnings = FALSE)
  .file_model <- paste0(
    "data-generated/model-cache/model-", .surv_name, "-", .sp, "-",
    if (mpa_dat_removed) "mpa-dat-removed-",
    paste(family$family, collapse = "-"), ".rds"
  )
  .file_ind <- gsub("model", "index", .file_model)
  model_info_file <- paste0("data-generated/model-cache/", .sp, "-", .surv_name, ".rds")

  if (shrunk) {
    .file_ind <- gsub(
      "mpa-dat-removed", "mpa-dat-removed-shrunk", .file_ind
    )
  }

  if (!file.exists(.file_ind)) {
    if (file.exists(model_info_file)) {
      model_info <- readRDS(model_info_file)
      model_already_decided <- TRUE
    } else {
      model_info <- list(
        spatiotemporal = list("iid", "iid"),
        family = sdmTMB::delta_gamma(),
        share_range = list(TRUE, TRUE)
      )
      model_already_decided <- FALSE
    }

    # order of trials:
    # - delta-Gamma, spatiotemporal: iid, iid
    # - delta-Gamma, spatiotemporal: off, iid
    # - Tweedie, spatiotemporal: iid
    # - delta-Gamma, spatiotemporal: off, off
    # - Tweedie, spatiotemporal: off

    # - delta-Gamma, spatiotemporal: iid, iid
    # - or whatever was fit last time!
    fit <- do_sdmTMB_fit(
      surv_dat,
      MPA_trend = MPA_trend,
      cutoff = cutoff,
      family = model_info$family,
      share_range = model_info$share_range,
      spatiotemporal = model_info$spatiotemporal,
      pred_grid = pred_grid,
      return_model = return_model,
      ...
    )
    s <- sanity(fit)
    all_ok <- all(unlist(s))

    if (!all_ok && !model_already_decided) {
      # - delta-Gamma, spatiotemporal: off, iid
      model_info$spatiotemporal <- list("off", "iid")

      fit <- do_sdmTMB_fit(
        surv_dat,
        MPA_trend = MPA_trend,
        cutoff = cutoff,
        family = model_info$family,
        pred_grid = pred_grid,
        return_model = return_model,
        spatiotemporal = model_info$spatiotemporal,
        share_range = model_info$share_range,
        ...
      )
      s <- sanity(fit)
      all_ok <- all(unlist(s))
    }

    if (!all_ok && !model_already_decided) {
      # - Tweedie, spatiotemporal: iid
      model_info$family <- sdmTMB::tweedie()
      model_info$spatiotemporal <- "iid"
      model_info$share_range <- TRUE

      fit <- do_sdmTMB_fit(
        surv_dat,
        MPA_trend = MPA_trend,
        cutoff = cutoff,
        family = model_info$family,
        pred_grid = pred_grid,
        return_model = return_model,
        spatiotemporal = model_info$spatiotemporal,
        share_range = model_info$share_range,
        ...
      )
      s <- sanity(fit)
      all_ok <- all(unlist(s))
    }

    if (!all_ok && !model_already_decided) {
      # - delta-Gamma, spatiotemporal: off, off
      model_info$family <- sdmTMB::delta_gamma()
      model_info$spatiotemporal <- list("off", "off")
      model_info$share_range <- list(TRUE, TRUE)

      fit <- do_sdmTMB_fit(
        surv_dat,
        MPA_trend = MPA_trend,
        cutoff = cutoff,
        family = model_info$family,
        pred_grid = pred_grid,
        return_model = return_model,
        spatiotemporal = model_info$spatiotemporal,
        share_range = model_info$share_range,
        ...
      )
      s <- sanity(fit)
      all_ok <- all(unlist(s))
    }

    if (!all_ok && !model_already_decided) {
      # - Tweedie, spatiotemporal: off
      model_info$spatiotemporal <- "off"
      model_info$family <- sdmTMB::tweedie()
      model_info$share_range <- TRUE

      fit <- do_sdmTMB_fit(
        surv_dat,
        MPA_trend = MPA_trend,
        cutoff = cutoff,
        family = model_info$family,
        pred_grid = pred_grid,
        return_model = return_model,
        spatiotemporal = model_info$spatiotemporal,
        share_range = model_info$share_range,
        ...
      )
      s <- sanity(fit)
      all_ok <- all(unlist(s))
    }

    if (!model_already_decided) saveRDS(model_info, model_info_file)

    if (!all_ok) {
      return(null_df)
      message("Didn't converge; discarding.")
    }

    cat(
      "Indexing",
      unique(surv_dat$survey_abbrev),
      unique(surv_dat$species_common_name), "\n"
    )
    pred <- try({
      predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
    })

    if (return_model) {
      return(pred$data)
    }
    if (inherits(pred, "try-error")) {
      return(null_df)
      message("Error on prediction; discarding.")
    }

    ind <- get_index(pred, area = pred_grid$area, bias_correct = TRUE) # 2 x 2 km
    ind$region <- "all"
    # ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
    #   geom_line() + geom_ribbon(alpha = 0.2)

    # if (length(unique(pred_grid$restricted)) > 1) {
    #   mpa_only <- pred_grid[pred_grid$restricted, , drop = FALSE]
    #   pred_mpa_only <- try({
    #     predict(fit, newdata = mpa_only, return_tmb_object = TRUE)
    #   })
    #   if (!inherits(pred, "try-error")) {
    #     ind2 <- get_index(pred_mpa_only, area = mpa_only$area, bias_correct = TRUE)
    #     ind2$region <- "mpa"
    #     ind <- dplyr::bind_rows(ind, ind2)
    #   } else {
    #     return(null_df)
    #   }
    # }

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
