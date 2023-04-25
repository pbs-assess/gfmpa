library(tidyverse)
library(sdmTMB)
theme_set(theme_light())

# spp <- "big skate"
# spp <- "pacific cod"
# survey <- "SYN WCHG"
SILENT <- TRUE

dir.create("figs/raw-data-maps", showWarnings = FALSE)
dir.create("figs/indexes", showWarnings = FALSE)
dir.create("figs/sim-data", showWarnings = FALSE)
dir.create("data-generated/sim-data", showWarnings = FALSE)
dir.create("data-generated/indexes", showWarnings = FALSE)

calc_indices <- function(spp, survey, force = FALSE) {
  survey <- unlist(strsplit(survey, "\\|"))
  cat(spp, "\n")
  cat(survey, "\n")
  spp_file <- gsub(" ", "-", spp)
  spp_file <- gsub("\\/", "-", spp_file)
  surv_file <- gsub(" ", "-", paste(survey, collapse = "-"))
  ggplot2::theme_set(ggplot2::theme_light())

  if (grepl("SYN", survey[[1]])) {
    surv_dat <- readRDS("data-generated/dat_to_fit.rds")
  } else {
    surv_dat <- readRDS("data-generated/dat_to_fit_hbll.rds")
  }

  grid <- readRDS("data-generated/grids-strata-restricted.rds")
  grid$area <- 4
  grid <- filter(grid, survey_abbrev %in% survey)

  surv_dat <- filter(surv_dat, survey_abbrev %in% survey, species_common_name == spp)
  if (nrow(surv_dat) == 0L) {
    return(NULL)
  }

  downdat <- readRDS("data-generated/downsampled-fitting-data.rds")
  surv_dat_down1 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 1)
  surv_dat_down2 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 2)
  surv_dat_down3 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 3)
  surv_dat_down4 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 4)
  surv_dat_down5 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 5)
  rm(downdat)

  prep_cols_syn <- function(d) {
    d$area_swept1 <- d$doorspread_m * d$tow_length_m
    d$area_swept2 <- d$doorspread_m * d$duration_min * d$speed_mpm
    d$area_swept <- ifelse(!is.na(d$tow_length_m),
      d$area_swept1, d$area_swept2
    )
    d <- dplyr::filter(d, !is.na(area_swept))
    d$offset <- log(d$area_swept * 0.00001)
    d$response <- d$catch_weight
    d
  }
  prep_cols_hbll <- function(d) {
    d$offset <- log(d$hook_count)
    d$response <- d$catch_count
    d
  }

  if (grepl("SYN", survey[[1]])) {
    surv_dat <- prep_cols_syn(surv_dat)
    # surv_dat_up <- prep_cols_syn(surv_dat_up)
    surv_dat_down1 <- prep_cols_syn(surv_dat_down1)
    surv_dat_down2 <- prep_cols_syn(surv_dat_down2)
    surv_dat_down3 <- prep_cols_syn(surv_dat_down3)
    surv_dat_down4 <- prep_cols_syn(surv_dat_down4)
    surv_dat_down5 <- prep_cols_syn(surv_dat_down5)
  } else {
    surv_dat <- prep_cols_hbll(surv_dat)
    # surv_dat_up <- prep_cols_hbll(surv_dat_up)
    surv_dat_down1 <- prep_cols_hbll(surv_dat_down1)
    surv_dat_down2 <- prep_cols_hbll(surv_dat_down2)
    surv_dat_down3 <- prep_cols_hbll(surv_dat_down3)
    surv_dat_down4 <- prep_cols_hbll(surv_dat_down4)
    surv_dat_down5 <- prep_cols_hbll(surv_dat_down5)
  }

  g <- ggplot(surv_dat, aes(X, Y, colour = restricted, size = response / exp(offset))) +
    geom_point(pch = 21) +
    facet_wrap(~year) +
    coord_fixed() +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "grey60")) +
    ggtitle(spp) +
    scale_size_area(max_size = 8)
  ggsave(paste0("figs/raw-data-maps/", spp_file, "-", surv_file, ".pdf"),
    width = 10, height = 10
  )

  if (!file.exists(paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds")) || force) {
    priors <- sdmTMBpriors(
      matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
      matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
    )
    # priors <- sdmTMBpriors()

    surv_dat_r <- filter(surv_dat, !restricted)

    cutoff <- if (survey == "SYN WCHG") 5 else 8

    mesh_all <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
    mesh_restr <- make_mesh(surv_dat_r, xy_cols = c("X", "Y"), cutoff = cutoff, mesh = mesh_all$mesh) # helps convergence not to use full mesh!?
    # mesh_up <- make_mesh(surv_dat_up, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    # mesh_down <- make_mesh(surv_dat_down, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    # mesh_down2 <- make_mesh(surv_dat_down2, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)

    mi <- list(
      spatiotemporal = list("iid", "iid"),
      family = sdmTMB::delta_gamma(),
      anisotropy = FALSE
    )

    # relax default gradient_thresh from 0.001 to 0.01:
    sanity <- function(...) sdmTMB::sanity(..., gradient_thresh = 0.01)

    # fit_restr2 <- try({
    #   sdmTMB(
    #     formula = response ~ 0 + as.factor(year),
    #     data = surv_dat_r,
    #     family = tweedie(),
    #     time = "year",
    #     spatiotemporal = "iid",
    #     offset = "offset",
    #     mesh = mesh_restr,
    #     anisotropy = FALSE,
    #     priors = priors,
    #     silent = SILENT,
    #     control = sdmTMBcontrol(newton_loops = 1L),
    #   )})

    cat("Fit restricted\n")

    #####################
#
#     mesh <- make_mesh(surv_dat_r, xy_cols = c("X", "Y"), cutoff = 5)
#     plot(mesh$mesh, asp = 1)
#
#     # sum(is.na(surv_dat_r$depth_m))
#     m_aniso <- sdmTMB(
#       formula = response ~ 0 + as.factor(year),
#       data = surv_dat_r,
#       family = mi$family,
#       time = "year",
#       spatiotemporal = mi$spatiotemporal,
#       offset = "offset",
#       mesh = mesh,
#       anisotropy = T,
#       # priors = priors,
#       silent = F,
#       control = sdmTMBcontrol(newton_loops = 1L),
#     )
    # m_iso <- update(m_aniso, anisotropy = FALSE)
    #
    # sanity(m_iso)
    # sanity(m_aniso)
    # AIC(m_iso, m_aniso)

    #####################

    # fit_restr_aniso <- try({
    #  sdmTMB(
    #     formula = response ~ 0 + as.factor(year),
    #     data = surv_dat_r,
    #     family = mi$family,
    #     time = "year",
    #     spatiotemporal = mi$spatiotemporal,
    #     offset = "offset",
    #     mesh = mesh_restr,
    #     anisotropy = TRUE,
    #     priors = priors,
    #     silent = SILENT
    #   )
    # })
    fit_restr_iso <- try({
     sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat_r,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_restr,
        anisotropy = FALSE,
        priors = priors,
        silent = SILENT
      )
    })

    # s_aniso <- all(unlist(sanity(fit_restr_aniso)))
    # s_iso <- all(unlist(sanity(fit_restr_iso)))
    # if (s_aniso && s_iso) {
    #   delta_AIC <- AIC(fit_restr_aniso) - AIC(fit_restr_iso)
    #   if (delta_AIC < -2) {
    #     mi$anisotropy <- TRUE
    #   } else {
    #     mi$anisotropy <- FALSE
    #   }
    # }
    # if (s_aniso && !s_iso) {
    #   mi$anisotropy <- TRUE
    # }
    # if (!s_aniso && s_iso) {
    #   mi$anisotropy <- FALSE
    # }
    # if (!s_aniso && !s_iso) {
    #   mi$anisotropy <- FALSE
    # }
    # if (mi$anisotropy) {
      # fit_restr <- fit_restr_aniso
    # } else {
      fit_restr <- fit_restr_iso
    # }
    ok <- all(unlist(sanity(fit_restr)))

    if (!ok) {
      mi$spatiotemporal <- list("off", "iid")
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      ok <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$family <- sdmTMB::tweedie()
      mi$spatiotemporal <- "iid"
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      ok <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$family <- sdmTMB::delta_gamma()
      mi$spatiotemporal <- list("off", "off")
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      ok <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$spatiotemporal <- "off"
      mi$family <- sdmTMB::tweedie()
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
    }
    sanity_restr <- all(unlist(sanity(fit_restr)))

    cat("Fit all\n")
    fit_all <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_all,
        anisotropy = mi$anisotropy,
        priors = priors,
        silent = SILENT
      )
    })
    sanity_all <- all(unlist(sanity(fit_all)))

    if (!sanity_all || !sanity_restr) {
      saveRDS(NULL, paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
      return(NULL)
    }

    cat("Fit downsamples\n")

    fit_down_model <- function(.dat) {
      mesh_down <- make_mesh(.dat, xy_cols = c("X", "Y"), cutoff = cutoff,  mesh = mesh_all$mesh)
      fit_down <- try({
        sdmTMB(
          formula = response ~ 0 + as.factor(year),
          data = .dat,
          family = mi$family,
          time = "year",
          spatiotemporal = mi$spatiotemporal,
          offset = "offset",
          mesh = mesh_down,
          anisotropy = mi$anisotropy,
          priors = priors,
          silent = SILENT
        )
      })
      sanity_down <- all(unlist(sanity(fit_down)))
      list(model = fit_down, sanity = sanity_down)
    }
    fit_down1 <- fit_down_model(surv_dat_down1)
    fit_down2 <- fit_down_model(surv_dat_down2)
    fit_down3 <- fit_down_model(surv_dat_down3)
    # fit_down4 <- fit_down_model(surv_dat_down4)
    # fit_down5 <- fit_down_model(surv_dat_down5)

    # up-sample -------------

    # Steps:
    # - figure out how many points have been restricted each year/stratum
    # - sample an equivalent number without replacement from the remaining blocks each year
    # - this is a quick hack to avoid and assigning survey blocks to all existing
    #   survey at locations, which is tricky because sometimes they fall very
    #   slightly outside of the block

    # select random locations:
    generate_up_sampled_dataset <- function(.seed = 1) {
      set.seed(.seed)
      up_sampled_dat <- surv_dat |>
        group_by(survey_abbrev, year, grouping_code) |>
        group_split() |>
        purrr::map_dfr(function(xx) {
          number_restricted <- sum(xx$restricted)
          available_grid <- grid |> filter(grouping_code %in% xx$grouping_code, !restricted)
          replace <- if (number_restricted > nrow(available_grid)) TRUE else FALSE
          out <- available_grid |> sample_n(size = number_restricted, replace = replace)
          out <- mutate(out, year = unique(xx$year))
          out <- mutate(out, up_sample = TRUE)
          out <- mutate(out, offset = 0)
          bind_rows(out, filter(xx, !restricted) |> mutate(up_sample = FALSE))
        })

      # simulate observations at those locations:
      if (isTRUE(mi$family$delta)) {
        pred1 <- predict(fit_all, newdata = up_sampled_dat, model = 1L, type = "response", nsim = 1L, offset = up_sampled_dat$offset)
        pred2 <- predict(fit_all, newdata = up_sampled_dat, model = 2L, type = "response", nsim = 1L, offset = up_sampled_dat$offset)
        s1 <- rbinom(nrow(pred1), size = 1L, prob = pred1)
        shape <- exp(fit_all$model$par[["ln_phi"]])
        scale <- pred2 / shape
        s2 <- rgamma(nrow(pred2), shape = shape, scale = scale)
        sim <- s1 * s2
      } else {
        pred <- predict(fit_all, newdata = up_sampled_dat, type = "response", nsim = 1L)
        p <- stats::plogis(fit_all$model$par[["thetaf"]]) + 1
        dispersion <- exp(fit_all$model$par[["ln_phi"]])
        sim <- fishMod::rTweedie(n = nrow(pred), p = p, mu = pred, phi = dispersion)
      }
      up_sampled_dat$simulated_response <- sim
      up_sampled_dat$upsample_seed <- .seed
      up_sampled_dat$response <- ifelse(
        up_sampled_dat$up_sample,
        up_sampled_dat$simulated_response, up_sampled_dat$response
      )
      saveRDS(up_sampled_dat, paste0("data-generated/sim-data/", spp_file, "-", surv_file, "-seed", .seed, ".rds"))

      if (.seed == 1) { # just plot first one
        g <- filter(up_sampled_dat, !up_sample) |>
          ggplot(aes(response + 1, simulated_response + 1)) +
          geom_point() +
          scale_x_log10() +
          scale_y_log10() +
          coord_fixed()
        ggsave(paste0("figs/sim-data/", spp_file, "-", surv_file, ".pdf"),
          width = 7, height = 7
        )
      }
      up_sampled_dat
    }
    up_sampled_dat1 <- generate_up_sampled_dataset(1)
    up_sampled_dat2 <- generate_up_sampled_dataset(2)
    up_sampled_dat3 <- generate_up_sampled_dataset(3)
    up_sampled_dat4 <- generate_up_sampled_dataset(4)
    up_sampled_dat5 <- generate_up_sampled_dataset(5)

    cat("Fit upsamples\n")
    fit_upsample <- function(dat) {
      mesh_up <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff, mesh = mesh_all$mesh)
      fit_up <- try({
        sdmTMB(
          formula = response ~ 0 + as.factor(year),
          data = dat,
          family = mi$family,
          time = "year",
          spatiotemporal = mi$spatiotemporal,
          offset = "offset",
          mesh = mesh_up,
          anisotropy = mi$anisotropy,
          priors = priors,
          silent = SILENT,
          control = sdmTMBcontrol(newton_loops = 1L),
        )
      })
      sanity_up <- all(unlist(sanity(fit_up)))
      list(model = fit_up, sanity = sanity_up)
    }
    fit_up1 <- fit_upsample(up_sampled_dat1)
    fit_up2 <- fit_upsample(up_sampled_dat2)
    fit_up3 <- fit_upsample(up_sampled_dat3)
    # fit_up4 <- fit_upsample(up_sampled_dat4)
    # fit_up5 <- fit_upsample(up_sampled_dat5)

    # ----------------------

    grid <- sdmTMB::replicate_df(grid, "year", time_values = unique(surv_dat$year))

    gr_full <- dplyr::filter(grid, survey_abbrev %in% survey)
    gr_mpa <- dplyr::filter(grid, survey_abbrev %in% survey, restricted == TRUE)
    gr_remaining <- dplyr::filter(grid, !restricted, survey_abbrev %in% survey)

    p <- predict(fit_all, newdata = gr_full, return_tmb_object = TRUE)
    ind <- get_index(p, bias_correct = TRUE)

    p <- predict(fit_all, newdata = gr_remaining, return_tmb_object = TRUE)
    ind_remain <- get_index(p, bias_correct = TRUE)

    pmpa <- predict(fit_all, newdata = gr_mpa, return_tmb_object = TRUE)
    indmpa <- get_index(pmpa, bias_correct = TRUE)

    pmpa_restr <- predict(fit_restr, newdata = gr_mpa, return_tmb_object = TRUE)
    indmpa_restr <- get_index(pmpa_restr, bias_correct = TRUE)

    pr <- predict(fit_restr, newdata = gr_full, return_tmb_object = TRUE)
    indr <- get_index(pr, bias_correct = TRUE)

    prs <- predict(fit_restr, newdata = gr_remaining, return_tmb_object = TRUE)
    indrs <- get_index(prs, bias_correct = TRUE)

    get_ind_down <- function(obj) {
      if (obj$sanity) {
        p_down <- predict(obj$model, newdata = gr_full, return_tmb_object = TRUE) #<
        return(get_index(p_down, bias_correct = TRUE))
      }
    }

    ind_down1 <- get_ind_down(fit_down1)
    ind_down2 <- get_ind_down(fit_down2)
    ind_down3 <- get_ind_down(fit_down3)
    # ind_down4 <- get_ind_down(fit_down4)
    # ind_down5 <- get_ind_down(fit_down5)

    get_ind_up <- function(obj) {
      if (obj$sanity) {
        p_up <- predict(obj$model, newdata = gr_remaining, return_tmb_object = TRUE) #<
        return(get_index(p_up, bias_correct = TRUE))
      }
    }
    ind_up1 <- get_ind_up(fit_up1)
    ind_up2 <- get_ind_up(fit_up2)
    ind_up3 <- get_ind_up(fit_up3)
    # ind_up4 <- get_ind_up(fit_up4)
    # ind_up5 <- get_ind_up(fit_up5)

    i <- bind_rows(
      mutate(ind, type = "Status quo"),
      mutate(indr, type = "Restricted"),
      mutate(indrs, type = "Restricted and shrunk"),
      mutate(ind_remain, type = "Status quo and shrunk"),
      mutate(indmpa, type = "MPA only"),
      mutate(indmpa_restr, type = "MPA only restricted")
    )
    if (fit_down1$sanity) i <- bind_rows(i, mutate(ind_down1, type = "Random down-sampled 1"))
    if (fit_down2$sanity) i <- bind_rows(i, mutate(ind_down2, type = "Random down-sampled 2"))
    if (fit_down3$sanity) i <- bind_rows(i, mutate(ind_down3, type = "Random down-sampled 3"))
    # if (fit_down4$sanity) i <- bind_rows(i, mutate(ind_down4, type = "Random down-sampled 4"))
    # if (fit_down5$sanity) i <- bind_rows(i, mutate(ind_down5, type = "Random down-sampled 5"))
    if (fit_up1$sanity) i <- bind_rows(i, mutate(ind_up1, type = "Random up-sampled and shrunk 1"))
    if (fit_up2$sanity) i <- bind_rows(i, mutate(ind_up2, type = "Random up-sampled and shrunk 2"))
    if (fit_up3$sanity) i <- bind_rows(i, mutate(ind_up3, type = "Random up-sampled and shrunk 3"))
    # if (fit_up4$sanity) i <- bind_rows(i, mutate(ind_up4, type = "Random up-sampled and shrunk 4"))
    # if (fit_up5$sanity) i <- bind_rows(i, mutate(ind_up5, type = "Random up-sampled and shrunk 5"))

    i$species_common_name <- spp
    i$survey_abbrev <- paste(survey, collapse = ", ")
    i$family <- paste(as.character(mi$family$family), collapse = "-")
    i$anisotropy <- as.character(mi$anisotropy)
    i$spatiotemporal <- paste(as.character(mi$spatiotemporal), collapse = ", ")
    saveRDS(i, paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
  } else {
    i <- readRDS(paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
  }

  if (!is.null(i)) {
    g <- i |>
      filter(!grepl("MPA", type)) |>
      filter(!grepl("Random up-sampled and shrunk [2-9]+", type)) |> # only visualize seed 1
      filter(!grepl("Random down-sampled [2-9]+", type)) |> # only visualize seed 1
      filter(type %in% c("Restricted and shrunk", "Status quo")) |>
      group_by(type) |>

      # mutate(lwr = lwr / exp(mean(log(est), na.rm = TRUE))) |>
      # mutate(upr = upr / exp(mean(log(est), na.rm = TRUE))) |>
      # mutate(est = est / exp(mean(log(est), na.rm = TRUE))) |>

      ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type, group = type)) +
      geom_pointrange(position = position_dodge(width = 1), pch = 21) +
      scale_colour_brewer(palette = "Dark2") +
      ylab("Index") +
      xlab("Year") +
      labs(colour = "Type") +
      ggtitle(spp)

    suppressWarnings({
      g <- g +
        geom_smooth(method = glm, se = FALSE, method.args = list(family = Gamma(link = "log")), formula = y ~ x)

      g1 <- g + coord_cartesian(expand = FALSE, ylim = c(0, NA))
      g2 <- g + scale_y_log10() + ylab("Index (log distributed)")
      g0 <- cowplot::plot_grid(g1, g2, nrow = 2L)
      ggsave(paste0("figs/indexes/", spp_file, "-", surv_file, ".pdf"),
        width = 11, height = 7, plot = g0
      )
    })
  }

  # in case force = FALSE and we don't have it in memory:
  if (!exists("fit_all") && "family" %in% names(i)) {
    mi <- list()
    if (i$family[1] == "tweedie") mi$family <- tweedie()
    if (grepl("delta", i$family[1])) mi$family <- delta_gamma()
    st <- as.list(strsplit(i$spatiotemporal, "-")[[1]])
    if (length(st) == 1L) st <- st[[1]]
    mi$spatiotemporal <- st
    mesh_all <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = if (survey == "SYN WCHG") 5 else 8)
    # priors <- sdmTMBpriors(
    #   matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
    #   matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
    # )
    priors <- sdmTMBpriors()
    grid <- sdmTMB::replicate_df(grid, "year", time_values = unique(surv_dat$year))
    gr_full <- dplyr::filter(grid, survey_abbrev %in% survey)
    fit_all <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_all,
        anisotropy = as.logical(i$anisotropy[1]),
        priors = priors,
        silent = SILENT
      )
    })
    sanity_all <- all(unlist(sanity(fit_all)))
    if (sanity_all) {
      p <- predict(fit_all, newdata = gr_full, return_tmb_object = TRUE)
    }
  }

  if (exists("sanity_all")) {
    if (sanity_all) {
      dir.create("fits/fitted-survey-maps", showWarnings = FALSE)
      if (!"est" %in% names(p$data)) {
        p$data$est <- log(plogis(p$data$est1) * exp(p$data$est2))
      }
      g <- ggplot(p$data, aes(X, Y, fill = est)) +
        geom_tile(width = 2.04, height = 2.04) +
        facet_wrap(~year) +
        scale_fill_viridis_c(option = "D") +
        geom_tile(data = filter(p$data, restricted), fill = "#00000035", width = 2, height = 2) +
        coord_fixed() +

        geom_point(data = filter(surv_dat, response > 0),
          mapping = aes(X, Y, colour = restricted, size = response / exp(offset)),
          pch = 21, inherit.aes = FALSE, alpha = 0.7) +
        geom_point(data = filter(surv_dat, response == 0),
          mapping = aes(X, Y, colour = restricted), pch = 4, inherit.aes = FALSE,
          alpha = 0.7, size = 0.6) +
        scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "white")) +
        ggtitle(spp) +
        scale_size_area(max_size = 8)

      ggsave(paste0("figs/fitted-survey-maps/", spp_file, "-", surv_file, ".pdf"),
        width = 14, height = 10, plot = g
      )
    }
  }

  return(NULL)
}

source("analysis/spp.R")

# syn_survs <- c("SYN WCHG", "SYN QCS|SYN HS", "SYN QCS", "SYN HS")
syn_survs <- c("SYN WCHG", "SYN QCS", "SYN HS")

library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
cores <- round(parallel::detectCores() / 2)
(cores <- parallel::detectCores() - 2L)
if (!is_rstudio && is_unix) plan(multicore, workers = cores) else plan(multisession, workers = cores)

to_fit <- expand_grid(spp = syn_highlights, survey = syn_survs)

## test: --------------------------------------------------
calc_indices(spp = "redstripe rockfish", survey = "SYN WCHG", force = T)
SILENT <- T
calc_indices(spp = "pacific ocean perch", survey = "SYN WCHG", force = T)

calc_indices(spp = syn_highlights[17], survey = syn_survs[3], force = F)
x <- readRDS("data-generated/indexes/pacific-cod-SYN-QCS.rds")
x |>
  filter(type != "MPA only", type != "MPA only restricted", type != "Restricted") |>
  ggplot(aes(year, se, colour = type)) +
  geom_point() +
  geom_line()
x |>
  filter(type != "MPA only", type != "MPA only restricted", type != "Restricted") |>
  ggplot(aes(year, est, colour = type)) +
  geom_point() +
  geom_line()
x <- readRDS("data-generated/sim-data/pacific-cod-SYN-QCS-seed1.rds")
glimpse(x)
x <- readRDS("data-generated/sim-data/pacific-cod-SYN-QCS-seed5.rds")
glimpse(x)
grid <- readRDS("data-generated/grids-strata-restricted.rds")
grid <- filter(grid, survey_abbrev %in% "SYN QCS")
pal <- c(as.character(colorBlindness::availableColors())[-1], c("grey60"))
x |> filter(year %in% c(2003, 2009, 2017)) |>
  ggplot(aes(X, Y, shape = up_sample)) +
  geom_tile(data = grid, mapping = aes(X, Y, colour = as.factor(grouping_code),
    fill = as.factor(grouping_code)), inherit.aes = FALSE, width = 2.01, height = 2.01) +
  geom_tile(data = filter(grid, restricted), mapping = aes(X, Y),
    inherit.aes = FALSE, width = 2, height = 2, fill = "grey60") +
  geom_point(alpha = 0.8, mapping = aes(size = response)) +
  facet_wrap(~year) +
  scale_shape_manual(values = c(21, 19)) +
  scale_size_area(max_size = 12) +
  coord_equal() +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(colour = "Stratum", fill = "Stratum", shape = "Up sampled",
    size = "Observed or\nsimulated catch") +
  theme(legend.position = "bottom")
ggsave("figs/upsample-example.pdf", width = 9, height = 4.2)
## end test -----------------------------

set.seed(123)
to_fit <- to_fit[sample(nrow(to_fit)),] # randomize for parallel
# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices, force = TRUE)

to_fit <- expand_grid(spp = hbll_highlights, survey = "HBLL OUT N")
## test:
# calc_indices(spp = hbll_highlights[1], survey = "HBLL OUT N")
# x <- readRDS("data-generated/indexes/arrowtooth-flounder-HBLL-OUT-N.rds")
# x |> filter(type != "MPA only", type != "MPA only restricted", type != "Restricted") |>
#   ggplot(aes(year, se, colour = type)) + geom_point() + geom_line()

# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices, force = TRUE)

plan(sequential)

f <- list.files("data-generated/indexes/", pattern = ".rds", full.names = TRUE)

ind <- purrr::map_dfr(f, readRDS)
ind$cv <- NULL
ind$cv <- sqrt(exp(ind$se^2) - 1)

table(ind$survey_abbrev, ind$type)

ocv <- ind |>
  filter(type == "Status quo") |>
  mutate(orig_cv = cv) |>
  select(year, orig_cv, species_common_name, survey_abbrev) |>
  distinct()

ind <- left_join(
  ind,
  ocv,
  by = join_by(year, species_common_name, survey_abbrev)
)

hbll <- ind[grepl("HBLL", ind$survey_abbrev), ]
syn <- ind[grepl("SYN", ind$survey_abbrev), ]
saveRDS(hbll, "data-generated/index-hbll-geo-clean.rds")
saveRDS(syn, "data-generated/index-syn-geo-clean.rds")

# g <- ind |>
#   # filter(grepl("SYN", survey_abbrev)) |>
#   ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type)) +
#   geom_pointrange(position = position_dodge(width = 0.55), pch = 21) +
#   scale_colour_manual(values =
#       c("Restricted" = "red", "Status quo" = "grey60", "Restricted and shrunk" = "purple")) +
#   ylab("Index") + xlab("Year") +
#   labs(colour = "Type") +
#   facet_grid(species_common_name~survey_abbrev, scales = "free_y")
# g <- g + scale_y_log10() + ylab("Index (log distributed)")
# ggsave("figs/giant-index-explore.pdf", width = 20, height = 50, limitsize = FALSE)
