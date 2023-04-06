library(tidyverse)
library(sdmTMB)
theme_set(theme_light())

# spp <- "big skate"
# spp <- "pacific cod"
# survey <- "SYN WCHG"
SILENT <- TRUE

dir.create("figs/raw-data-maps", showWarnings = FALSE)
dir.create("figs/indexes", showWarnings = FALSE)
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

  # updat <- readRDS("data-generated/upsampled-fitting-data.rds")
  downdat <- readRDS("data-generated/downsampled-fitting-data.rds")

  surv_dat <- filter(surv_dat, survey_abbrev %in% survey, species_common_name == spp)
  # surv_dat_up <- filter(updat, survey_abbrev %in% survey, species_common_name == spp, upsample_seed == 1)
  surv_dat_down <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 1)
  surv_dat_down2 <- filter(downdat, survey_abbrev %in% survey, species_common_name == spp, downsample_seed == 2)

  if (nrow(surv_dat) == 0L) {
    return(NULL)
  }

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
    surv_dat_down <- prep_cols_syn(surv_dat_down)
    surv_dat_down2 <- prep_cols_syn(surv_dat_down2)
  } else {
    surv_dat <- prep_cols_hbll(surv_dat)
    # surv_dat_up <- prep_cols_hbll(surv_dat_up)
    surv_dat_down <- prep_cols_hbll(surv_dat_down)
    surv_dat_down2 <- prep_cols_hbll(surv_dat_down2)
  }

  g <- ggplot(surv_dat, aes(X, Y, colour = restricted, size = response / exp(offset))) +
    geom_point(pch = 21) +
    facet_wrap(~year) +
    coord_fixed() +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "grey60")) +
    ggtitle(spp)
  ggsave(paste0("figs/raw-data-maps/", spp_file, "-", surv_file, ".pdf"),
    width = 10, height = 10
  )

  if (!file.exists(paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds")) || force) {
    priors <- sdmTMBpriors(
      matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
      matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
    )

    surv_dat_r <- filter(surv_dat, !restricted)

    mesh_all <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = 10)
    mesh_restr <- make_mesh(surv_dat_r, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    # mesh_up <- make_mesh(surv_dat_up, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    mesh_down <- make_mesh(surv_dat_down, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    mesh_down2 <- make_mesh(surv_dat_down2, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)

    mi <- list(
      spatiotemporal = list("iid", "iid"),
      family = sdmTMB::delta_gamma()
    )

    fit_restr <- try({
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
        silent = SILENT,
        control = sdmTMBcontrol(newton_loops = 1L),
      )
    })
    s <- sanity(fit_restr)
    ok <- all(unlist(s))

    if (!ok) {
      mi$spatiotemporal <- list("off", "iid")
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      s <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$family <- sdmTMB::tweedie()
      mi$spatiotemporal <- "iid"
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      s <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$family <- sdmTMB::delta_gamma()
      mi$spatiotemporal <- list("off", "off")
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
      s <- all(unlist(sanity(fit_restr)))
    }
    if (!ok) {
      mi$spatiotemporal <- "off"
      mi$family <- sdmTMB::tweedie()
      fit_restr <- try({
        update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)
      })
    }
    sanity_restr <- all(unlist(sanity(fit_restr)))

    fit_all <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_all,
        anisotropy = FALSE,
        priors = priors,
        silent = SILENT,
        control = sdmTMBcontrol(newton_loops = 1L),
      )
    })
    sanity_all <- all(unlist(sanity(fit_all)))

    # fit_up <- try({sdmTMB(
    #   formula = response ~ 0 + as.factor(year),
    #   data = surv_dat_up,
    #   family = mi$family,
    #   time = "year",
    #   spatiotemporal = mi$spatiotemporal,
    #   offset = "offset",
    #   mesh = mesh_up,
    #   anisotropy = FALSE,
    #   priors = priors,
    #   silent = SILENT,
    #   control = sdmTMBcontrol(newton_loops = 1L),
    # )})
    # sanity_up <- all(unlist(sanity(fit_up)))

    fit_down <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat_down,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_down,
        anisotropy = FALSE,
        priors = priors,
        silent = SILENT,
        control = sdmTMBcontrol(newton_loops = 1L),
      )
    })
    sanity_down <- all(unlist(sanity(fit_down)))

    fit_down2 <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = surv_dat_down2,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_down2,
        anisotropy = FALSE,
        priors = priors,
        silent = SILENT,
        control = sdmTMBcontrol(newton_loops = 1L),
      )
    })
    sanity_down2 <- all(unlist(sanity(fit_down2)))

    if (!sanity_all || !sanity_restr) {
      saveRDS(NULL, paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
      return(NULL)
    }

    # up-sample -------------

    # - figure out how many points have been restricted each year/stratum
    # - sample an equivalent number without replacement from the remaining blocks each year
    # - if any are within 2 km of an existing sample, try again!?
    # - this is a quick hack to avoid and assigning survey blocks to all existing
    #   survey at locations, which is tricky because sometimes they fall very
    #   slightly outside of the block
    # number_restricted <- surv_dat |>
    #   group_by(survey_abbrev, year) |>
    #   summarise(n = sum(restricted))
    #
    # available_grid <- filter(grid, !restricted)
    # unavailable_grid <- filter(grid, restricted)

    up_sampled_dat <- surv_dat |>
      group_by(survey_abbrev, year, grouping_code) |>
      group_split() |>
      purrr::map_dfr(function(xx) {
        number_restricted <- sum(xx$restricted)
        available_grid <- grid |> filter(grouping_code %in% xx$grouping_code, !restricted)
        replace <- if (number_restricted > nrow(available_grid)) TRUE else FALSE
        out <- available_grid |> sample_n(size = number_restricted, replace = replace)
        out <- mutate(out, year = unique(xx$year))
        # print("stratum")
        # print(unique(xx$grouping_code))
        # print(out)
        out <- mutate(out, up_sample = TRUE)
        out <- mutate(out, offset = 0)
        bind_rows(out, filter(xx, !restricted) |> mutate(up_sample = FALSE))
      })
    # glimpse(up_sampled_dat)

    # pred <- predict(fit_all, newdata = up_sampled_dat, type = "response")
    # if (isTRUE(mi$family$delta)) {
    #   s1 <- rbinom(nrow(pred), size = 1L, prob = pred$est1)
    #   shape <- exp(fit_all$model$par[["ln_phi"]])
    #   scale <- pred$est2 / shape
    #   s2 <- rgamma(nrow(pred), shape = shape, scale = scale)
    #   sim <- s1 * s2
    # } else {
    #   p <- stats::plogis(fit_all$model$par[["thetaf"]]) + 1
    #   dispersion <- exp(fit_all$model$par[["ln_phi"]])
    #   sim <- fishMod::rTweedie(n = nrow(pred), p = p, mu = pred$est, phi = dispersion)
    # }

    if (isTRUE(mi$family$delta)) {
      pred1 <- predict(fit_all, newdata = up_sampled_dat, model = 1L, type = "response", nsim = 1L)
      pred2 <- predict(fit_all, newdata = up_sampled_dat, model = 2L, type = "response", nsim = 1L)
      s1 <- rbinom(nrow(pred1), size = 1L, prob = pred1)
      shape <- exp(fit_all$model$par[["ln_phi"]])
      scale <- pred2 / shape
      s2 <- rgamma(nrow(pred2), shape = shape, scale = scale)
      sim <- s1 * s2
    } else {
      pred <- predict(fit_all, newdata = up_sampled_dat, type = "response", nsim = 1L)
      # new_par <- sdmTMB:::rmvnorm_prec(fit_all$tmb_obj$env$last.par.best, fit_all$sd_report, 1L)
      # phi_i <- names(fit_all$tmb_obj$env$last.par.best) == "ln_phi"
      # phi <- exp(new_par[phi_i])
      # theta_i <- names(fit_all$tmb_obj$env$last.par.best) == "thetaf"
      # phi <- exp(new_par[phi_i])
      # p <- plogis(new_par[theta_i]) + 1
      # dispersion <- exp(new_par[phi_i])

      p <- stats::plogis(fit_all$model$par[["thetaf"]]) + 1
      dispersion <- exp(fit_all$model$par[["ln_phi"]])
      sim <- fishMod::rTweedie(n = nrow(pred), p = p, mu = pred, phi = dispersion)
    }
    up_sampled_dat$simulated_response <- sim

    filter(up_sampled_dat, !up_sample) |>
      ggplot(aes(response + 1, simulated_response + 1)) +
      geom_point() +
      scale_x_log10() +
      scale_y_log10()
    group_by(up_sampled_dat, up_sample, year) |>
      summarise(v = var(response/exp(offset)), vsim = var(simulated_response/exp(offset)))
    group_by(up_sampled_dat, up_sample, year) |>
      summarise(v = var(response/exp(offset)), vsim = var(simulated_response/exp(offset))) |>
      summarise(mean_v = mean(v), mean_vsim = mean(vsim))

    up_sampled_dat$response <- ifelse(up_sampled_dat$up_sample, up_sampled_dat$simulated_response, up_sampled_dat$response)
    stopifnot(sum(is.na(up_sampled_dat$response)) == 0L)

    mesh_up <- make_mesh(up_sampled_dat, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)
    fit_up <- try({
      sdmTMB(
        formula = response ~ 0 + as.factor(year),
        data = up_sampled_dat,
        family = mi$family,
        time = "year",
        spatiotemporal = mi$spatiotemporal,
        offset = "offset",
        mesh = mesh_up,
        anisotropy = FALSE,
        priors = priors,
        silent = SILENT,
        control = sdmTMBcontrol(newton_loops = 1L),
      )
    })
    sanity_up <- all(unlist(sanity(fit_up)))

    # xx <- filter(surv_dat, year == 2005, survey_abbrev == "SYN HS")
    # number_restricted <- sum(xx$restricted)
    #
    # x <- sample_n(available_grid, size = number_restricted, replace = FALSE)
    # surv_dat
    #
    # plot(xx$X, xx$Y, col = "black")
    # points(x$X, x$Y, col = "red")
    # points(nox$X, nox$Y, col = "#0000FF10", pch = 4)
    #

    # x <- sample_n(available_grid, size = 2000, replace = TRUE)
    # nox <- sample_n(unavailable_grid, size = 2000, replace = TRUE)
    # plot(xx$X, xx$Y, col = "black")
    # points(x$X, x$Y, col = "red")
    # points(nox$X, nox$Y, col = "#0000FF10", pch = 4)

    # ----------------------

    grid <- sdmTMB::replicate_df(grid, "year", time_values = unique(surv_dat$year))

    gr_full <- dplyr::filter(grid, survey_abbrev %in% survey)
    gr_mpa <- dplyr::filter(grid, survey_abbrev %in% survey, restricted == TRUE)
    gr_remaining <- dplyr::filter(grid, !restricted, survey_abbrev %in% survey)

    p <- predict(fit_all, newdata = gr_full, return_tmb_object = TRUE)
    ind <- get_index(p, bias_correct = TRUE)

    pmpa <- predict(fit_all, newdata = gr_mpa, return_tmb_object = TRUE)
    indmpa <- get_index(pmpa, bias_correct = TRUE)

    pmpa_restr <- predict(fit_restr, newdata = gr_mpa, return_tmb_object = TRUE)
    indmpa_restr <- get_index(pmpa_restr, bias_correct = TRUE)

    pr <- predict(fit_restr, newdata = gr_full, return_tmb_object = TRUE)
    indr <- get_index(pr, bias_correct = TRUE)

    prs <- predict(fit_restr, newdata = gr_remaining, return_tmb_object = TRUE)
    indrs <- get_index(prs, bias_correct = TRUE)

    if (sanity_down) {
      p_down <- predict(fit_down, newdata = gr_remaining, return_tmb_object = TRUE)
      ind_down <- get_index(p_down, bias_correct = TRUE)

      p_down <- predict(fit_down, newdata = gr_full, return_tmb_object = TRUE)
      ind_down_full <- get_index(p_down, bias_correct = TRUE)
    }
    if (sanity_down2) {
      p_down2 <- predict(fit_down2, newdata = gr_remaining, return_tmb_object = TRUE)
      ind_down2 <- get_index(p_down2, bias_correct = TRUE)

      p_down2 <- predict(fit_down2, newdata = gr_full, return_tmb_object = TRUE)
      ind_down2_full <- get_index(p_down2, bias_correct = TRUE)
    }
    if (sanity_up) {
      p_up <- predict(fit_up, newdata = gr_remaining, return_tmb_object = TRUE)
      ind_up <- get_index(p_up, bias_correct = TRUE)
    }

    # p_up <- predict(fit_up, newdata = gr_remaining, return_tmb_object = TRUE)
    # ind_up <- get_index(p_up, bias_correct = TRUE)

    i <- bind_rows(
      mutate(ind, type = "Status quo"),
      mutate(indr, type = "Restricted"),
      mutate(indrs, type = "Restricted and shrunk"),
      mutate(indmpa, type = "MPA only"),
      # mutate(ind_up, type = "Restricted, shrunk, up-sampled"),
      mutate(indmpa_restr, type = "MPA only restricted")
    )

    if (sanity_down) {
      i <- bind_rows(
        i,
        mutate(ind_down, type = "Random down-sampled and shrunk")
      )
      i <- bind_rows(
        i,
        mutate(ind_down_full, type = "Random down-sampled")
      )
    }
    if (sanity_down2) {
      i <- bind_rows(
        i,
        mutate(ind_down2, type = "Random down-sampled and shrunk 2")
      )
      i <- bind_rows(
        i,
        mutate(ind_down2_full, type = "Random down-sampled 2")
      )
    }
    if (sanity_up) {
      i <- bind_rows(
        i,
        mutate(ind_up, type = "Random up-sampled and shrunk")
      )
    }

    i$species_common_name <- spp
    i$survey_abbrev <- paste(survey, collapse = ", ")
    i$family <- paste(as.character(mi$family$family), collapse = "-")
    i$spatiotemporal <- paste(as.character(mi$spatiotemporal), collapse = ", ")
    saveRDS(i, paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
  } else {
    i <- readRDS(paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
  }

  if (!is.null(i)) {
    g <- i |>
      filter(!grepl("MPA", type)) |>
      # filter(!grepl("Random", type)) |>
      ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type)) +
      geom_pointrange(position = position_dodge(width = 1), pch = 21) +
      # scale_colour_manual(values =
      #     c("Restricted" = "red", "Status quo" = "grey60", "Restricted and shrunk" = "purple")) +
      scale_colour_brewer(palette = "Dark2") +
      ylab("Index") +
      xlab("Year") +
      labs(colour = "Type") +
      ggtitle(spp)

    g1 <- g + coord_cartesian(expand = FALSE, ylim = c(0, NA))
    g2 <- g + scale_y_log10() + ylab("Index (log distributed)")
    g0 <- cowplot::plot_grid(g1, g2, nrow = 2)

    ggsave(paste0("figs/indexes/", spp_file, "-", surv_file, ".pdf"),
      width = 10, height = 7
    )
  }

  return(NULL)
}

source("analysis/spp.R")

syn_survs <- c("SYN WCHG", "SYN QCS|SYN HS", "SYN QCS", "SYN HS")

library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
cores <- round(parallel::detectCores() / 2) + 2
if (!is_rstudio && is_unix) plan(multicore, workers = cores) else plan(multisession, workers = cores)

to_fit <- expand_grid(spp = syn_highlights, survey = syn_survs)

## test:
calc_indices(spp = syn_highlights[17], survey = syn_survs[4], force = TRUE)
x <- readRDS("data-generated/indexes/arrowtooth-flounder-SYN-WCHG.rds")
x |>
  filter(type != "MPA only", type != "MPA only restricted", type != "Restricted") |>
  ggplot(aes(year, se, colour = type)) +
  geom_point() +
  geom_line()

# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices)

to_fit <- expand_grid(spp = hbll_highlights, survey = "HBLL OUT N")
## test:
# calc_indices(spp = hbll_highlights[1], survey = "HBLL OUT N")
# x <- readRDS("data-generated/indexes/arrowtooth-flounder-HBLL-OUT-N.rds")
# x |> filter(type != "MPA only", type != "MPA only restricted", type != "Restricted") |>
#   ggplot(aes(year, se, colour = type)) + geom_point() + geom_line()

# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices)

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
