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

calc_indices <- function(spp, survey) {
  cat(spp, "\n")
  cat(survey, "\n")
  spp_file <- gsub(" ", "-", spp)
  spp_file <- gsub("\\/", "-", spp_file)
  surv_file <- gsub(" ", "-", survey)
  if (file.exists(paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))) {
    return(NULL)
  }
  ggplot2::theme_set(ggplot2::theme_light())

  if (grepl("SYN", survey)) {
    surv_dat <- readRDS("data-generated/dat_to_fit.rds")
    grid <- readRDS("data-generated/syn-grid-w-restr.rds")
    grid$area <- 4
  } else {
    surv_dat <- readRDS("data-generated/dat_to_fit_hbll.rds")
    grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
    grid$survey_abbrev <- "HBLL OUT N"
    grid$area <- 4
  }

  surv_dat <- filter(surv_dat, survey_abbrev == survey, species_common_name == spp)
  if (nrow(surv_dat) == 0L) return(NULL)

  if (grepl("SYN", survey)) {
    surv_dat$area_swept1 <- surv_dat$doorspread_m * surv_dat$tow_length_m
    surv_dat$area_swept2 <- surv_dat$doorspread_m * surv_dat$duration_min * surv_dat$speed_mpm
    surv_dat$area_swept <- ifelse(!is.na(surv_dat$tow_length_m),
      surv_dat$area_swept1, surv_dat$area_swept2
    )
    surv_dat <- dplyr::filter(surv_dat, !is.na(area_swept))
    surv_dat$offset <- log(surv_dat$area_swept * 0.00001)
    surv_dat$response <- surv_dat$catch_weight
  } else {
    surv_dat$offset <- log(surv_dat$hook_count)
    surv_dat$response <- surv_dat$catch_count
  }

  g <- ggplot(surv_dat, aes(X, Y, colour = restricted, size = response / exp(offset))) +
    geom_point(pch = 21) +
    facet_wrap(~year) +
    coord_fixed() +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "grey60")) +
    ggtitle(spp)
  ggsave(paste0("figs/raw-data-maps/", spp_file, "-", surv_file, ".pdf"), width = 10, height = 10)

  priors <- sdmTMBpriors(
    matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
    matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
  )

  surv_dat_r <- filter(surv_dat, !restricted)

  mesh_all <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = 8)
  mesh_restr <- make_mesh(surv_dat_r, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)

  mi <- list(
    spatiotemporal = list("iid", "iid"),
    family = sdmTMB::delta_gamma()
  )

  fit_restr <- try({sdmTMB(
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
  )})
  s <- sanity(fit_restr)

  if (!s$all_ok) {
    mi$spatiotemporal <- list("off", "iid")
    fit_restr <- try({update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)})
    s <- sanity(fit_restr)
  }
  if (!s$all_ok) {
    mi$family <- sdmTMB::tweedie()
    mi$spatiotemporal <- "iid"
    fit_restr <- try({update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)})
    s <- sanity(fit_restr)
  }
  if (!s$all_ok) {
    mi$family <- sdmTMB::delta_gamma()
    mi$spatiotemporal <- list("off", "off")
    fit_restr <- try({update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)})
    s <- sanity(fit_restr)
  }
  if (!s$all_ok) {
    mi$spatiotemporal <- "off"
    mi$family <- sdmTMB::tweedie()
    fit_restr <- try({update(fit_restr, spatiotemporal = mi$spatiotemporal, family = mi$family)})
  }
  sanity_restr <- sanity(fit_restr)

  fit_all <- try({sdmTMB(
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
  )})
  sanity(fit_all)
  sanity_all <- sanity(fit_all)

  if (!sanity_all$all_ok || !sanity_restr$all_ok) {
    return(NULL)
  }

  fit_all
  fit_restr

  gr <- dplyr::filter(grid, year %in% surv_dat$year, survey_abbrev == survey)
  grs <- dplyr::filter(grid, year %in% surv_dat$year, !restricted, survey_abbrev == survey)

  p <- predict(fit_all, newdata = gr, return_tmb_object = TRUE)
  ind <- get_index(p, bias_correct = TRUE)

  pr <- predict(fit_restr, newdata = gr, return_tmb_object = TRUE)
  indr <- get_index(pr, bias_correct = TRUE)

  prs <- predict(fit_restr, newdata = grs, return_tmb_object = TRUE)
  indrs <- get_index(prs, bias_correct = TRUE)

  i <- bind_rows(
    mutate(ind, type = "Status quo"),
    mutate(indr, type = "Restricted"),
    mutate(indrs, type = "Restricted and shrunk")
  )

  g <- ggplot(i, aes(year, est, ymin = lwr, ymax = upr, colour = type)) +
    geom_pointrange(position = position_dodge(width = 0.55), pch = 21) +
    scale_colour_manual(values =
        c("Restricted" = "red", "Status quo" = "grey60", "Restricted and shrunk" = "purple")) +
    ylab("Index") + xlab("Year") +
    labs(colour = "Type") +
    ggtitle(spp)

  g1 <- g + coord_cartesian(expand = FALSE, ylim = c(0.01, NA))
  g2 <- g + scale_y_log10() + ylab("Index (log distributed)")
  g0 <- cowplot::plot_grid(g1, g2, nrow = 2)

  ggsave(paste0("figs/indexes/", spp_file, "-", surv_file, ".pdf"), width = 7, height = 7)

  i$species_common_name <- spp
  i$survey_abbrev <- survey
  i$cv <- mutate(i, cv = sqrt(exp(se^2) - 1))

  saveRDS(i, paste0("data-generated/indexes/", spp_file, "-", surv_file, ".rds"))
}

source("analysis/spp.R")

syn_survs <- c("SYN WCHG", "SYN QCS", "SYN HS")

library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 5L) else plan(multisession, workers = 5L)

to_fit <- expand_grid(spp = syn_highlights, survey = syn_survs)
# calc_indices(spp = syn_highlights[1], survey = syn_survs[1])
# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices)

to_fit <- expand_grid(spp = hbll_highlights, survey = "HBLL OUT N")
# calc_indices(spp = hbll_highlights[1], survey = "HBLL OUT N")
# purrr::pmap(to_fit, calc_indices)
furrr::future_pmap(to_fit, calc_indices)


