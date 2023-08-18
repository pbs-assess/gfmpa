setwd(here::here()) # for background RStudio jobs
library(tidyverse)
library(sdmTMB)
# theme_set(theme_light())
source("analysis/theme.R")
library(future)
plan(multisession, workers = 8L)

metrics_wide <- readRDS("data-generated/metrics-wide2.rds")
m <- metrics_wide |>
  filter(!is.na(prop_mpa)) |>
  filter(prop_mpa > 0.15) |>
  filter(est_type == "geostat") |>
  filter(survey_abbrev %in% c("SYN WCHG", "HBLL OUT N"))
m <- m |>
  mutate(species_common_name = tolower(species_common_name)) |>
  select(species_common_name, survey_abbrev, orig_cv_mean) |>
  distinct() |>
  arrange(survey_abbrev, species_common_name)
m$species_common_name[m$species_common_name == "rougheye/blackspotted rockfish"] <- "rougheye/blackspotted rockfish complex"

N_ITER <- 8 * 15

out_list <- list()
g_list <- list()

for (s_i in seq_len(nrow(m))) {
  spp <- m[s_i, "species_common_name", drop = TRUE]
  survey <- m[s_i, "survey_abbrev", drop = TRUE]
  cat(spp, "\n")
  cat(survey, "\n")

  if (grepl("SYN", survey)) {
    surv_dat <- readRDS("data-generated/dat_to_fit.rds")
  } else {
    surv_dat <- readRDS("data-generated/dat_to_fit_hbll.rds")
  }

  surv_dat <- dplyr::filter(surv_dat, survey_abbrev %in% survey, species_common_name == spp)

  # figure out what model was fit before:
  f <- list.files("data-generated/indexes", pattern = ".rds", full.names = TRUE)
  ind <- purrr::map_dfr(f, readRDS)
  ind$cv <- NULL
  ind$cv <- sqrt(exp(ind$se^2) - 1)

  ind_spp <- filter(ind, species_common_name == spp, type == "Status quo", survey_abbrev == survey) |>
    select(family, anisotropy, spatiotemporal) |>
    distinct()

  mi <- list()
  if (ind_spp$family[1] == "tweedie") mi$family <- tweedie()
  delta_model <- FALSE
  if (grepl("binomial", ind_spp$family[1])) {
    mi$family <- delta_gamma()
    delta_model <- TRUE
  }
  st <- as.list(strsplit(ind_spp$spatiotemporal, ", ")[[1]])
  if (length(st) == 1L) st <- st[[1]]
  mi$spatiotemporal <- st

  ########

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

  if (grepl("SYN", survey)) {
    surv_dat <- prep_cols_syn(surv_dat)
  } else {
    surv_dat <- prep_cols_hbll(surv_dat)
  }

  priors <- sdmTMBpriors(
    matern_s = pc_matern(range_gt = 20, sigma_lt = 5),
    matern_st = pc_matern(range_gt = 20, sigma_lt = 5)
  )

  cutoff <- if ("SYN WCHG" %in% survey) 5 else 8
  mesh_all <- make_mesh(surv_dat, xy_cols = c("X", "Y"), cutoff = cutoff)

  fit <-
    sdmTMB(
      response ~ 0 + as.factor(year),
      data = surv_dat,
      family = mi$family,
      time = "year",
      priors = priors,
      spatiotemporal = mi$spatiotemporal,
      offset = "offset",
      mesh = mesh_all
    )
  sanity(fit)

  datr <- filter(surv_dat, !restricted)
  meshr <- make_mesh(datr, xy_cols = c("X", "Y"), mesh = mesh_all$mesh)

  do_sim_check <- function(i, change_per_year = log(0.9)) {
    p <- select(surv_dat, X, Y, year, restricted)
    p$year_covariate <- p$year - min(p$year)

    suppressMessages({
      b1 <- tidy(fit)
      b2 <- tidy(fit, "ran_pars")
    })
    b <- bind_rows(b1, b2)

    omega_s <- get_pars(fit)$omega_s
    epsilon_st <- get_pars(fit)$epsilon_st

    if (delta_model) {
      s_bin <- sdmTMB_simulate(
        formula = ~ 1 + year_covariate,
        data = p,
        mesh = mesh_all,
        family = binomial(),
        time = "year",
        sigma_O = b$estimate[b$term == "sigma_O"],
        sigma_E = if ("sigma_E" %in% b$term) b$estimate[b$term == "sigma_E"] else NULL,
        range = b$estimate[b$term == "range"],
        fixed_re = list(omega_s = omega_s[, 1, drop = FALSE], epsilon_st = NULL, zeta_s = NULL),
        B = c(mean(b[grep("year", b$term), "estimate", drop = TRUE]), 0),
        seed = 42 * i
      )

      suppressMessages({
        b1 <- tidy(fit, model = 2)
        b2 <- tidy(fit, "ran_pars", model = 2)
      })
      b <- bind_rows(b1, b2)

      s_pos <- sdmTMB_simulate(
        formula = ~ 1 + year_covariate,
        data = p,
        mesh = mesh_all,
        family = Gamma(link = "log"),
        time = "year",
        sigma_O = b$estimate[b$term == "sigma_O"],
        sigma_E = if ("sigma_E" %in% b$term) b$estimate[b$term == "sigma_E"] else NULL,
        range = b$estimate[b$term == "range"],
        phi = b$estimate[b$term == "phi"],
        fixed_re = list(omega_s = omega_s[, 2, drop = FALSE], epsilon_st = NULL, zeta_s = NULL),
        B = c(mean(b[grep("year", b$term), "estimate", drop = TRUE]), change_per_year),
        seed = 421 * i
      )

      # note offset is not set
      s1 <- select(s_bin, year, X, Y, mu_bin = mu, obs_bin = observed)
      s2 <- select(s_pos, mu_pos = mu, obs_pos = observed)

      s <- bind_cols(s1, s2) |> bind_cols(select(surv_dat, restricted))
      s <- mutate(s, obs = obs_bin * obs_pos, mu = mu_bin * mu_pos)
    } else { # tweedie
      s <- sdmTMB_simulate(
        formula = ~ 1 + year_covariate,
        data = p,
        mesh = mesh_all,
        family = tweedie(),
        time = "year",
        sigma_O = b$estimate[b$term == "sigma_O"],
        sigma_E = b$estimate[b$term == "sigma_E"],
        tweedie_p = b$estimate[b$term == "tweedie_p"],
        range = b$estimate[b$term == "range"],
        phi = b$estimate[b$term == "phi"],
        fixed_re = list(omega_s = omega_s, epsilon_st = NULL, zeta_s = NULL),
        B = c(mean(b[grep("year", b$term), "estimate", drop = TRUE]), change_per_year),
        seed = 421 * i
      )
      s <- select(s, year, X, Y, mu, obs = observed) |> bind_cols(select(surv_dat, restricted))
    }
    sr <- filter(s, !restricted)

    s$year_zero <- s$year - mean(s$year)

    datr <- filter(s, !restricted)

    sanity_true <- function(x) {
      a <- try({
        sanity(x)
      })
      if (class(a) == "try-error") {
        return(FALSE)
      }
      if (isFALSE(a)) {
        return(FALSE)
      }
      if (length(a) > 1L) {
        a$sigmas_ok <- NULL
        a$se_magnitude_ok <- NULL
        a$all_ok <- NULL
        a$range_ok <- NULL
        return(all(unlist(a)))
      } else {
        return(FALSE)
      }
    }

    if (delta_model) {
      mi_temp <- mi

      fit_squo <- try({
        sdmTMB(list(obs ~ 1, obs ~ year_zero),
          data = s, family = mi_temp$family,
          spatial = "on", spatiotemporal = list("off", "off"), time = "year", mesh = mesh_all,
          priors = priors
        )
      })

      fit_rest <- try({
        sdmTMB(list(obs ~ 1, obs ~ year_zero),
          data = datr, family = mi_temp$family,
          spatial = "on", spatiotemporal = list("off", "off"), time = "year", mesh = meshr,
          priors = priors
        )
      })
    } else { # tweedie
      mi_temp <- mi
      fit_squo <- try({
        sdmTMB(obs ~ year_zero,
          data = s, family = mi_temp$family,
          spatial = "on", spatiotemporal = "off", time = "year", mesh = mesh_all,
          priors = priors
        )
      })

      fit_rest <- try({
        sdmTMB(obs ~ year_zero,
          data = datr, family = mi_temp$family,
          spatial = "on", spatiotemporal = "off", time = "year", mesh = meshr,
          priors = priors
        )
      })
    }

    aq <- group_by(s, year) |>
      summarise(m = mean(mu)) |>
      mutate(type = "squo")
    ar <- group_by(datr, year) |>
      summarise(m = mean(mu)) |>
      mutate(type = "rest")
    bind_rows(aq, ar) |> ggplot(aes(year, m, colour = type)) +
      geom_line()

    if (class(fit_squo) != "try-error" && class(fit_rest) != "try-error") {
      if (sanity_true(fit_rest) && sanity_true(fit_squo)) {
        ret <- bind_rows(
          mutate(tidy(fit_squo, conf.int = TRUE, model = if (delta_model) 2L else 1L), type = "status quo"),
          mutate(tidy(fit_rest, conf.int = TRUE, model = if (delta_model) 2L else 1L), type = "restricted and shrunk")
        ) |>
          filter(grepl("year", term)) |>
          mutate(true_effect = change_per_year, i = i)
        return(ret)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }

  # figure out decline to test:
  yrs <- sort(unique(surv_dat$year))
  yrs <- yrs - min(yrs)

  cv <- m[s_i, "orig_cv_mean", drop = TRUE]

  FRAC_TEST <- if (cv > 0.5) {
    0.5 # same as next...
  } else if (cv <= 0.5 & cv > 0.25) {
    0.5
  } else if (cv <= 0.25 & cv > 0.15) {
    0.7
  } else if (cv <= 0.15 & cv > 0) {
    0.7
  } else {
    stop("???")
  }

  get_change_per_year <- function(x, ny = max(yrs)) {
    log(x^(1 / ny))
  }
  out <- furrr::future_map_dfr(
    seq_len(N_ITER), do_sim_check,
    change_per_year = get_change_per_year(FRAC_TEST),
    .options = furrr::furrr_options(seed = TRUE)
  )
  if (nrow(out) > 0L) {
    out$species_common_name <- spp
    out$survey_abbrev <- survey
    out$fract_tested <- FRAC_TEST
    out_list[[s_i]] <- out
    group_by(out, type) |>
      filter(!is.na(conf.high)) |>
      summarise(
        power = mean(conf.high < 0),
        coverage = mean(conf.high > true_effect & conf.low < true_effect),
        mean_se = mean(std.error),
        m_error = mean(estimate / true_effect),
        s_error = mean(estimate > 0),
        s_error2 = mean(conf.low > 0)
      ) |>
      knitr::kable(digits = 2L)

    g <- out |>
      mutate(sig = conf.high < 0) |>
      ggplot(aes(i, estimate, ymin = conf.low, ymax = conf.high, colour = sig)) +
      geom_pointrange(pch = 21) +
      geom_hline(yintercept = out$true_effect[1], lty = 2) +
      geom_hline(yintercept = 0, lty = 1) +
      coord_flip() +
      scale_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "red")) +
      xlab("Iteration") +
      facet_wrap(~type) +
      ggtitle(spp)

    g_list[[s_i]] <- g
  } else {
    out_list[[s_i]] <- NULL
    g_list[[s_i]] <- NULL
  }
}

plan(sequential)

saveRDS(out_list, file = "data-generated/power-cached-output-2023-06-06-eps-null-spatial.rds")
out_list <- readRDS("data-generated/power-cached-output-2023-06-06-eps-null-spatial.rds")

d <- dplyr::bind_rows(out_list)

cvs <- m |>
  select(species_common_name, survey_abbrev, orig_cv_mean) |>
  distinct() |>
  mutate(species_common_name = tolower(species_common_name))

# https://stackoverflow.com/questions/72741439/how-to-show-arrows-in-backward-and-forward-directions-in-a-ggplot2-legend
draw_key_arrow_left <- function(data, params, size, dir) {
  if (is.null(data$linetype)) {
    data$linetype <- 0
  } else {
    data$linetype[is.na(data$linetype)] <- 0
  }
  grid::segmentsGrob(0.9, 0.5, 0.1, 0.5,
    gp = grid::gpar(
      col = alpha(data$colour %||% data$fill %||% "black", data$alpha),
      fill = alpha(data$colour %||% data$fill %||% "black", data$alpha),
      lwd = (data$size %||% 0.5) * .pt,
      lty = data$linetype %||% 1,
      lineend = "butt"
    ),
    arrow = params$arrow
  )
}

make_power_fig <- function(survs, return_data = FALSE) {
  x <- d |>
    group_by(species_common_name, survey_abbrev) |>
    filter(!is.na(conf.high)) |>
    mutate(n_converged = sum(type == "status quo")) |>
    filter(n_converged > 10) |>
    filter(survey_abbrev %in% survs) |>
    group_by(species_common_name, survey_abbrev, type) |>
    summarise(
      power = mean(conf.high < 0),
      coverage = mean(conf.high > true_effect & conf.low < true_effect),
      fract_tested = mean(fract_tested)
    ) |>
    ungroup() |>
    group_by(species_common_name, survey_abbrev) |>
    mutate(power_diff = power[type == "restricted and shrunk"] - power[type == "status quo"]) |>
    left_join(cvs) |>
    ungroup() |>
    arrange(survey_abbrev, fract_tested, power_diff, species_common_name, type) |>
    select(-coverage)

  x <- select(x, -power_diff) |>
    tidyr::pivot_wider(names_from = type, values_from = power)
  table_dat <<- x |> select(-orig_cv_mean)

  cat("Power less than 0.3 and omitted for space\n")
  under <- filter(x, `status quo` <= 0.2)
  print(paste0(under$species_common_name))
  cat("\n")

  x <- filter(x, `status quo` > 0.2)

  # viridisLite::plasma(3, end = 0.85)
  pal <- c("20%" = "#0D0887FF", "30%" = "#B7318AFF", "50%" = "#FEBA2CFF")
  pal <- c("20%" = "#0D0887FF", "30%" = "#B7318AFF", "50%" = "#0D0887FF")

  x <- x |>
    mutate(sp = stringr::str_to_title(species_common_name)) |>
    mutate(sp = paste0(sp, " (", sdmTMB:::mround(orig_cv_mean, 2), ")")) |>
    filter(!grepl("Blacksp", sp)) |>
    mutate(sp = gsub("Complex", "", sp)) |>
    mutate(sp = forcats::fct_inorder(sp)) |>
    mutate(sp = forcats::fct_rev(sp))
  g <- ggplot(x, aes(
    xend = `restricted and shrunk`,
    x = `status quo`,
    y = sp,
    yend = sp,
    colour = as.factor(paste0((1 - fract_tested) * 100, "%"))
  )) +
    geom_segment(arrow = arrow(length = unit(6, "pt")), key_glyph = "arrow_left") +
    scale_colour_manual(values = pal, guide = guide_legend(reverse = TRUE)) +
    theme(axis.title.y = element_blank()) +
    xlab("Power to detect decline") +
    scale_x_continuous(expand = c(0, 0), limits = c(min(x$`restricted and shrunk`) - 0.02, 1)) +
    labs(colour = "Simulated decline") +
    theme(
      legend.position = c(0.28, 0.17),
      plot.margin = margin(11 / 2, 11 / 2 + 4, 11 / 2, 11 / 2 - 2),
      axis.title = element_text(size = 10)
    )
  if (!return_data) {
    return(g)
  } else {
    return(x)
  }
}

lims <- scale_x_continuous(expand = c(0, 0), limits = c(0.15, 1))
g1 <- make_power_fig("SYN WCHG") +
  geom_hline(yintercept = c(13.5), lty = 2, col = "grey50") +
  lims +
  ggtitle("SYN WCHG") +
  theme(legend.position = "none")

# NUMBERS FOR PAPER:
table_dat |>
  filter(species_common_name %in% c("redstripe rockfish", "redbanded rockfish")) |>
  as.data.frame()
g1
g2 <- make_power_fig("HBLL OUT N") +
  geom_hline(yintercept = c(13.5), lty = 2, col = "grey50") +
  lims +
  ggtitle("HBLL OUT N") +
  theme(legend.position = c(0.29, 0.17))
g2
g <- cowplot::plot_grid(g1, g2, ncol = 1L, align = "v", rel_heights = c(1.37, 1))
g
ggsave("figs/power-june6-eps-null-spatial.png", width = 4.2, height = 6.4)
ggsave("figs/power-june6-eps-null-spatial.pdf", width = 4.2, height = 6.4)

keep <- make_power_fig(c("HBLL OUT N", "SYN WCHG"), return_data = TRUE) |>
  select(species_common_name, survey_abbrev) |>
  distinct()

make_power_fig2 <- function(dd) {
  dd |>
    semi_join(keep) |>
    mutate(sig = conf.high < 0) |>
    ggplot(aes(i, estimate, ymin = conf.low, ymax = conf.high, colour = sig)) +
    geom_pointrange(pch = 21) +
    geom_hline(yintercept = 0, lty = 1) +
    coord_flip() +
    scale_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "red")) +
    xlab("Iteration") +
    facet_wrap(species_common_name ~ type)
}

# dogfish HBLL:
# so yeah, there's a bias based on the sampling locations for dogfish which
# makes the true decline steeper if you ignore the MPA area
# It's not a pure power analysis. It is whether you
# would detect a significant decline if you imposed a given true decline in the
# full survey domain. Sometimes restricting the survey imposes a bias itself
# making the realized trend steeper or shallower.

# I guess, the simplest explanation is that if a population becomes
# increasingly concentrated in the MPAs, obviously we will be prone to
# overestimating its decline. Pretty intuitive actually.

# In outlining the discussion of the GF-NSB MPA project I realized it would be
# useful to do a 'small' power analysis to be able to ascribe meaning to the
# precision effects. Here's where I currently am. I'm running another version
# now. Start of the arrow is the power with the full dataset (probability of a
# continuous year predictor detecting a significant negative slope). Pointy part
# of arrow is the power if removing the NSB areas. Different stocks get
# different simulated declines depending on index CV. The simulation comes from
# the models fitted to the full dataset with estimated random field values and
# observation error.
# It started as a back of the envelope thing for a species or two. It's
# ballooned a bit. But I think it's useful in turning these abstract % increases
# in CV into something concrete.
# Sometimes 2 wrongs make a right: removing the NSB biases the index negative
# making it easier to detect a decline despite fewer samples (dogfish in HBLL is
# a clear one, not shown in this version).
