setwd(here::here()) # for background RStudio jobs
library(tidyverse)
library(sdmTMB)
# theme_set(theme_light())
source("analysis/theme.R")
library(future)
plan(multisession, workers = 10L)


# survey <- "SYN WCHG"
# spp <- "canary rockfish"
# # spp <- "pacific ocean perch"
# # spp <- "redbanded rockfish"
# # spp <- "redstripe rockfish" # https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2021/2021_014-eng.html
# spp <- "yellowmouth rockfish"
# spp <- "north pacific spiny dogfish"
# spp <- "rougheye/blackspotted rockfish complex"
# spp <- "widow rockfish" # https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41020224.pdf
# spp <- "yellowtail rockfish"
# spp <- "shortraker rockfish"
# spp <- "silvergray rockfish"
# spp <- "lingcod"

metrics_wide <- readRDS("data-generated/metrics-wide2.rds")
m <- metrics_wide |>
  filter(!is.na(prop_mpa)) |>
  filter(prop_mpa > 0.15) |>
  filter(est_type == "geostat") |>
  filter(survey_abbrev %in% c("SYN WCHG", "HBLL OUT N"))
m <- m |>
  # filter(m, !grepl("Greenstripe", species_common_name)) |>
  mutate(species_common_name = tolower(species_common_name)) |>
  select(species_common_name, survey_abbrev, orig_cv_mean) |>
  distinct() |>
  arrange(survey_abbrev, species_common_name)
m$species_common_name[m$species_common_name == "rougheye/blackspotted rockfish"] <- "rougheye/blackspotted rockfish complex"

N_ITER <- 50L

# spp_list <- tolower(sort(unique(m$species_common_name)))
# spp_list <- spp_list[!grepl("greenstripe", spp_list)] # few observations, main model fails
# spp_list[spp_list == "rougheye/blackspotted rockfish"] <- "rougheye/blackspotted rockfish complex"

out_list <- list()
g_list <- list()

for (s_i in seq_len(nrow(m))) {
# for (s_i in 1:2) {
  spp <- m[s_i, "species_common_name", drop = TRUE]
  survey <- m[s_i, "survey_abbrev", drop = TRUE]
  cat(spp, "\n")
  cat(survey, "\n")

  if (grepl("SYN", survey)) {
    surv_dat <- readRDS("data-generated/dat_to_fit.rds")
  } else {
    surv_dat <- readRDS("data-generated/dat_to_fit_hbll.rds")
  }

  # group_by(surv_dat, species_common_name) |> summarise(mpos = mean(frac_pos)) |> arrange(-mpos) |> filter(species_common_name %in% spp_list) |> as.data.frame()

  surv_dat <- dplyr::filter(surv_dat, survey_abbrev %in% survey, species_common_name == spp)

  ########

  # figure out what model was fit before:
  f <- list.files("data-generated/indexes", pattern = ".rds", full.names = TRUE)
  ind <- purrr::map_dfr(f, readRDS)
  ind$cv <- NULL
  ind$cv <- sqrt(exp(ind$se^2) - 1)

  ind_spp <- filter(ind, species_common_name == spp, type == "Status quo", survey_abbrev == survey) |> select(family, anisotropy, spatiotemporal) |> distinct()

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
    # s <- simulate(fit, nsim = 1L)

    p <- select(surv_dat, X, Y, year, restricted)
    p$year_covariate <- p$year - min(p$year)

    suppressMessages({
      b1 <- tidy(fit)
      b2 <- tidy(fit, "ran_pars")
    })
    b <- bind_rows(b1, b2)

    omega_s <- get_pars(fit)$omega_s
    epsilon_st <- get_pars(fit)$epsilon_st

    if (grepl("binomial", ind_spp$family[1])) {
      s_bin <- sdmTMB_simulate(
        formula = ~ 1 + year_covariate,
        data = p,
        mesh = mesh_all,
        family = binomial(),
        time = "year",
        sigma_O = b$estimate[b$term == "sigma_O"],
        # sigma_E = b$estimate[b$term == "sigma_E"],
        sigma_E = if ("sigma_E" %in% b$term) b$estimate[b$term == "sigma_E"] else NULL,
        range = b$estimate[b$term == "range"],
        fixed_re = list(omega_s = omega_s[, 1, drop = FALSE], epsilon_st = epsilon_st[,,1,drop=FALSE], zeta_s = NULL),
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
        fixed_re = list(omega_s = omega_s[, 2, drop = FALSE], epsilon_st = epsilon_st[,,2,drop=FALSE], zeta_s = NULL),
        B = c(mean(b[grep("year", b$term), "estimate", drop = TRUE]), change_per_year),
        seed = 421 * i
      )

      # note offset is not set
      s1 <- select(s_bin, year, X, Y, mu_bin = mu, obs_bin = observed)
      s2 <- select(s_pos, mu_pos = mu, obs_pos = observed)

      s <- bind_cols(s1, s2) |> bind_cols(select(surv_dat, restricted))
      s <- mutate(s, obs = obs_bin * obs_pos, mu = mu_bin * mu_pos)
      # hist(s$obs)
      # head(s)

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
        fixed_re = list(omega_s = omega_s, epsilon_st = epsilon_st, zeta_s = NULL),
        B = c(mean(b[grep("year", b$term), "estimate", drop = TRUE]), change_per_year),
        seed = 421 * i
      )
      s <- select(s, year, X, Y, mu, obs = observed) |> bind_cols(select(surv_dat, restricted))
    }
    sr <- filter(s, !restricted)

    # library(ggplot2)
    #
    # group_by(s_bin, year) |> summarise(m = mean(mu))
    # group_by(s_pos, year) |> summarise(m = mean(mu))
    #
    # ggplot(s, aes(year, mu)) + geom_point()
    # ggplot(s, aes(year, obs)) + geom_point()


    # s$obs_scaled <- s$obs / 1000
    s$year_zero <- s$year - mean(s$year)

    # a <- sdmTMB(obs ~ year_zero, data = s, family = tweedie(link = "log"),
    #   spatial = "on", spatiotemporal = "iid", time = "year", mesh = mesh_all)
    #
    # datr <- filter(s, !restricted)
    # ar <- sdmTMB(obs ~ year_zero, data = datr, family = tweedie(link = "log"),
    #   spatial = "on", spatiotemporal = "iid", time = "year", mesh = meshr)

    fit_squo <- try({
      sdmTMB(list(obs ~ 1, obs ~ year_zero),
        data = s, family = mi$family,
        spatial = "on", spatiotemporal = mi$spatiotemporal, time = "year", mesh = mesh_all
        # priors = priors
      )
    })

    datr <- filter(s, !restricted)
    fit_rest <- try({
      sdmTMB(list(obs ~ 1, obs ~ year_zero),
        data = datr, family = mi$family,
        spatial = "on", spatiotemporal = mi$spatiotemporal, time = "year", mesh = meshr,
        # priors = priors
      )
    })

    if (class(fit_squo) != "try-error" && class(fit_rest) != "try-error") {
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
  }

  # figure out decline to test:
  yrs <- sort(unique(surv_dat$year))
  yrs <- yrs - min(yrs)

  # these are remaining fractions
  # e.g., 0.7 means a 30% decline over 16 years
  # FRAC_TEST <- switch(spp,
  #   "redbanded rockfish" = 0.7,
  #   "pacific ocean perch" = 0.7,
  #   "canary rockfish" = 0.4,
  #   "redstripe rockfish" = 0.6,
  #   "yellowmouth rockfish" = 0.5,
  #   "north pacific spiny dogfish" = 0.6,
  #   "rougheye/blackspotted rockfish complex" = 0.5,
  #   "widow rockfish" = 0.5,
  #   "yellowtail rockfish" = 0.5,
  #   0.5
  # )
  # FRAC_TEST <- 0.5

  cv <- m[s_i, "orig_cv_mean", drop = TRUE]

  FRAC_TEST <- if (cv > 0.5) {
    0.3
  } else if (cv <= 0.5 & cv > 0.3) {
    0.5
  } else {
    0.7
  }

  # FRAC_TEST <- 0.98

  # decl <- 0.98^16
  # decl^(1/16) = 0.98
  # x = log(decl^(1/16))

  # print(spp)
  # print(FRAC_TEST)

  get_change_per_year <- function(x, ny = max(yrs)) {
    log(x^(1 / ny))
  }

  # x <- exp(get_change_per_year(FRAC_TEST) * yrs)
  # plot(yrs, x, type = "o")
  # 1 - x

  out <- furrr::future_map_dfr(
    # out <- purrr::map_dfr(
    seq_len(N_ITER), do_sim_check,
    change_per_year = get_change_per_year(FRAC_TEST),
    .options = furrr::furrr_options(seed = TRUE)
  )

  out$species_common_name <- spp
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
    # geom_pointrange(position = position_dodge(width = 0.5)) +
    geom_pointrange(pch = 21) +
    geom_hline(yintercept = out$true_effect[1], lty = 2) +
    geom_hline(yintercept = 0, lty = 1) +
    coord_flip() +
    scale_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "red")) +
    xlab("Iteration") +
    facet_wrap(~type) +
    ggtitle(spp)

  g_list[[s_i]] <- g
}

plan(sequential)

saveRDS(out_list, file = "data-generated/power-cached-output-2023-06-05.rds")

d <- dplyr::bind_rows(out_list)

d |>
  mutate(sig = conf.high < 0) |>
  ggplot(aes(i, estimate, ymin = conf.low, ymax = conf.high, colour = sig)) +
  # geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_pointrange(pch = 21) +
  geom_hline(yintercept = out$true_effect[1], lty = 2) +
  geom_hline(yintercept = 0, lty = 1) +
  coord_flip() +
  scale_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "red")) +
  xlab("Iteration") +
  facet_wrap(species_common_name ~ type)

cvs <- m |>
  select(species_common_name, orig_cv_mean) |>
  distinct() |>
  mutate(species_common_name = tolower(species_common_name))

x <- group_by(d) |>
  filter(!is.na(conf.high)) |>
  group_by(species_common_name, type) |>
  summarise(
    power = mean(conf.high < 0),
    coverage = mean(conf.high > true_effect & conf.low < true_effect),
    fract_tested = mean(fract_tested)
    # mean_se = mean(std.error),
    # m_error = mean(estimate / true_effect),
    # s_error = mean(estimate > 0),
    # s_error2 = mean(conf.low > 0)
  ) |>
  ungroup() |>
  group_by(species_common_name) |>
  # filter(power[type == "restricted and shrunk"] < 1) |>
  mutate(power_diff = power[type == "restricted and shrunk"] - power[type == "status quo"]) |>
  left_join(cvs) |>
  ungroup() |>
  arrange(fract_tested, power_diff, species_common_name, type) |>
  select(-coverage)

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

select(x, -power_diff) |>
  tidyr::pivot_wider(names_from = type, values_from = power) |>
  mutate(sp = stringr::str_to_title(species_common_name)) |>
  filter(!grepl("Blacksp", sp)) |>
  mutate(sp = gsub("Complex", "", sp)) |>
  mutate(sp = forcats::fct_inorder(sp)) |>
  mutate(sp = forcats::fct_rev(sp)) |>
  ggplot(aes(
    xend = `restricted and shrunk`,
    x = `status quo`,
    y = sp,
    yend = sp,
    colour = as.factor(paste0((1 - fract_tested) * 100, "%"))
  )) +
  scale_colour_viridis_d(end = 0.85, option = "C") +
  geom_segment(arrow = arrow(length = unit(6, "pt")), key_glyph = "arrow_left") +
  # theme_bw() +
  theme(axis.title.y = element_blank()) +
  xlab("Power to detect decline") +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.35, 1)) +
  geom_hline(yintercept = c(14.5, 21.5), lty = 2, col = "grey50") +
  labs(colour = "Simulated decline") +
  theme(
    legend.position = c(0.28, 0.12),
    plot.margin = margin(11 / 2, 11 / 2 + 2, 11 / 2, 11 / 2),
    axis.title = element_text(size = 10)
  )

ggsave("figs/power.png", width = 4.2, height = 4.4)
ggsave("figs/power.pdf", width = 4.2, height = 4.4)
