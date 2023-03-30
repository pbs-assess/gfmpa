library(ggplot2)
library(dplyr)
library(future)
plan(multisession, workers = 8L)

# calculate design-based biomass estimate from output of get_survey_sets()
calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] |>
    group_by(year, survey_abbrev, area_km2, grouping_code) |>
    summarise(density = mean(density_kgpm2 * 1e6), .groups = "drop_last") |>
    group_by(year) |>
    summarise(
      survey_abbrev = dat$survey_abbrev[1],
      species_common_name = dat$species_common_name[1],
      est = sum(density * area_km2), .groups = "drop"
    )
}

calc_bio_one_year <- function(dat, i = seq_len(nrow(dat))) {
  d <- dat[i, ] |>
    group_by(grouping_code, area_km2) |>
    summarise(density = mean(density_kgpm2), .groups = "drop")
  sum(d$density * d$area_km2)
}

bootstrap_one_year <- function(x, reps) {
  out <- tryCatch({
  b <- boot::boot(x,
    statistic = calc_bio_one_year, strata = x$grouping_code,
    R = reps
  )
  suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
  tibble::tibble(
    survey_abbrev = x$survey_abbrev[1],
    year = x$year[1],
    est = mean(b$t),
    lwr = bci$percent[[4]],
    upr = bci$percent[[5]],
    cv = sd(b$t) / mean(b$t)
  )
  }, error = function(e) {
    tibble::tibble(
      survey_abbrev = x$survey_abbrev[1],
      year = x$year[1],
      est = mean(b$t),
      lwr = NA_real_,
      upr = NA_real_,
      cv = NA_real_
    )
  })
  out
}

boot_biomass_furrr <- function(dat, reps = 500L) {
  set.seed(12830)
  dat |>
    group_split(year) %>%
    furrr::future_map_dfr(bootstrap_one_year,
      reps = reps,
      .options = furrr::furrr_options(seed = TRUE)
    )
}

boot_biomass_purrr <- function(dat, reps = 500L) {
  set.seed(12830)
  dat |>
    group_split(year) %>%
    purrr::map_dfr(bootstrap_one_year,
      reps = reps
    )
}

surv_dat <- readRDS("data-generated/dat_to_fit.rds")
hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
hbll$density_kgpm2 <- hbll$density_ppkm2 / 1e6 # hack; 1e6 gets reversed
surv_dat <- bind_rows(surv_dat, hbll)

areas <- readRDS("data-generated/stratum-areas.rds") |>
  rename(status_quo_area = area, post_nsb_area = area_nsb)
surv_dat <- left_join(surv_dat, areas)
surv_dat$area_km2 <- NULL # filled later

dat_status_quo_list <- surv_dat |>
  rename(area_km2 = status_quo_area) |>
  group_by(species_common_name, survey_abbrev) |>
  group_split()

dat_nsb_list <- surv_dat |>
  rename(area_km2 = post_nsb_area) |>
  filter(!restricted) |>
  group_by(species_common_name, survey_abbrev) |>
  group_split()

# boot_status_quo <- purrr::map_dfr(dat_status_quo_list, function(.x) {
#   cat(.x$species_common_name[1], .x$survey_abbrev[1], "\n")
#   out <- boot_biomass_furrr(.x, reps = 1000L)
#   out$species_common_name <- .x$species_common_name[1]
#   select(out, survey_abbrev, species_common_name, everything())
# })
#
# boot_nsb <- purrr::map_dfr(dat_nsb_list, function(.x) {
#   out <- boot_biomass_furrr(.x, reps = 1000L)
#   out$species_common_name <- .x$species_common_name[1]
#   select(out, survey_abbrev, species_common_name, everything())
# })

# -----------

# boot_biomass_purrr(dat_status_quo_list[[2]], reps = 10L)
#
# xx <- system.time(boot_biomass_purrr(dat_status_quo_list[[1]], reps = 900L))
# cat("Expectation:",
#   round(length(dat_status_quo_list) * xx[[3]] / 60 / 6, 2),
#   "minutes\n"
# )

tictoc::tic()
boot_status_quo <- furrr::future_map_dfr(dat_status_quo_list, function(.x) {
  out <- boot_biomass_purrr(.x, reps = 900L)
  out$species_common_name <- .x$species_common_name[1]
  select(out, survey_abbrev, species_common_name, everything())
}, .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)
tictoc::toc()

boot_status_quo$est_type <- "bootstrap"
boot_status_quo$type <- "Status quo"
saveRDS(boot_status_quo, "data-generated/stratified-random-design-boot.rds")

tictoc::tic()
boot_nsb <- furrr::future_map_dfr(dat_nsb_list, function(.x) {
  out <- boot_biomass_purrr(.x, reps = 900L)
  out$species_common_name <- .x$species_common_name[1]
  select(out, survey_abbrev, species_common_name, everything())
}, .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)
tictoc::toc()
boot_nsb$est_type <- "bootstrap"
boot_nsb$type <- "Restricted and shrunk"
saveRDS(boot_nsb, "data-generated/stratified-random-design-boot-nsb.rds")

# Design-based estimators:

# working through Cochran 1977:
# https://ia801409.us.archive.org/35/items/Cochran1977SamplingTechniques_201703/Cochran_1977_Sampling%20Techniques.pdf
# p. 95
# https://umaine.edu/chenlab/wp-content/uploads/sites/185/2016/08/Xu-et-al-2015.pdf

# testing:
library(survey)
d <- dat_status_quo_list[[2]]
d <- filter(d, year == min(d$year))
d$dens <- d$density_kgpm2 * 1e6
mydesign <- svydesign(id = ~ 1, strata = ~ grouping_code, data = d, fpc = ~ area_km2)
s <- svymean(~ dens, design = mydesign)

tot_area <- distinct(select(d, area_km2)) |> pull(area_km2) |> sum()

d |>
  group_by(year, survey_abbrev, area_km2, grouping_code) |>
  summarise(density = mean(density_kgpm2 * 1e6), .groups = "drop_last") |>
  group_by(year) |>
  summarise(
    est = sum(density * area_km2), .groups = "drop"
  )

s[[1]] * tot_area

svytotal(~dens, design = mydesign)

n_h <- group_by(d, grouping_code) |>
  summarise(n_h = n()) |> pull(n_h)

N_h <- group_by(d, grouping_code) |>
  summarise(area = area_km2[1]) |> pull(area)

var_h <- group_by(d, grouping_code) |>
  summarise(v = var(dens)) |> pull(v)

cochran_var <- function(N_h, n_h, var_h, total = FALSE) {
  N <- sum(N_h)
  e1 <- 1 / N^2
  e2 <- 0
  for (i in seq_along(N_h)) {
    e2 <- e2 + N_h[i] * (N_h[i] - n_h[i]) * (var_h[i] / n_h[i])
  }
  if (total) {
    return(e2)
  } else {
    return(e1 * e2)
  }
}

svymean(~dens, design = mydesign)
sqrt(cochran_var(N_h, n_h, var_h))

svytotal(~dens, design = mydesign)
W_h <- N_h / sum(N_h)
co <- cochran_var(N_h, n_h, var_h)
sqrt(co * sum(N_h)^2)
sqrt(cochran_var(N_h, n_h, var_h, total = TRUE))

run_cochran_var <- function(dat) {
  n_h <- group_by(dat, grouping_code) |>
    summarise(n_h = n()) |> pull(n_h)
  N_h <- group_by(dat, grouping_code) |>
    summarise(area = area_km2[1]) |> pull(area)
  var_h <- group_by(dat, grouping_code) |>
    summarise(v = var(dens)) |> pull(v)
  sqrt(cochran_var(N_h, n_h, var_h, total = TRUE))
}
run_cochran_var(d)

run_strat_stats <- function(dat) {
  dat <- mutate(dat, dens = density_kgpm2 * 1e6)
  .out <- tryCatch({
  mydesign <- survey::svydesign(id = ~ 1,
    strata = ~ grouping_code, data = dat, fpc = ~ area_km2)
  x <- survey::svytotal(~ dens, design = mydesign)
  data.frame(est = x[[1]], se = as.numeric(sqrt(attr(x, "var"))))
  }, error = function(e) {
    data.frame(est = NA_real_, se = NA_real_)
  })
}

out <- surv_dat |>
  rename(area_km2 = status_quo_area) |> #<
  group_by(species_common_name, survey_abbrev, year) |>
  group_split() |>
  furrr::future_map_dfr(~ {
    tibble(run_strat_stats(.x),
      species_common_name = .x$species_common_name[1],
      survey_abbrev = .x$survey_abbrev[1],
      year = .x$year[1])
  }) |>
  mutate(cv = se / est, lwr = est - qnorm(0.975) * se, upr = est + qnorm(0.975)) |>
  select(survey_abbrev, species_common_name, year, everything())
sum(is.na(out$est))

out_nsb <- surv_dat |>
  filter(!restricted) |> #<
  rename(area_km2 = post_nsb_area) |> #<
  group_by(species_common_name, survey_abbrev, year) |>
  group_split() |>
  furrr::future_map_dfr(~ {
    tibble(run_strat_stats(.x),
      species_common_name = .x$species_common_name[1],
      survey_abbrev = .x$survey_abbrev[1],
      year = .x$year[1])
  }) |>
  mutate(cv = se / est, lwr = est - qnorm(0.975) * se, upr = est + qnorm(0.975)) |>
  select(survey_abbrev, species_common_name, year, everything())
sum(is.na(out_nsb$est))

out$est_type <- "cochran"
out_nsb$est_type <- "cochran"
out$type <- "Status quo"
out_nsb$type <- "Restricted and shrunk"
saveRDS(out, "data-generated/stratified-random-design-var.rds")
saveRDS(out_nsb, "data-generated/stratified-random-design-var-nsb.rds")

ind <- bind_rows(list(boot_status_quo, boot_nsb, out, out_nsb))

ocv <- ind |>
  filter(type == "Status quo") |>
  mutate(orig_cv = cv) |>
  select(year, orig_cv, species_common_name, survey_abbrev, est_type) |>
  distinct()

ind <- left_join(
  ind,
  ocv,
  by = join_by(year, species_common_name, survey_abbrev, est_type)
)

saveRDS(ind, "data-generated/stratified-random-design-all.rds")

plan(sequential)
