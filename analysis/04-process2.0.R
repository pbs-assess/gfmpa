library(dplyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
hbll$est_type <- "geostat"
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
syn$est_type <- "geostat"
design <- readRDS(file = "data-generated/stratified-random-design-all.rds")
s <- bind_rows(list(hbll, syn, design))

# nrow(s)
# table(s$type)
s <- filter(s, type != "MPA only restricted")
# table(s$est_type)
# table(s$family)
# table(s$spatiotemporal)
# table(s$survey_abbrev)
# table(s$survey_abbrev)

s$species_common_name <- stringr::str_to_title(s$species_common_name)
s <- mutate(s, species_common_name = gsub(
  "Rougheye/Blackspotted Rockfish Complex",
  "Rougheye/Blackspotted Rockfish", species_common_name
))

filter(s, species_common_name == "Aleutian Skate", survey_abbrev == "SYN WCHG")

prop_mpa <- s |>
  filter(type %in% c("Status quo", "MPA only"), est_type == "geostat") |>
  group_by(species_common_name, survey_abbrev) |>
  summarise(prop_mpa = mean(est[type == "MPA only"] / est[type == "Status quo"]))
stopifnot(sum(is.na(prop_mpa$prop_mpa)) == 0L)

squo_est <- s |>
  filter(type %in% "Status quo") |>
  select(species_common_name, survey_abbrev, year, est_type, orig_est = est) |>
  distinct()

s <- left_join(s, prop_mpa) |>
  left_join(squo_est)

# geo-mean center
s <- s |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  mutate(lwr = lwr / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(upr = upr / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(est = est / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(orig_est = orig_est / exp(mean(log(orig_est), na.rm = TRUE))) |>
  ungroup()

s <- s |> mutate(re = (est - orig_est) / orig_est)

re_slopes <- s |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  group_split() |>
  purrr::map_dfr(function(x) {
    x$decade <- x$year / 10
    fit <- try({lm(re ~ decade, data = x)}, silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(data.frame(slope_re = NA_real_, se_slope_re = NA_real_,
        species_common_name = x$species_common_name[1], survey_abbrev = x$survey_abbrev[1],
        type = x$type[1], est_type = x$est_type[1]))
    } else {
      return(data.frame(slope_re = coef(fit)[[2]], se_slope_re = summary(fit)$coefficients[2, 2],
        species_common_name = x$species_common_name[1], survey_abbrev = x$survey_abbrev[1],
        type = x$type[1], est_type = x$est_type[1]))
    }
  })
re_slopes <- mutate(re_slopes,
  slope_re_med = slope_re,
  slope_re_lwr = slope_re - 1.96 * se_slope_re,
  slope_re_upr = slope_re + 1.96 * se_slope_re
  ) |>
  select(-slope_re, -se_slope_re)

metrics <- s |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  summarise(
    cv_med = median(cv, na.rm = TRUE),
    cv_lwr = quantile(cv, probs = 0.1, na.rm = TRUE),
    cv_upr = quantile(cv, probs = 0.9, na.rm = TRUE),
    mare_med = median(abs(re), na.rm = TRUE),
    mare_lwr = quantile(abs(re), probs = 0.1, na.rm = TRUE),
    mare_upr = quantile(abs(re), probs = 0.9, na.rm = TRUE),
    cv_perc_med = mean(((cv - orig_cv) / orig_cv) * 100, na.rm = TRUE),
    cv_perc_lwr = quantile(((cv - orig_cv) / orig_cv) * 100, probs = 0.1, na.rm = TRUE),
    cv_perc_upr = quantile(((cv - orig_cv) / orig_cv) * 100, probs = 0.9, na.rm = TRUE)
  ) |>
  ungroup()

metrics <- left_join(metrics, re_slopes)


