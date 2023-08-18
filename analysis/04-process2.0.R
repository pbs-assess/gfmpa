library(dplyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)

survey_data <- readRDS("data-generated/dat_to_fit.rds")
hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
survey_data <- bind_rows(survey_data, hbll)
survey_data <- mutate(survey_data, density = if_else(grepl("HBLL", survey_abbrev), density_ppkm2, density_kgpm2))
prop_mpa_set <- group_by(survey_data, survey_abbrev, species_common_name) |>
  summarise(prop_mpa_set = sum(density[restricted == TRUE], na.rm = TRUE) / sum(density, na.rm = TRUE), .groups = "drop") |>
  mutate(species_common_name = stringr::str_to_title(species_common_name)) |>
  mutate(species_common_name = gsub(
    "Rougheye/Blackspotted Rockfish Complex",
    "Rougheye/Blackspotted Rockfish", species_common_name
  ))

prop_pos <- group_by(survey_data, survey_abbrev, species_common_name) |>
  summarise(prop_pos_set = mean(density > 0), .groups = "drop") |>
  mutate(species_common_name = stringr::str_to_title(species_common_name)) |>
  mutate(species_common_name = gsub(
    "Rougheye/Blackspotted Rockfish Complex",
    "Rougheye/Blackspotted Rockfish", species_common_name
  ))

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
hbll$est_type <- "geostat"
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
syn$est_type <- "geostat"
design <- readRDS(file = "data-generated/stratified-random-design-all.rds")
s <- bind_rows(list(hbll, syn, design))
s <- mutate(s,
  lwr_50 = exp(log(est) - qnorm(0.75) * se),
  upr_50 = exp(log(est) + qnorm(0.75) * se)
)

table(s$type)
table(s$est_type)
s <- filter(s, type != "MPA only restricted")

s$species_common_name <- stringr::str_to_title(s$species_common_name)
s <- mutate(s, species_common_name = gsub(
  "Rougheye/Blackspotted Rockfish Complex",
  "Rougheye/Blackspotted Rockfish", species_common_name
))

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
  left_join(squo_est) |>
  left_join(prop_mpa_set) |>
  left_join(prop_pos)

# geo-mean center
s <- s |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  mutate(lwr_50 = lwr_50 / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(upr_50 = upr_50 / exp(mean(log(est), na.rm = TRUE))) |>
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
    fit <- try(
      {
        lm(re ~ decade, data = x)
      },
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      return(data.frame(
        slope_re = NA_real_, se_slope_re = NA_real_,
        species_common_name = x$species_common_name[1], survey_abbrev = x$survey_abbrev[1],
        type = x$type[1], est_type = x$est_type[1]
      ))
    } else {
      return(data.frame(
        slope_re = coef(fit)[[2]], se_slope_re = summary(fit)$coefficients[2, 2],
        species_common_name = x$species_common_name[1], survey_abbrev = x$survey_abbrev[1],
        type = x$type[1], est_type = x$est_type[1]
      ))
    }
  })
re_slopes <- mutate(re_slopes,
  slope_re_med = slope_re,
  slope_re_med_abs = abs(slope_re),
  slope_re_lwr = slope_re - qnorm(0.75) * se_slope_re,
  slope_re_upr = slope_re + qnorm(0.75) * se_slope_re
) |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  mutate(slope_re_lwr_abs = case_when(
    slope_re_lwr < 0 & slope_re_upr < 0 ~ abs(slope_re_lwr),
    slope_re_lwr > 0 & slope_re_upr > 0 ~ abs(slope_re_lwr),
    slope_re_lwr < 0 & slope_re_upr > 0 ~ 0, # spans 0, take 0
    .default = NA_real_
  )) |>
  mutate(slope_re_upr_abs = case_when(
    slope_re_lwr < 0 & slope_re_upr < 0 ~ abs(slope_re_upr),
    slope_re_lwr > 0 & slope_re_upr > 0 ~ abs(slope_re_upr),
    slope_re_lwr < 0 & slope_re_upr > 0 ~ max(abs(slope_re_upr), abs(slope_re_lwr)), # spans 0, take max abs value
    .default = NA_real_
  )) |>
  select(-slope_re, -se_slope_re) |>
  ungroup()

qs <- c(0.25, 0.75)

metrics <- s |>
  group_by(species_common_name, survey_abbrev, type, est_type) |>
  summarise(
    orig_cv_mean = mean(orig_cv, na.rm = TRUE),
    cv_med = median(cv, na.rm = TRUE),
    cv_lwr = quantile(cv, probs = qs[1], na.rm = TRUE),
    cv_upr = quantile(cv, probs = qs[2], na.rm = TRUE),
    mare_med = median(abs(re), na.rm = TRUE),
    mare_lwr = quantile(abs(re), probs = qs[1], na.rm = TRUE),
    mare_upr = quantile(abs(re), probs = qs[2], na.rm = TRUE),
    cv_perc_med = median(((cv - orig_cv) / orig_cv) * 100, na.rm = TRUE),
    cv_perc_lwr = quantile(((cv - orig_cv) / orig_cv) * 100, probs = qs[1], na.rm = TRUE),
    cv_perc_upr = quantile(((cv - orig_cv) / orig_cv) * 100, probs = qs[2], na.rm = TRUE),
    coverage = mean(orig_est < upr_50 & orig_est > lwr_50)
  ) |>
  ungroup()

metrics <- left_join(metrics, prop_mpa)
metrics <- left_join(metrics, re_slopes)

lwr <- select(metrics, species_common_name:est_type, ends_with("lwr")) |>
  tidyr::pivot_longer(cols = ends_with("lwr"), values_to = "lwr", names_to = "measure") |>
  mutate(measure = gsub("([a-z]+)_lwr", "\\1", measure))
upr <- select(metrics, species_common_name:est_type, ends_with("upr")) |>
  tidyr::pivot_longer(cols = ends_with("upr"), values_to = "upr", names_to = "measure") |>
  mutate(measure = gsub("([a-z]+)_upr", "\\1", measure))
med <- select(metrics, species_common_name:est_type, ends_with("med")) |>
  tidyr::pivot_longer(cols = ends_with("med"), values_to = "est", names_to = "measure") |>
  mutate(measure = gsub("([a-z]+)_med", "\\1", measure))
nrow(lwr)
nrow(upr)
nrow(med)
metrics_long <- left_join(lwr, upr) |>
  left_join(med) |>
  left_join(prop_mpa) |>
  left_join(prop_mpa_set) |>
  left_join(prop_pos) |>
  mutate(measure_clean = case_when(
    measure == "cv_perc" ~ "% increase CV\n(precision loss)",
    measure == "slope_re" ~ "RE trend\n(trend bias)",
    measure == "mare" ~ "MARE\n(accuracy loss)",
    .default = measure
  ))

metrics_long$survey_abbrev <- factor(metrics_long$survey_abbrev,
  levels =
    c("SYN WCHG", "HBLL OUT N", "SYN HS", "SYN QCS", "SYN QCS, SYN HS")
)

saveRDS(metrics_long, "data-generated/metrics-long2.rds")
saveRDS(metrics, "data-generated/metrics-wide2.rds")

metrics |>
  group_by(est_type, type) |>
  filter(type != "MPA only") |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  summarise(
    mean_cv = mean(cv_med, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_cv_change = mean(cv_perc_med, na.rm = TRUE),
    mean_mare = mean(mare_med, na.rm = TRUE),
    mean_slope = mean(slope_re_med, na.rm = TRUE)
  ) |>
  # filter(survey_abbrev == "SYN WCHG") |>
  # filter(survey_abbrev == "HBLL OUT N") |>
  # filter(!grepl("down", type)) |>
  filter(!grepl("cochran", est_type)) |>
  filter(type %in% c(
    "Restricted and shrunk",
    "Random up-sampled and shrunk 1",
    "Random up-sampled and shrunk 2",
    "Random up-sampled and shrunk 3",
    "Random up-sampled and shrunk 4",
    "Random up-sampled and shrunk 5",
    "Restricted and shrunk",
    "Status quo",
    "Random down-sampled and shrunk 1",
    "Random down-sampled and shrunk 2"
  )) |>
  arrange(mean_mare) |>
  knitr::kable(digits = 2)

# Average conclusions
# CV can get worse with upsampling and improved with downsampling!?
# where and why!??
#
# all treatments are more uncertain than status quo
#
# geostat definitely better than design based
#
# within geostat, all pretty similar
#
# all treatments have loss of accuracy
# same or at least need many more seeds to detect diff

# stock-wise percent change in precision, random upsampling does mostly reduce loss of precision

# ON AVERAGE... about 1/4 as much loss of precision...
# could drop extrapolation from up and downsampled...

# add:
# 3-5 down-samples
# 3-5 upsamples-samples

# look at coverage??
# increasingly confident about a wrong index??

# answer:
# 50% coverage does go down with upsampling

metrics |>
  group_by(est_type, type) |>
  summarise(mean_cv_med = mean(mare_med)) |>
  arrange(mean_cv_med)
