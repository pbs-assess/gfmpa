library(dplyr)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)

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
    cv_perc_upr = quantile(((cv - orig_cv) / orig_cv) * 100, probs = 0.9, na.rm = TRUE),
    coverage = mean(orig_est < upr_50 & orig_est > lwr_50)
  ) |>
  ungroup()

metrics <- left_join(metrics, re_slopes) |>
  mutate(abs_slope_re = abs(slope_re_med))

glimpse(metrics)

metrics |> group_by(est_type, type) |>
  filter(type != "MPA only") |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  summarise(
    mean_cv = mean(cv_med, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_cv_change = mean(cv_perc_med, na.rm = TRUE),
    mean_mare = mean(mare_med, na.rm = TRUE),
    mean_slope = mean(abs_slope_re, na.rm = TRUE)
  ) |>
  # filter(survey_abbrev == "SYN WCHG") |>
  # filter(survey_abbrev == "HBLL OUT N") |>
  # filter(!grepl("down", type)) |>
  filter(!grepl("cochran", est_type)) |>
  filter(type %in% c("Restricted and shrunk", "Random up-sampled and shrunk", "Restricted and shrunk", "Status quo", "Random down-sampled and shrunk", "Random down-sampled and shrunk 2")) |>
  # arrange(-mean_slope) |>
  arrange(mean_cv_change) |>
  # arrange(mean_cv) |>
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

metrics |> group_by(est_type, type) |>
  summarise(mean_cv_med = mean(mare_med)) |>
  arrange(mean_cv_med)

# remove bottom 10% quantile of biomass
# plot omega + intercept with 'v' slopes zeta things
# as a check - plot prediction at beginning/middle/end
# all juvenile, mature male, mature female

# is 5% and 95% good maturity thresholds!?
# include maturity!? which is middle 90%
# potentialy add pre 2003 lengths to density plot for US West Coast
