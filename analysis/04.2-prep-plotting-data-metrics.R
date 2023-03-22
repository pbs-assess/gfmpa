library(dplyr)

index <- readRDS("data-generated/index-filtered.rds")
cvdata1 <- readRDS("data-generated/hbll-cv-w-lm-slopes.rds")
cvdata2 <- readRDS("data-generated/syn-cv-w-lm-slopes.rds")

lu <- tibble(
  "Restriction type" = c("re_restr", "re_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)
lu_cv2 <- tibble(
  "Restriction type" = c("cv_change_restr", "cv_change_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

cv2 <- index |>
  group_by(species_common_name, survey_abbrev, year) |>
  summarise(
    cv_change_restr = round((cv[type == "Restricted"] - cv[type == "Status quo"]) /
        cv[type == "Status quo"], 8),
    cv_change_shrunk = round((cv[type == "Restricted and shrunk"] - cv[type == "Status quo"]) /
        cv[type == "Status quo"], 8),
    .groups = "drop_last"
  )

cv_long2 <- cv2 |>
  tidyr::pivot_longer(starts_with("cv"),
    names_to = "Restriction type", values_to = "CV change"
  )

x <- index |>
  group_by(species_common_name, survey_abbrev, type) |>
  mutate(est = est / exp(mean(log(est)))) |>
  group_by(species_common_name, survey_abbrev, year) |>
  summarise(
    re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
    re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"],
    .groups = "drop_last"
  )

x_long <- x |>
  tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") |>
  left_join(lu, by = join_by(`Restriction type`))

dd1b <- cv_long2 |>
  group_by(survey_abbrev, species_common_name, `Restriction type`) |>
  summarise(
    lwr = quantile(`CV change` * 100, 0.025),
    upr = quantile(`CV change` * 100, 0.975),
    est = median(`CV change` * 100),
    .groups = "drop_last"
  ) |>
  ungroup() |>
  group_by(species_common_name) |>
  # filter(!survey_abbrev %in% c("SYN HS", "SYN WCHG")) |> # order based only on illustrated surveys?
  mutate(
    est_avg = mean(est, na.rm = TRUE),
    measure = "% increase CV\n(precision loss)"
  ) |>
  left_join(lu_cv2, by = join_by(`Restriction type`))

dd2 <- x_long |>
  group_by(survey_abbrev, species_common_name, `Restriction type`) |>
  summarise(
    lwr = quantile(abs(re), 0.025),
    upr = ifelse(quantile(abs(re), 0.975) > 0.5, 0.5, quantile(abs(re), 0.975)),
    est = median(abs(re)),
    .groups = "drop_last"
  ) |>
  # filter(!survey_abbrev %in% c("SYN HS", "SYN WCHG")) |> # order based only on illustrated surveys?
  mutate(
    est_avg = mean(est, na.rm = TRUE),
    measure = "MARE\n(accuracy loss)"
  )

# combine precision, accuracy and bias data from all surveys ----
cvdata1 <- filter(cvdata1, !is.na(restr_clean))
cvdata2 <- filter(cvdata2, !is.na(restr_clean))

cvdata <- bind_rows(cvdata1, cvdata2)
cvdata <- mutate(cvdata, species_common_name = gsub("Rougheye/Blackspotted Rockfish Complex", "Rougheye/Blackspotted Rockfish", species_common_name))

# exclude any that were filtered on orig CV > 1
incl <- select(index, species_common_name, survey_abbrev) |> distinct()
cvdata <- semi_join(cvdata, incl)

saveRDS(cvdata, file = "data-generated/metrics-wide.rds")

dd3 <- cvdata |>
  group_by(survey_abbrev, species_common_name, `Restriction type`) |>
  summarise(
    lwr = (median(slope_re) - qnorm(0.975) * mean(se_slope_re)) / 1,
    upr = (median(slope_re) + qnorm(0.975) * mean(se_slope_re)) / 1,
    est = median((slope_re)) / 1,
    .groups = "drop_last"
  ) |>
  # filter(!survey_abbrev %in% c("SYN HS", "SYN WCHG")) |> # order based only on illustrated surveys?
  mutate(
    est_avg = mean(abs(est), na.rm = TRUE) / 1,
    measure = "RE trend\n(trend bias)"
  ) |>
  select(
    survey_abbrev, species_common_name, `Restriction type`,
    lwr, upr,
    est_avg,
    est, measure
  )

dd <- bind_rows(dd2, dd3) |>
  left_join(lu, by = join_by(`Restriction type`)) |>
  bind_rows(dd1b)

dd$survey_abbrev <- factor(dd$survey_abbrev, levels = c("SYN QCS, SYN HS", "HBLL OUT N", "SYN WCHG" ))

saveRDS(dd, file = "data-generated/metrics-long.rds")
