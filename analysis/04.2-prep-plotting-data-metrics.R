library(dplyr)

index <- readRDS("data-generated/index-filtered.rds")
index$est_type <- "geostat"
index_design <- readRDS("data-generated/stratified-random-design-all.rds")
index <- bind_rows(index, index_design)

cvdata1 <- readRDS("data-generated/hbll-cv-w-lm-slopes.rds")
cvdata1$est_type <- "geostat"
cvdata2 <- readRDS("data-generated/syn-cv-w-lm-slopes.rds")
cvdata2$est_type <- "geostat"
cvdata_boot <- readRDS("data-generated/design-boot-cv-w-lm-slopes.rds")
cvdata_cochran <- readRDS("data-generated/design-cochran-cv-w-lm-slopes.rds")

lu <- tibble(
  "Restriction type" = c("re_restr", "re_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)
lu_cv2 <- tibble(
  "Restriction type" = c("cv_change_restr", "cv_change_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

cv2 <- index |>
  group_by(species_common_name, survey_abbrev, year, est_type) |>
  summarise(
    cv_change_restr = if ("Restricted" %in% type) round((cv[type == "Restricted"] - cv[type == "Status quo"]) /
        cv[type == "Status quo"], 8) else NA_real_,
    cv_change_shrunk = round((cv[type == "Restricted and shrunk"] - cv[type == "Status quo"]) /
        cv[type == "Status quo"], 8),
    .groups = "drop_last"
  )

cv_long2 <- cv2 |>
  tidyr::pivot_longer(starts_with("cv"),
    names_to = "Restriction type", values_to = "CV change"
  )

x <- index |>
  group_by(species_common_name, survey_abbrev, est_type, type) |>
  mutate(est = est / exp(mean(log(est)))) |>
  group_by(species_common_name, survey_abbrev, est_type, year) |>
  summarise(
    re_restr = if ("Restricted" %in% type) (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"] else NA_real_,
    re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"],
    .groups = "drop_last"
  )

x_long <- x |>
  tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") |>
  left_join(lu, by = join_by(`Restriction type`))

# precision loss -----------------------------------------------------
dd1b <- cv_long2 |>
  group_by(survey_abbrev, species_common_name, est_type, `Restriction type`) |>
  summarise(
    lwr = quantile(`CV change` * 100, 0.025, na.rm = TRUE),
    upr = quantile(`CV change` * 100, 0.975, na.rm = TRUE),
    est = median(`CV change` * 100, na.rm = TRUE),
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

# MARE -----------------------------------------------------
dd2 <- x_long |>
  group_by(survey_abbrev, species_common_name, est_type, `Restriction type`) |>
  summarise(
    lwr = quantile(abs(re), 0.025, na.rm = TRUE),
    upr = ifelse(quantile(abs(re), 0.975, na.rm = TRUE) > 0.9, 0.9, quantile(abs(re), 0.975, na.rm = TRUE)),
    est = median(abs(re), na.rm = TRUE),
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
cvdata_boot <- filter(cvdata_boot, !is.na(restr_clean))
cvdata_cochran <- filter(cvdata_cochran, !is.na(restr_clean))

cvdata <- bind_rows(list(cvdata1, cvdata2, cvdata_boot, cvdata_cochran))
cvdata <- mutate(cvdata, species_common_name = gsub("Rougheye/Blackspotted Rockfish Complex", "Rougheye/Blackspotted Rockfish", species_common_name))

# exclude any that were filtered on orig CV > 1
# incl <- select(index, species_common_name, survey_abbrev, est_type) |> distinct()
# cvdata <- semi_join(cvdata, incl)

saveRDS(cvdata, file = "data-generated/metrics-wide.rds")

# trend in RE ---------------------------------------------------------------
dd3 <- cvdata |>
  group_by(survey_abbrev, species_common_name, est_type, `Restriction type`) |>
  summarise(
    lwr = (median(slope_re) - qnorm(0.975) * mean(se_slope_re)) / 1,
    upr = (median(slope_re) + qnorm(0.975) * mean(se_slope_re)) / 1,
    est = median((slope_re)) / 1,
    .groups = "drop_last"
  ) |>
  mutate(
    est_avg = mean(abs(est), na.rm = TRUE) / 1,
    measure = "RE trend\n(trend bias)"
  ) |>
  select(
    survey_abbrev, species_common_name, est_type, `Restriction type`,
    lwr, upr,
    est_avg,
    est, measure
  )

dd <- bind_rows(dd2, dd3) |>
  left_join(lu, by = join_by(`Restriction type`)) |>
  bind_rows(dd1b)

dd$survey_abbrev <- factor(dd$survey_abbrev, levels = c("SYN QCS, SYN HS", "SYN QCS", "SYN HS", "HBLL OUT N", "SYN WCHG" ))

saveRDS(dd, file = "data-generated/metrics-long.rds")
