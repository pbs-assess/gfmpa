# get slopes of relative error
library(dplyr)
library(ggplot2)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())
# options(dplyr.summarise.inform = FALSE)

for (survey in c("HBLL", "SYN", "design-boot", "design-cochran")) {
  # for now using delta-gamma?
  if (survey == "HBLL") index <- readRDS(file = "data-generated/index-hbll-geo-clean.rds")
  if (survey == "SYN") index <- readRDS(file = "data-generated/index-syn-geo-clean.rds")
  if (survey == "design-cochran") {
    index <- readRDS(file = "data-generated/stratified-random-design-all.rds") |> filter(grepl("cochran", est_type))
  }
  if (survey == "design-boot") {
    index <- readRDS(file = "data-generated/stratified-random-design-all.rds") |> filter(grepl("boot", est_type))
  }

  index$species_common_name <- stringr::str_to_title(index$species_common_name)
  index <- mutate(index, species_common_name = gsub(
    "Rougheye/Blackspotted Rockfish Complex",
    "Rougheye/Blackspotted Rockfish", species_common_name
  ))

  index <- index |> dplyr::filter(type != "MPA only restricted")

  has_prop_mpa <- sum(grepl("MPA only", index$type)) > 0
  # get cv ratios and prop mpa
  cv2 <- index |>
    group_by(species_common_name, survey_abbrev, year) |>
    mutate(
      prop_mpa = if (has_prop_mpa) est[type == "MPA only"] / est[type == "Status quo"] else NA_real_,
      cv_ratio_restr = if (has_prop_mpa) cv[type == "Restricted"] /
        cv[type == "Status quo"] else NA_real_,
      cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
        cv[type == "Status quo"],
      cv_shrunk = cv[type == "Restricted and shrunk"]
    )

  cv_long <- cv2 |>
    select(-prop_mpa) |>
    tidyr::pivot_longer(starts_with("cv"),
      names_to = "Restriction type", values_to = "CV ratio"
    )

  prop_mpa <- cv2 |>
    select(species_common_name, survey_abbrev, year, prop_mpa) |>
    group_by(species_common_name, survey_abbrev) |>
    summarise(
      prop_mpa = mean(prop_mpa), .groups = "drop_last"
    )

  lu <- tibble(
    "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  cvratio <- cv_long |>
    group_by(survey_abbrev, species_common_name, `Restriction type`) |>
    summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`), .groups = "drop_last") |>
    left_join(lu)

  # get raw cvs
  lu <- tibble(
    "type" = c("Restricted", "Restricted and shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  cvraw <- index |>
    filter(type != "MPA only") |>
    group_by(species_common_name, survey_abbrev, type) |>
    summarise(
      cv_orig = mean(orig_cv),
      cv_max = max(cv),
      cv_mean = mean(cv),
      .groups = "drop_last"
    ) |>
    # filter(type != "Status quo") |>
    left_join(lu)


  # # bad year in SYN, removed from plot so removing here too
  # if (survey == "SYN") {
  #   index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
  # }

  re <- index |>
    group_by(species_common_name, survey_abbrev, type) |>
    mutate(est = est / exp(mean(log(est)))) |>
    ungroup() |>
    group_by(species_common_name, survey_abbrev) |>
    group_split() |>
    purrr::map_dfr(function(.x) {
      data.frame(
        species_common_name = unique(.x$species_common_name),
        survey_abbrev = unique(.x$survey_abbrev),
        year = .x$year,
        prop_mpa = if (has_prop_mpa) .x$est[.x$type == "Status quo"] else NA_real_,
        re_restr = if (has_prop_mpa) (.x$est[.x$type == "Restricted"] -
            .x$est[.x$type == "Status quo"]) / .x$est[.x$type == "Status quo"] else NA_real_,
        re_shrunk = (.x$est[.x$type == "Restricted and shrunk"] -
            .x$est[.x$type == "Status quo"]) / .x$est[.x$type == "Status quo"]
      )
    })

  lu <- tibble(
    "Restriction type" = c("re_restr", "re_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  re_long <- re |>
    tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") |>
    left_join(lu)

  mare <- re_long |>
    group_by(survey_abbrev, species_common_name, `Restriction type`) |>
    summarise(lwr = min(re), upr = max(re), est = median(abs(re)), .groups = "drop_last") |>
    left_join(lu)

  # rename variables for other previously calculated indices
  cvratio2 <- cvratio |>
    rename(cv_ratio = est, cv_lwr = lwr, cv_upr = upr) |>
    select(-`Restriction type`)
  mare2 <- mare |>
    rename(mare = est, mare_lwr = lwr, mare_upr = upr) |>
    mutate(spp_by_survey = paste(species_common_name, survey_abbrev))

  # mean(re_long$year)
  re_long <- re_long |>
    mutate(
      decade = (year - 2012) / 10,
      re100 = re,
      spp_by_survey = paste(species_common_name, survey_abbrev)
    )

  survs <- unique(re_long$survey_abbrev)

  re_long |>
    filter(`Restriction type` == "re_restr") |>
    ggplot(aes(year, re100, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

  re_long |>
    filter(`Restriction type` == "re_shrunk") |>
    ggplot(aes(year, re100, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

  index |>
    filter(type == "Restricted and shrunk") |>
    ggplot(aes(year, cv, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

  index |>
    filter(type == "Status quo") |>
    ggplot(aes(year, cv, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

coefs <- re_long |> group_by(survey_abbrev, species_common_name) |>
    group_split() |>
    purrr::map_dfr(function(dat) {

      d1 <- filter(dat, `Restriction type` == "re_restr")
      if (sum(is.na(d1$re)) == 0L) {
        m1 <- tryCatch(lm(re100 ~ 1 + decade, data = d1), error = function(err) NA)
      }

      d2 <- filter(dat, `Restriction type` == "re_shrunk")
      if (sum(is.na(d2$re)) == 0L) {
        m2 <- tryCatch(lm(re100 ~ 1 + decade, data = d2), error = function(err) NA)
      }

      if (sum(is.na(d1$re)) == 0L) {
        m3 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d1), error = function(err) NA)
      }
      if (sum(is.na(d2$re)) == 0L) {
        m4 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d2), error = function(err) NA)
      }

      tibble(
        survey_abbrev = dat$survey_abbrev[1],
        species_common_name = dat$species_common_name[1],
        `(Intercept)` = c(
          tryCatch(m1$coefficients[[1]], error = function(err) NA),
          tryCatch(m2$coefficients[[1]], error = function(err) NA)
        ),
        slope_re = c(
          tryCatch(m1$coefficients[["decade"]], error = function(err) NA),
          tryCatch(m2$coefficients[["decade"]], error = function(err) NA)
        ),
        se_slope_re = c(
          tryCatch(summary(m1)$coefficients[2, 2], error = function(err) NA),
          tryCatch(summary(m2)$coefficients[2, 2], error = function(err) NA)
        ),
        slope_mpa = c(
          tryCatch(m3$coefficients[["decade"]], error = function(err) NA),
          tryCatch(m4$coefficients[["decade"]], error = function(err) NA)
        ),
        restr_clean = c("Same survey domain", "Shrunk survey domain")
      )
    })

cvraw <- filter(cvraw, !is.na(restr_clean))
cvratio2 <- filter(cvratio2, !is.na(restr_clean))

# join everything together
cvdata <- left_join(cvratio2, mare2) |>
    left_join(coefs) |>
    left_join(cvraw) |>
    left_join(prop_mpa)

  # grab trend for each species:
  trends_squo <- filter(index, type == "Status quo") |>
    group_by(species_common_name, survey_abbrev) |>
    group_split() |>
    purrr::map_dfr(function(.x) {
      .x$est[.x$est == 0] <- min(.x$est, na.rm = TRUE) # b/c log NOTE this TODO
      if (sum(is.na(.x$est)) == 0L && sum(is.infinite(.x$est)) == 0L && sum(.x$est == 0) == 0L) {
        m <- lm(log(est) ~ I(year * 10), data = .x)
        data.frame(species_common_name = unique(.x$species_common_name), survey_abbrev = unique(.x$survey_abbrev), slope_squo_decade = coef(m)[[2]])
      } else {
        data.frame(species_common_name = unique(.x$species_common_name), survey_abbrev = unique(.x$survey_abbrev), slope_squo_decade = NA_real_)
      }
    })

  cvdata <- left_join(cvdata, trends_squo)

  if (survey == "HBLL") {
    cvdata$est_type <- "geostat"
    saveRDS(cvdata, "data-generated/hbll-cv-w-lm-slopes.rds")
  }
  if (survey == "SYN") {
    cvdata$est_type <- "geostat"
    saveRDS(cvdata, "data-generated/syn-cv-w-lm-slopes.rds")
  }
  if (survey == "design-boot") {
    cvdata$est_type <- "bootstrap"
    saveRDS(cvdata, "data-generated/design-boot-cv-w-lm-slopes.rds")
  }
  if (survey == "design-cochran") {
    cvdata$est_type <- "cochran"
    saveRDS(cvdata, "data-generated/design-cochran-cv-w-lm-slopes.rds")
  }

  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "mare_interp_ratio"))+ geom_hline(yintercept = 1))
  #
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "mare_interp_ratio"))+ geom_hline(yintercept = 1))
}
