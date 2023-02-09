# get slopes of relative error
library(dplyr)
library(ggplot2)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

for (survey in c("HBLL", "SYN")) {
  # for now using delta-gamma?
  if (survey == "HBLL") y <- readRDS(file = "data-generated/index-hbll-geo-clean-binomial_gamma.rds")
  if (survey == "SYN") y <- readRDS(file = "data-generated/index-syn-geo-clean-binomial_gamma.rds")

  index$species_common_name <- stringr::str_to_title(index$species_common_name)

  # index <- filter(index, !species_common_name == "Shortraker Rockfish")
  # remove Shortspine as only HBLL spp below 10% occurance that converges
  # index <- filter(index, !(survey_abbrev == "HBLL OUT N" & species_common_name == "Shortspine Thornyhead"))
  index <- mutate(index, species_common_name = gsub(
    "Rougheye/Blackspotted Rockfish Complex",
    "Rougheye/Blackspotted Rockfish", species_common_name
  ))

  # get cv ratios and prop mpa
  cv2 <- index %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      prop_mpa = est[type == "MPA only"] / est[type == "Status quo"],
      cv_ratio_restr = cv[type == "Restricted"] /
        cv[type == "Status quo"],
      cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
        cv[type == "Status quo"]
    )

  cv_long <- cv2 %>%
    select(-prop_mpa) %>%
    tidyr::pivot_longer(starts_with("cv"),
      names_to = "Restriction type", values_to = "CV ratio"
    )

  prop_mpa <- cv2 %>%
    select(species_common_name, survey_abbrev, year, prop_mpa) %>%
    group_by(species_common_name, survey_abbrev) %>%
    summarise(
      prop_mpa = mean(prop_mpa)
    )

  lu <- tibble(
    "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  cvratio <- cv_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
    left_join(lu)



  # get raw cvs
  lu <- tibble(
    "type" = c("Restricted", "Restricted and shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  cvraw <- index %>%
    filter(type != "MPA only") %>%
    group_by(species_common_name, survey_abbrev, type) %>%
    summarise(
      cv_orig = mean(orig_cv),
      cv_max = max(cv),
      cv_mean = mean(cv)
    ) %>%
    filter(type != "Status quo") %>%
    left_join(., lu)


  # bad year in SYN, removed from plot so removing here too
  if (survey == "SYN") {
    index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
  }

  re <- index %>%
    group_by(species_common_name, survey_abbrev, type) %>%
    mutate(est = est / exp(mean(log(est)))) %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      prop_mpa = est[type == "MPA only"] / est[type == "Status quo"],
      re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
      re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
    )

  lu <- tibble(
    "Restriction type" = c("re_restr", "re_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  re_long <- re %>%
    tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") %>%
    left_join(lu)

  mare <- re_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    left_join(lu)


  # rename variables for other previously calculated indices
  cvratio2 <- cvratio %>%
    rename(cv_ratio = est, cv_lwr = lwr, cv_upr = upr) %>%
    select(-`Restriction type`)
  mare2 <- mare %>%
    rename(mare = est, mare_lwr = lwr, mare_upr = upr) %>%
    mutate(spp_by_survey = paste(species_common_name, survey_abbrev))



  # mean(re_long$year)
  re_long <- re_long %>%
    mutate(
      decade = (year - 2012) / 10,
      re100 = re,
      spp_by_survey = paste(species_common_name, survey_abbrev)
    )

  survs <- unique(re_long$survey_abbrev)

  re_long %>%
    filter(`Restriction type` == "re_restr") %>%
    ggplot(aes(year, re100, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

  index %>%
    filter(type == "Restricted and shrunk") %>%
    ggplot(aes(year, cv, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")

  index %>%
    filter(type == "Status quo") %>%
    ggplot(aes(year, cv, colour = species_common_name)) +
    geom_line() +
    facet_wrap(~survey_abbrev, dir = "v")


  # # calculate slope from fig 2 relative error plots using random effects
  # library(lme4)
  #
  # by_surv <- list()
  # for (j in seq_along(survs)){
  #
  # restr <- re_long %>% filter(survey_abbrev == survs[j]) %>%
  #   filter(`Restriction type` == "re_restr")
  # m_restr <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = restr)
  # m_restr
  #
  # shrunk <- re_long %>% filter(survey_abbrev == survs[j]) %>%
  #   filter(`Restriction type` == "re_shrunk")
  # m_shrunk  <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = shrunk)
  # m_shrunk
  #
  # # These values are a combination of the fixed effects and the variance components
  # restr_coefs <- coef(m_restr)$species_common_name %>%
  #   tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
  # restr_coefs$restr_clean <- "Same survey domain"
  # restr_coefs$survey_abbrev <- survs[j]
  #
  # shrunk_coefs <- coef(m_shrunk)$species_common_name %>%
  #   tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
  # shrunk_coefs$restr_clean <- "Shrunk survey domain"
  # shrunk_coefs$survey_abbrev <- survs[j]
  #
  # by_surv[[j]] <- bind_rows(restr_coefs, shrunk_coefs)
  # }
  # coefs <- do.call(rbind, by_surv)
  #
  # # join everything together
  # cvdata <- left_join(cvratio2, mare2) %>% left_join(., coefs) %>%
  #   left_join(., cvraw) %>% left_join(., prop_mpa)
  #
  # if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-re-slopes.rds")
  # if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-re-slopes.rds")



  # OR use simple linear models

  survs <- unique(re_long$survey_abbrev)
  spp <- unique(re_long$species_common_name)
  sps <- list()

  # get_slopes <- function(dat) {
  #   d1 <- filter(dat, `Restriction type` == "re_restr")
  #   m1 <- tryCatch(lm(re100 ~ 1 + decade, data = d1), error = function(err) NA)
  #
  #   d2 <- filter(dat, `Restriction type` == "re_shrunk")
  #   m2 <- tryCatch(lm(re100 ~ 1 + decade, data = d2), error = function(err) NA)
  #
  #   m3 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d1), error = function(err) NA)
  #   m4 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d2), error = function(err) NA)
  #
  #   tibble(
  #     survey_abbrev = c(survs[j], survs[j]),
  #     species_common_name = c(spp[i], spp[i]),
  #     `(Intercept)` = c(
  #       tryCatch(m1$coefficients[[1]], error = function(err) NA),
  #       tryCatch(m2$coefficients[[1]], error = function(err) NA)
  #     ),
  #     slope_re = c(
  #       tryCatch(m1$coefficients[["decade"]], error = function(err) NA),
  #       tryCatch(m2$coefficients[["decade"]], error = function(err) NA)
  #     ),
  #     se_slope_re = c(
  #       tryCatch(summary(m1)$coefficients[2,2], error = function(err) NA),
  #       tryCatch(summary(m2)$coefficients[2,2], error = function(err) NA)
  #     ),
  #     slope_mpa = c(
  #       tryCatch(m3$coefficients[["decade"]], error = function(err) NA),
  #       tryCatch(m4$coefficients[["decade"]], error = function(err) NA)
  #     ),
  #     restr_clean = c("Same survey domain", "Shrunk survey domain")
  #   )
  # }
  # coefs <- group_by(re_long, survey_abbrev, species_common_name) %>%
  #   group_split() %>%
  #   purrr::map_dfr(get_slopes)


  for (i in seq_along(spp)) {
    for (j in seq_along(survs)) {
      dat <- re_long %>%
        filter(survey_abbrev == survs[j]) %>%
        filter(species_common_name == spp[i])

      d1 <- filter(dat, `Restriction type` == "re_restr")
      m1 <- tryCatch(lm(re100 ~ 1 + decade, data = d1), error = function(err) NA)

      d2 <- filter(dat, `Restriction type` == "re_shrunk")
      m2 <- tryCatch(lm(re100 ~ 1 + decade, data = d2), error = function(err) NA)

      m3 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d1), error = function(err) NA)
      m4 <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d2), error = function(err) NA)

      k <- length(spp) * (j - 1) + i

      sps[[k]] <- tibble(
        survey_abbrev = c(survs[j], survs[j]),
        species_common_name = c(spp[i], spp[i]),
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
    }
  }
  coefs <- do.call(rbind, sps)

  # join everything together
  cvdata <- left_join(cvratio2, mare2) %>%
    left_join(., coefs) %>%
    left_join(., cvraw) %>%
    left_join(., prop_mpa)

  if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-lm-slopes.rds")
  if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-lm-slopes.rds")


  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "mare_interp_ratio"))+ geom_hline(yintercept = 1))
  #
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
  # (g <- ggplot(d) + geom_point(aes_string("cv_orig", "mare_interp_ratio"))+ geom_hline(yintercept = 1))
}
