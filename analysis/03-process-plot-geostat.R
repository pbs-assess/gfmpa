library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)
dir.create("figs", showWarnings = FALSE)

survey <- "HBLL"
fam <- "binomial_gamma"

for (survey in c("HBLL", "SYN")) {
  cat("survey: ", survey, "\n")
  cat("family: ", fam, "\n")

  if (survey == "HBLL") index <- readRDS("data-generated/index-hbll-geo-clean.rds")
  if (survey == "SYN") index <- readRDS("data-generated/index-syn-geo-clean.rds")

  index$species_common_name <- stringr::str_to_title(index$species_common_name)

  g <- ggplot(index, aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
    geom_line(linewidth = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Type", fill = "Type") +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2")

  if (survey == "HBLL") {
    g <- g +
      facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
      ylab("Relative numbers of fish")
  }
  if (survey == "SYN") {
    g <- g +
      facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") +
      ylab("Relative biomass")
  }

  if (survey == "HBLL") ggsave(paste0("figs/index-hbll-geo-restricted.pdf"), width = 12, height = 8)
  if (survey == "SYN") ggsave(paste0("figs/index-syn-geo-restricted.pdf"), width = 9, height = 60, limitsize = FALSE)

  syn_highlights <- c(
    "Arrowtooth Flounder",
    "Dover Sole",
    "Flathead Sole",
    "Longnose Skate",
    "North Pacific Spiny Dogfish",
    "Pacific Halibut",
    "Petrale Sole",
    "Redbanded Rockfish",
    "Rex Sole",
    "Blackspotted Rockfish",
    "Sandpaper Skate",
    "Shortspine Thornyhead",
    "Yellowtail Rockfish",
    "Widow Rockfish",
    "Walleye Pollock",
    "Pacific Ocean Perch",
    "Pacific Cod",
    "Redstripe Rockfish"
  )

  # focus on some
  if (survey == "SYN") {
    g <- index %>%
      filter(species_common_name %in% syn_highlights) %>%
      ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
      geom_line(lwd = 0.9) +
      geom_ribbon(alpha = 0.2, colour = NA) +
      labs(x = "Year", colour = "Type", fill = "Type") +
      scale_color_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2") +
      facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") +
      ylab("Relative biomass") +
      theme(strip.text.y = element_text(size = 7))
    ggsave("figs/index-syn-geo-restricted-highlight2.pdf", width = 8, height = 15)
  }

  x <- index %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      cv_ratio_restr = cv[type == "Restricted"] /
        cv[type == "Status quo"],
      cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
        cv[type == "Status quo"]
    )

  x_long <- x %>%
    tidyr::pivot_longer(starts_with("cv"),
      names_to = "Restriction type", values_to = "CV ratio"
    )

  # x_long %>%
  #   ggplot(aes(`CV ratio`)) +
  #   facet_wrap(~survey_abbrev) +
  #   geom_histogram() +
  #   geom_vline(xintercept = 1, col = "red") +
  #   facet_wrap(vars(`Restriction type`)) +
  #   coord_cartesian(expand = FALSE)
  #
  # if (survey == "HBLL") ggsave("figs/index-geo-hbll-cv-ratios.pdf", width = 8, height = 4)
  # if (survey == "SYN") ggsave("figs/index-geo-syn-cv-ratios.pdf", width = 8, height = 4)

  x_long %>%
    group_by(`Restriction type`, survey_abbrev) %>%
    summarise(mean_ratio = mean(`CV ratio`)) %>%
    knitr::kable(digits = 2)

  lu <- tibble(
    "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )
  g <- x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
    left_join(lu) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean, ymin = lwr, ymax = upr)) +
    geom_hline(yintercept = 1, lty = 2, col = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    coord_flip() +
    xlab("") +
    ylab("Ratio of index CV\n(restricted/original)") +
    labs(colour = "Restricted survey domain") +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "top") +
    theme(panel.grid.major.y = element_line(colour = "grey90"))

  if (survey == "SYN") {
    g <- g + scale_y_log10() + facet_wrap(~survey_abbrev)
    g
  }
  if (survey == "HBLL") ggsave("figs/index-geo-hbll-cv-ratio-dotplot.pdf", width = 7, height = 7)
  if (survey == "SYN") ggsave("figs/index-geo-syn-cv-ratio-dotplot.pdf", width = 10, height = 8)

  # # bad year in SYN
  # if (survey == "SYN") {
  #   index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
  # }

  x <- index %>%
    group_by(species_common_name, survey_abbrev, type) %>%
    mutate(est = est / exp(mean(log(est)))) %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
      re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
    )

  lu <- tibble(
    "Restriction type" = c("re_restr", "re_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  x_long <- x %>%
    tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") %>%
    left_join(lu)

  g <- ggplot(x_long, aes(year, re, colour = restr_clean)) +
    geom_line() +
    geom_hline(yintercept = 0, lty = 2) +
    scale_color_brewer(palette = "Set2") +
    ylab("Relative error") +
    xlab("Year") +
    labs(colour = "Survey domain treatment")

  if (survey == "HBLL") {
    g <- g + facet_wrap(~species_common_name, scales = "free_y", ncol = 5)
  }
  if (survey == "SYN") {
    g <- g + facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") + theme(strip.text.y = element_text(size = 7))
  }

  if (survey == "HBLL") ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 12, height = 8)
  if (survey == "SYN") ggsave("figs/index-syn-geo-restricted-re.pdf", width = 9, height = 60, limitsize = FALSE)

  # focus on some
  if (survey == "SYN") {
    g <- x_long %>%
      filter(species_common_name %in% syn_highlights) %>%
      ggplot(aes(year, re, colour = restr_clean)) +
      geom_line() +
      geom_hline(yintercept = 0, lty = 2) +
      facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") +
      scale_color_brewer(palette = "Set2") +
      ylab("Relative error") +
      xlab("Year") +
      labs(colour = "Survey domain treatment") +
      theme(strip.text.y = element_text(size = 7))
    ggsave("figs/index-syn-geo-restricted-re-highlights.pdf", width = 8, height = 15)
  }

  g <- x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    left_join(lu) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean)) +
    geom_point(position = position_dodge(width = 0)) +
    coord_flip() +
    xlab("") +
    ylab("Median absolute relative error (MARE)\n(restricted compared to original)") +
    labs(colour = "Restricted survey domain") +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "top") +
    theme(panel.grid.major.y = element_line(colour = "grey90")) +
    guides(colour = guide_legend(nrow = 2L))

  if (survey == "HBLL") ggsave("figs/index-geo-hbll-mare-dotplot.pdf", width = 6, height = 7)

  if (survey == "SYN") {
    g <- g + facet_wrap(~survey_abbrev)
    ggsave("figs/index-geo-syn-mare-dotplot.pdf", width = 8, height = 8)
  }

  x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    group_by(`Restriction type`, survey_abbrev) %>%
    summarise(mean = mean(est)) %>%
    knitr::kable(digits = 2)

  x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    group_by(`Restriction type`, survey_abbrev, species_common_name) %>%
    summarise(mean = mean(est)) %>%
    arrange(mean) %>%
    knitr::kable(digits = 2)
}
