library(dplyr)
library(ggplot2)
ggplot2::theme_set(ggsidekick::theme_sleek())

index <- readRDS("data-generated/index-filtered.rds")
metrics <- readRDS("data-generated/metrics-wide.rds")

m <- select(metrics, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()
index <- index |> left_join(m, by = join_by(species_common_name, survey_abbrev))

dodge_width <- 1

make_ts_plot <- function(survey_keep, ncol = NULL) {
  g <- index |>
    filter(survey_abbrev %in% survey_keep) |>
    mutate(species_common_name_lower = tolower(species_common_name)) |>
    filter(type %in% c("Status quo", "Restricted and shrunk")) |>
    group_by(type, species_common_name) |>
    mutate(lwr = lwr / exp(mean(log(est)))) |>
    mutate(upr = upr / exp(mean(log(est)))) |>
    mutate(est = est / exp(mean(log(est)))) |>
    mutate(upr = ifelse(upr > 2 * max(est), 2 * max(est), upr)) |>
    mutate(species_common_name = gsub("Rougheye/Blackspotted Rockfish", "Rougheye/Black. Rockfish", species_common_name)) |>
    arrange(-prop_mpa, species_common_name) |>
    ggplot(aes(year, est,
      ymin = lwr, ymax = upr,
      colour = type
    )) +
    geom_linerange(position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width), pch = 21, size = 1.8) +
    geom_point(position = position_dodge(width = dodge_width), pch = 20, size = 1.8, alpha = 0.2) +
    coord_cartesian(
      expand = FALSE,
      xlim = range(index$year) + c(-0.5, 0.5),
      ylim = c(0, NA)
    ) +
    labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
    scale_colour_manual(values = c("grey50", "red")) +
    facet_wrap(~ forcats::fct_inorder(species_common_name), scales = "free_y", ncol = ncol) +
    ylab("Relative abundance or biomass") +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    theme(
      strip.text = element_text(colour = "black"),
      legend.position = "top", axis.text.y = element_text(size = 8)
    )
  g
}

g <- make_ts_plot("SYN QCS, SYN HS", ncol = NULL) + ylab("Relative biomass")
ggsave("figs/ts-qcs-hs.pdf", width = 12, height = 9)
g <- make_ts_plot("SYN WCHG") + ylab("Relative biomass")
ggsave("figs/ts-wchg.pdf", width = 12, height = 8)
g <- make_ts_plot("HBLL OUT N") + ylab("Relative abundance")
ggsave("figs/ts-hbll.pdf", width = 11, height = 7.5)
