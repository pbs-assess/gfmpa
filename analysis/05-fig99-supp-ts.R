library(dplyr)
library(ggplot2)
source("analysis/theme.R")

# index <- readRDS("data-generated/index-filtered.rds")

# index <- readRDS("data-generated/metrics-long2.rds")

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
index <- bind_rows(hbll, syn)

index <- mutate(index, species_common_name = gsub(
  "rougheye/blackspotted rockfish complex",
  "rougheye/blackspotted rockfish", species_common_name
)) |>
  filter(!(year == 2014 & survey_abbrev == "SYN WCHG"))

metrics <- readRDS("data-generated/metrics-wide2.rds")

prop <- metrics |>
  select(
    prop_mpa, species_common_name,
    survey_abbrev
  ) |>
  distinct() |>
  mutate(species_common_name = tolower(species_common_name))

index$prop_mpa <- NULL
index <- left_join(index, prop)

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
      xlim = range(index$year) + c(-0.5, 0.5)
      # ylim = c(0, NA)
    ) +
    labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
    scale_colour_manual(values = c("grey50", "red")) +
    facet_wrap(~ forcats::fct_inorder(stringr::str_to_title(species_common_name)), scales = "free_y", ncol = ncol) +
    ylab("Relative abundance or biomass") +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    theme(
      strip.text = element_text(colour = "black"),
      legend.position = "top", axis.text.y = element_text(size = 8)
    ) +
    # scale_y_log10() +
    geom_smooth(method = "gam", se = FALSE, alpha = 0.1, formula = y ~ s(x, k = 8), method.args = list(family = Gamma(link = "log")), lwd = 0.8)
    # geom_smooth(method = "loess", se = F, alpha = 0.5, formula = y ~ x, lwd = 0.8 )
  g
}

g <- make_ts_plot("SYN QCS", ncol = NULL) + ylab("Relative biomass")
ggsave("figs/ts-qcs.pdf", width = 12, height = 9)

g <- make_ts_plot("SYN HS", ncol = NULL) + ylab("Relative biomass")
ggsave("figs/ts-hs.pdf", width = 12, height = 9)

g <- make_ts_plot("SYN QCS, SYN HS", ncol = NULL) + ylab("Relative biomass")
ggsave("figs/ts-qcs-hs.pdf", width = 12, height = 9)

g <- make_ts_plot("SYN WCHG") + ylab("Relative biomass")
ggsave("figs/ts-wchg.pdf", width = 12, height = 8)

g <- make_ts_plot("HBLL OUT N") + ylab("Relative abundance")
ggsave("figs/ts-hbll.pdf", width = 11, height = 7.5)
