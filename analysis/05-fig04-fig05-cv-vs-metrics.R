library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long2.rds")
metrics_wide <- readRDS("data-generated/metrics-wide2.rds")

metrics_long <- filter(metrics_long, !survey_abbrev %in% c("SYN QCS", "SYN HS"))
metrics_wide <- filter(metrics_wide, !survey_abbrev %in% c("SYN QCS", "SYN HS"))

# m <- select(metrics_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa, slope_squo_decade) |>
#   distinct()
#
# metrics_long <- metrics_long |>
#   left_join(m, by = join_by(species_common_name, survey_abbrev))
#
# metrics_long <- metrics_long |>
#   mutate(survey_abbrev = factor(survey_abbrev,
#     levels = c(
#       "SYN WCHG",
#       "HBLL OUT N",
#       "SYN QCS, SYN HS"))
#   )

g <- metrics_wide |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
  # filter(grepl("mare", measure)) |>
  ggplot(
    aes(orig_cv_mean, mare_med, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("CV of 'Status quo' index") +
  ylab("MARE (accuracy loss)") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = restricted_cols) +
  scale_x_log10() +
  stat_smooth(se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
  theme(
    legend.position = c(0.25, 0.85),
    strip.placement = "outside"
  ) +
  # coord_cartesian(ylim = c(0, 0.9), xlim = c(0.08, 0.85)) +
  coord_cartesian(ylim = c(0, 0.5), xlim = c(0.08, 0.8)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.4, 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .01)))
g
ggsave("figs/cv-status-quo-mare.pdf", width = 4.2, height = 3.6)

g <- metrics_long |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
  filter(measure != "cv") |>
  # mutate(est = ifelse(grepl("trend", measure), abs(est), est)) |>
  # mutate(measure = ifelse(grepl("trend", measure), "| RE trend |\n(absolute trend bias)", measure)) |>
  # mutate(cv_perc_med = ifelse(cv_perc_med <= 0, 0.01, cv_perc_med)) |>
  # filter(est = ifelse(grepl("CV", measure) & est <= 0) |>
  # mutate(est = ifelse(grepl("CV", measure) & est <= 0, 999, est)) |>
  # mutate(measure = factor(measure, levels =
  #     c("% increase CV\n(precision loss)",
  #       "MARE\n(accuracy loss)",
  #       "| RE trend |\n(absolute trend bias)"
  #       ))) |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, .04))) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = restricted_cols) +
  # scale_colour_brewer(palette = "Set2") +
  facet_grid(
    rows = vars(measure_clean),
    switch = "y",
    scales = "free_y"
  ) +
  scale_x_log10() +
  # scale_y_log10() +
  stat_smooth(se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.25, 0.92),
    strip.placement = "outside"
  )

g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")

ggsave("figs/prop-mpa-vs-metrics.pdf", width = 3.70, height = 6.5)
ggsave("figs/prop-mpa-vs-metrics.png", width = 3.70, height = 6.5)

g +   facet_wrap(~measure_clean, scales = "free_y", ncol = 3
) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.1, 0.8),
    strip.placement = "outside"
  )
ggsave("figs/prop-mpa-vs-metrics-wide.pdf", width = 10, height = 3.5)
ggsave("figs/prop-mpa-vs-metrics-wide.png", width = 10, height = 3.5)

# print("prop CV below zero")
# x <- metrics_long |>
#   filter(restr_clean %in% "Shrunk survey domain") |>
#   filter(grepl("CV", measure) & est <= 0)
# min(x$est)

# metrics_long |>
#   filter(restr_clean %in% "Shrunk survey domain") |>
#   # mutate(est = ifelse(grepl("trend", measure), abs(est), est)) |>
#   mutate(measure = ifelse(grepl("trend", measure), "RE trend\n(absolute trend bias)", measure)) |>
#   # mutate(est = ifelse(grepl("CV", measure) & est <= 0, 0.01, est)) |>
#   # filter(est = ifelse(grepl("CV", measure) & est <= 0) |>
#   # mutate(est = ifelse(grepl("CV", measure) & est <= 0, 999, est)) |>
#   mutate(measure = factor(measure, levels =
#       c("% increase CV\n(precision loss)",
#         "MARE\n(accuracy loss)",
#         "RE trend\n(absolute trend bias)"
#       ))) |>
#   ggplot(
#     aes(slope_squo_decade, est, colour = survey_abbrev)
#   ) +
#   geom_point(pch = 21, alpha = 1) +
#   scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0.015, .04))) +
#   geom_point(pch = 19, alpha = 0.2) +
#   xlab("Slope status quo per decade") +
#   guides(shape = "none") +
#   scale_colour_manual(name = "Survey", values = restricted_cols) +
#   facet_grid(
#     rows = vars(measure),
#     switch = "y",
#     scales = "free_y"
#   )
