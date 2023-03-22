library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long.rds")
metrics_wide <- readRDS("data-generated/metrics-wide.rds")

m <- select(metrics_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()

metrics_long <- metrics_long |>
  left_join(m, by = join_by(species_common_name, survey_abbrev))

metrics_long <- metrics_long |>
  mutate(survey_abbrev = factor(survey_abbrev,
    levels = c(
      "SYN WCHG",
      "HBLL OUT N",
      "SYN QCS, SYN HS"))
  )

g <- metrics_long |>
  filter(restr_clean %in% "Shrunk survey domain") |>
  filter(grepl("MARE", measure)) |>
  ggplot(
    aes(cv_orig, est, colour = survey_abbrev)
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
  coord_cartesian(ylim = c(0, NA)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.4, 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .01)))
ggsave("figs/cv-status-quo-mare.pdf", width = 4.2, height = 3.6)

g <- metrics_long |>
  filter(restr_clean %in% "Shrunk survey domain") |>
  mutate(est = ifelse(grepl("trend", measure), abs(est), est)) |>
  mutate(measure = ifelse(grepl("trend", measure), "| RE trend |\n(absolute trend bias)", measure)) |>
  mutate(est = ifelse(grepl("CV", measure) & est <= 0, 0.01, est)) |>
  # filter(est = ifelse(grepl("CV", measure) & est <= 0) |>
  # mutate(est = ifelse(grepl("CV", measure) & est <= 0, 999, est)) |>
  mutate(measure = factor(measure, levels =
      c("% increase CV\n(precision loss)",
        "MARE\n(accuracy loss)",
        "| RE trend |\n(absolute trend bias)"
        ))) |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0.015, .04))) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = restricted_cols) +
  facet_grid(
    rows = vars(measure),
    switch = "y",
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.25, 0.92),
    strip.placement = "outside"
  )

g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")

ggsave("figs/prop-mpa-vs-metrics.pdf", width = 3.70, height = 6)

print("prop CV below zero")
x <- metrics_long |>
  filter(restr_clean %in% "Shrunk survey domain") |>
  filter(grepl("CV", measure) & est <= 0)
min(x$est)
