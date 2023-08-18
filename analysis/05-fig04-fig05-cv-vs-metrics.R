library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long2.rds")
metrics_wide <- readRDS("data-generated/metrics-wide2.rds")

metrics_long <- filter(metrics_long, !survey_abbrev %in% c("SYN QCS, SYN HS"))
metrics_wide <- filter(metrics_wide, !survey_abbrev %in% c("SYN QCS, SYN HS"))

g <- metrics_wide |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
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
  stat_smooth(
    se = TRUE, alpha = 0.2, method = "glm",
    method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE
  ) +
  theme(
    legend.position = c(0.23, 0.8),
    strip.placement = "outside"
  ) +
  coord_cartesian(ylim = c(0, 0.5), xlim = c(0.08, 0.8)) +
  scale_x_continuous(breaks = c(0.1, 0.2, 0.4, 0.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .01)))
g
cv_fig <- g # for assembly later...
ggsave("figs/cv-status-quo-mare.pdf", width = 4.2, height = 3.6)

g <- metrics_long |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
  filter(measure != "cv") |>
  mutate(est = if_else(measure == "slope_re", abs(est), est)) |>
  mutate(est = if_else(measure == "cv_perc" & est < 0, 0.01, est)) |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, .04))) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = restricted_cols) +
  facet_wrap(
    ~measure_clean,
    nrow = 1L,
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(
    se = TRUE, alpha = 0.15, method = "glm",
    method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.25, 0.918),
    strip.placement = "outside"
  ) +
  coord_cartesian(xlim = c(0.06, max(metrics_long$prop_mpa, na.rm = TRUE) * 1.04), expand = FALSE)

g + facet_wrap(~measure_clean, scales = "free_y", ncol = 3) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.1, 0.8),
    strip.placement = "outside"
  )
ggsave("figs/prop-mpa-vs-metrics-wide.pdf", width = 10, height = 3.5)
ggsave("figs/prop-mpa-vs-metrics-wide.png", width = 10, height = 3.5)

# combine these two into one... -------------

make_panel <- function(dat, xvar, yvar, xlab = "CV of 'Status quo' index", ylab = "MARE (accuracy loss)") {
  g <- dat |>
    filter(type %in% "Restricted and shrunk") |>
    filter(est_type %in% "geostat") |>
    ggplot(
      aes({{ xvar }}, {{ yvar }}, colour = survey_abbrev)
    ) +
    geom_point(pch = 21, alpha = 1) +
    geom_point(pch = 19, alpha = 0.2) +
    xlab(xlab) +
    ylab(ylab) +
    guides(shape = "none") +
    scale_colour_manual(name = "Survey", values = restricted_cols) +
    scale_x_log10() +
    stat_smooth(se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
    theme(
      legend.position = c(0.23, 0.7),
      strip.placement = "outside"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, .01)))
  g
}
g1 <- make_panel(metrics_wide, orig_cv_mean, mare_med) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "a") +
  coord_cartesian(ylim = c(0, 0.4))
g2 <- make_panel(metrics_wide, prop_mpa, mare_med, xlab = "Proportion stock in MPA", ylab = "MARE (accuracy loss)") +
  guides(colour = "none") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "b") +
  coord_cartesian(ylim = c(0, 0.4))

g3 <- metrics_wide |>
  mutate(cv_perc_med = ifelse(cv_perc_med < 0, 0.01, cv_perc_med)) |>
  make_panel(prop_mpa, cv_perc_med, xlab = "Proportion stock in MPA", ylab = "% increase CV (precision loss)") +
  guides(colour = "none") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "c") +

  coord_cartesian(ylim = c(0, 50))
cowplot::plot_grid(g1, g2, g3, nrow = 1, align = "h", axis = "b")
ggsave("figs/metrics-cross-plot1.pdf", width = 9, height = 2.75)

# -------------

# design-based instead!?

metrics_long2 <- readRDS("data-generated/metrics-long2.rds")
metrics_long2 <- filter(metrics_long2, !survey_abbrev %in% c("SYN QCS, SYN HS"))

g <- metrics_long2 |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "bootstrap") |>
  mutate(est = if_else(measure == "slope_re", abs(est), est)) |>
  mutate(est = if_else(measure == "cv_perc" & est < 0, 0.01, est)) |>
  filter(measure != "cv") |>
  ggplot(
    aes(prop_mpa_set, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = RColorBrewer::brewer.pal(4, "Set2")) +
  facet_grid(
    rows = vars(measure_clean),
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(
    se = FALSE, alpha = 0.15, method = "glm",
    method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.25, 0.918),
    strip.placement = "outside"
  ) +
  coord_cartesian(xlim = c(0.05, max(metrics_long$prop_mpa_set, na.rm = TRUE)))

g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")
g
ggsave("figs/prop-mpa-vs-metrics-design.pdf", width = 3.70, height = 6.5)
ggsave("figs/prop-mpa-vs-metrics-design.png", width = 3.70, height = 6.5)
