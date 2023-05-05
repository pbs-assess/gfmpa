library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long2.rds")
metrics_wide <- readRDS("data-generated/metrics-wide2.rds")

# metrics_long <- filter(metrics_long, !survey_abbrev %in% c("SYN QCS", "SYN HS"))
# metrics_wide <- filter(metrics_wide, !survey_abbrev %in% c("SYN QCS", "SYN HS"))

metrics_long <- filter(metrics_long, !survey_abbrev %in% c("SYN QCS, SYN HS"))
metrics_wide <- filter(metrics_wide, !survey_abbrev %in% c("SYN QCS, SYN HS"))

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

# restricted_cols <- RColorBrewer::brewer.pal(4, "Set2")
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
    legend.position = c(0.23, 0.8),
    strip.placement = "outside"
  ) +
  # coord_cartesian(ylim = c(0, 0.9), xlim = c(0.08, 0.85)) +
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
  # scale_colour_brewer(palette = "Set2") +
  facet_wrap(
    ~measure_clean,
    nrow = 1L,
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(se = TRUE, alpha = 0.15, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.25, 0.918),
    strip.placement = "outside"
  ) +
  coord_cartesian(xlim = c(0.06, max(metrics_long$prop_mpa, na.rm = TRUE)* 1.04), expand = FALSE)

# g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = c("b", "c"))
# g

# ggsave("figs/prop-mpa-vs-metrics.pdf", width = 3.70, height = 6.5)
# ggsave("figs/prop-mpa-vs-metrics.png", width = 3.70, height = 6.5)

# cv_fig + tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "a")

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

# combine these two into one... -------------

make_panel <- function(dat, xvar, yvar, xlab = "CV of 'Status quo' index", ylab = "MARE (accuracy loss)", .est_type = "geostat") {
  g <- dat |>
    filter(type %in% "Restricted and shrunk") |>
    filter(est_type %in% .est_type) |>
    # filter(grepl("mare", measure)) |>
    ggplot(
      aes({{xvar}}, {{yvar}}, colour = survey_abbrev)
    ) +
    geom_point(pch = 21, alpha = 1) +
    geom_point(pch = 19, alpha = 0.2) +
    xlab(xlab) +
    ylab(ylab) +
    guides(shape = "none") +
    scale_colour_manual(name = "Survey", values = restricted_cols) +
    scale_x_log10() +
    geom_smooth(se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
    geom_smooth(mapping = aes(colour = NULL, group = NULL), se = TRUE, alpha = 0.2, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE, colour = "grey20") +
    theme(
      legend.position = c(0.23, 0.67),
      strip.placement = "outside"
    ) +
    # coord_cartesian(ylim = c(0, 0.9), xlim = c(0.08, 0.85)) +
    # coord_cartesian(ylim = c(0, 0.5), xlim = c(0.08, 0.8)) +
    # scale_x_continuous(breaks = c(0.1, 0.2, 0.4, 0.8)) +
    scale_y_continuous(expand = expansion(mult = c(0, .01)))
  g
}
g1 <- make_panel(metrics_wide, orig_cv_mean, mare_med) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "a") +
  coord_cartesian(ylim = c(0, .4), expand = FALSE)
g2 <- make_panel(metrics_wide, prop_mpa, mare_med, xlab = "Proportion stock in MPA", ylab = "MARE (accuracy loss)") +
  guides(colour = "none") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "b") +
  coord_cartesian(ylim = c(0, .4), expand = FALSE)
g3 <- metrics_wide |>
  mutate(cv_perc_med = ifelse(cv_perc_med < 0, 0.01, cv_perc_med)) |>
  make_panel(prop_mpa, cv_perc_med, xlab = "Proportion stock in MPA", ylab = "% increase CV (precision loss)") +
  guides(colour = "none") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "c") +
  coord_cartesian(ylim = c(0, 50), expand = FALSE)
cowplot::plot_grid(g1, g2, g3, nrow = 1, align = "h", axis = "b")

ggsave("figs/metrics-cross-plot1.pdf", width = 8.5, height = 2.55)

# design-based instead!?

metrics_long2 <- readRDS("data-generated/metrics-long2.rds")
metrics_long2 <- filter(metrics_long2, !survey_abbrev %in% c("SYN QCS, SYN HS"))

tokeep <- metrics_long2 |>
  filter(est_type %in% c("bootstrap", "geostat")) |>
  filter(type %in% "Restricted and shrunk") |>
  select(est_type, species_common_name, survey_abbrev) |>
  distinct() |>
  group_by(species_common_name, survey_abbrev) |>
  summarise(n_ests = n()) |>
  filter(n_ests == 2L) |>
  select(n_ests)

met_dat <- metrics_long2 |>
  semi_join(tokeep) |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% c("bootstrap", "geostat")) |>
  mutate(est_type = gsub("bootstrap", "Design-based", est_type)) |>
  mutate(est_type = gsub("geostat", "Geostatistical", est_type)) |>
  mutate(est = if_else(measure == "slope_re", abs(est), est)) |>
  mutate(est = if_else(measure == "cv_perc" & est < 0, 0.01, est)) |>
  filter(measure != "cv")

g <- met_dat |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  # scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, .04))) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  # scale_colour_manual(name = "Survey", values = restricted_cols) +
  scale_colour_manual(name = "Survey", values = RColorBrewer::brewer.pal(4, "Set2")) +
  # scale_colour_brewer(palette = "Set2") +
  facet_grid(
    rows = vars(measure_clean),
    cols = vars(est_type),
    # switch = "y",
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(se = FALSE, alpha = 0.15, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE) +
  stat_smooth(mapping = aes(colour = NULL), se = FALSE, alpha = 0.15, method = "glm", method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE, colour = "grey20") +
  # stat_smooth(se = TRUE, alpha = 0.2, method = "glm", formula = y ~ x, show.legend = FALSE) +
  # geom_smooth(se = FALSE) +
  theme(
    # axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    legend.position = c(0.15, 0.918),
    strip.placement = "outside"
  ) +
  ggrepel::geom_text_repel(aes(label = species_common_name), size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = c(0.015, max(metrics_long$prop_mpa, na.rm = TRUE))) +
  ylab("Metric value")

g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")
# g
ggsave("figs/prop-mpa-vs-metrics-design2.pdf", width = 7.5, height = 9.5)
ggsave("figs/prop-mpa-vs-metrics-design2.png", width = 7.5, height = 9.5)

# stats? ----------------------

glimpse(met_dat)

d <- filter(met_dat, measure == "mare")
plot(log(d$prop_mpa), log(d$est))
plot(log(d$prop_mpa), log(d$est))

fit <- glm(est ~ log(prop_mpa_set + 0.001) * est_type, data = d, family = Gamma(link = "log"))
arm::display(fit)

d$est_type <- factor(d$est_type, levels = c("Geostatistical", "Design-based"))
fit <- glm(est ~ log(prop_mpa) * est_type, data = d, family = Gamma(link = "log"))
arm::display(fit)

fit <- glm(est ~ log(prop_mpa_set + 0.001) * est_type, data = filter(met_dat, measure == "slope_re"), family = Gamma(link = "log"))
arm::display(fit)

d <- filter(met_dat, measure == "cv_perc")
plot(log(d$prop_mpa), log(d$est))
plot(log(d$prop_mpa), log(d$est))

fit <- glm(est ~ log(prop_mpa + 0.001) * est_type, data = filter(met_dat, measure == "cv_perc"), family = Gamma(link = "log"))
arm::display(fit)

get_slope <- function(.measure = "mare",
  levels = c("Geostatistical", "Design-based")) {
  d <- filter(met_dat, measure == .measure)
  d$est_type <- factor(d$est_type, levels = levels)
  fit <- glm(est ~ log(prop_mpa) * est_type, data = d, family = Gamma(link = "log"))
  b <- coef(fit)[[2]]
  suppressMessages(ci <- confint(fit))
  data.frame(b = b, lwr = ci[2,1], upr = ci[2,2], measure = .measure, base_level = levels[1], stringsAsFactors = FALSE)
}

slopes <- list()
slopes[[1]] <- get_slope("mare")
slopes[[2]] <- get_slope("slope_re")
slopes[[3]] <- get_slope("cv_perc")

lv <- c("Design-based", "Geostatistical")
slopes[[4]] <- get_slope("mare", lv)
slopes[[5]] <- get_slope("slope_re", lv)
slopes[[6]] <- get_slope("cv_perc", lv)

slopes <- bind_rows(slopes)
slopes[[1]] <- sdmTMB:::mround(slopes[[1]], 2)
slopes[[2]] <- sdmTMB:::mround(slopes[[2]], 2)
slopes[[3]] <- sdmTMB:::mround(slopes[[3]], 2)
slopes

slopes <- mutate(slopes, text = paste0(b, "\\% (95\\% CI: ", lwr, "--", upr, ")"))
slopes

saveRDS(slopes, "data-generated/metrics-slopes-table.rds")
