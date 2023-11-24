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

# plotting --------------------

g <- metrics_long |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
  filter(measure != "cv") |>
  mutate(est = if_else(measure == "slope_re", abs(est), est)) |>
  # mutate(est = if_else(measure == "cv_perc" & est < 0, 0.01, est)) |>
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
  # stat_smooth(
  #   se = TRUE, alpha = 0.15, method = "glm",
  #   method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE
  # ) +
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

make_panel <- function(
    dat, xvar, yvar, ribbon_select,
    xlab = "CV of 'Status quo' index", ylab = "MARE (accuracy loss)",
    .est_type = "geostat") {

  g <- dat |>
    filter(type %in% "Restricted and shrunk") |>
    filter(est_type %in% .est_type) |>
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
    geom_smooth(
      se = TRUE, alpha = 0.2, method = "glm",
      method.args = list(family = Gamma(link = "log")),
      formula = y ~ x, show.legend = FALSE
    ) +
    geom_smooth(
      mapping = aes(colour = NULL, group = NULL),
      se = TRUE, alpha = 0.2, method = "glm",
      method.args = list(family = Gamma(link = "log")),
      formula = y ~ x, show.legend = FALSE, colour = "grey20"
    ) +
    theme(
      legend.position = c(0.23, 0.67),
      strip.placement = "outside"
    )

    # scale_y_continuous(expand = expansion(mult = c(0, .01)))
  g
}
g1 <- make_panel(metrics_wide, orig_cv_mean, mare_med, "mare") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "a") +
  coord_cartesian(ylim = c(0, .4), expand = FALSE)
g2 <- make_panel(metrics_wide, prop_mpa, mare_med, "slope_re",
  xlab = "Proportion stock in MPA", ylab = "MARE (accuracy loss)"
) +
  guides(colour = "none") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "b") +
  coord_cartesian(ylim = c(0, .4), expand = FALSE)

# third one is with normal, model on variance!
met_dat <- metrics_long |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% "geostat") |>
  filter(measure != "cv") |>
  mutate(est = if_else(measure == "slope_re", abs(est), est))
met_dat <- filter(met_dat, measure == "cv_perc")

m <- glmmTMB::glmmTMB(est ~ prop_mpa, dispformula = ~ prop_mpa, data = met_dat)
nd <- data.frame(prop_mpa = seq(min(met_dat$prop_mpa), max(met_dat$prop_mpa), length.out = 300))
p <- predict(m, newdata = nd, se.fit = TRUE)
rd <- data.frame(
  prop_mpa = nd$prop_mpa,
  measure = met_dat$measure[1], mid = (p$fit),
  lwr = (p$fit - p$se.fit * 1.96),
  upr = (p$fit + p$se.fit * 1.96), stringsAsFactors = FALSE
)
rd2 <- group_by(met_dat, survey_abbrev) |>
  group_split() |>
  purrr::map_dfr(function(.x) {
    m <- glmmTMB::glmmTMB(est ~ prop_mpa, dispformula = ~ prop_mpa, data = .x)
    nd <- data.frame(prop_mpa = seq(min(.x$prop_mpa), max(.x$prop_mpa), length.out = 300))
    p <- predict(m, newdata = nd, se.fit = TRUE)
    data.frame(
      prop_mpa = nd$prop_mpa,
      survey_abbrev = .x$survey_abbrev[1],
      measure = .x$measure[1],
      mid = (p$fit),
      lwr = (p$fit - p$se.fit * 1.96),
      upr = (p$fit + p$se.fit * 1.96), stringsAsFactors = FALSE
    )
  })
row.names(rd2) <- NULL

g3 <- met_dat |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  ylab("% increase CV (precision loss") +
  guides(shape = "none", colour = "none") +
  scale_colour_manual(name = "Survey", values = restricted_cols) +
  scale_x_log10() +
  geom_ribbon(data = rd,
    mapping = aes(x = prop_mpa, y = mid, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = rd,
    mapping = aes(x = prop_mpa, y = mid), inherit.aes = FALSE, alpha = 1, lwd = 1.1) +
  geom_ribbon(data = rd2,
    mapping = aes(x = prop_mpa, y = mid, ymin = lwr, ymax = upr, group = survey_abbrev), inherit.aes = FALSE, alpha = 0.2) +
  geom_line(data = rd2,
    mapping = aes(x = prop_mpa, y = mid, colour = survey_abbrev), inherit.aes = FALSE, alpha = 1, lwd = 1.1) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "c")

g3


# g3 <- metrics_wide |>
#   # mutate(cv_perc_med = ifelse(cv_perc_med < 0, 0.01, cv_perc_med)) |>
#   make_panel(prop_mpa, cv_perc_med, "cv_perc",
#     xlab = "Proportion stock in MPA", ylab = "% increase CV (precision loss)"
#   ) +
#   guides(colour = "none") +
#   tagger::tag_facets(tag_prefix = "(", position = "tl", tag_pool = "c") +
#   # coord_cartesian(ylim = c(0, 50), expand = FALSE)
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
  # mutate(est = if_else(measure == "cv_perc" & est < 0, 0.01, est)) |>
  filter(measure != "cv")

g <- met_dat |>
  ggplot(
    aes(prop_mpa, est, colour = survey_abbrev)
  ) +
  geom_point(pch = 21, alpha = 1) +
  geom_point(pch = 19, alpha = 0.2) +
  xlab("Proportion stock in MPA") +
  guides(shape = "none") +
  scale_colour_manual(name = "Survey", values = RColorBrewer::brewer.pal(4, "Set2")) +
  facet_grid(
    rows = vars(measure_clean),
    cols = vars(est_type),
    scales = "free_y"
  ) +
  scale_x_log10() +
  stat_smooth(
    se = FALSE, alpha = 0.15, method = "glm",
    method.args = list(family = Gamma(link = "log")), formula = y ~ x, show.legend = FALSE
  ) +
  stat_smooth(
    mapping = aes(colour = NULL), se = FALSE, alpha = 0.15, method = "glm",
    method.args = list(family = Gamma(link = "log")),
    formula = y ~ x, show.legend = FALSE, colour = "grey20"
  ) +
  theme(
    axis.title.x = element_text(size = 10),
    legend.position = c(0.15, 0.918),
    strip.placement = "outside"
  ) +
  ggrepel::geom_text_repel(aes(label = species_common_name), size = 2.5, alpha = 0.6) +
  coord_cartesian(xlim = c(0.015, max(metrics_long$prop_mpa, na.rm = TRUE))) +
  ylab("Metric value")

## add on Gaussian fits for first row!
temp <- filter(met_dat, measure == "cv_perc")
line_dat <- group_by(temp, est_type) |>
  group_split() |>
  purrr::map_dfr(function(.x) {
    m <- glmmTMB::glmmTMB(est ~ prop_mpa, dispformula = ~ prop_mpa, data = .x)
    nd <- data.frame(prop_mpa = seq(min(.x$prop_mpa, na.rm = TRUE), max(.x$prop_mpa, na.rm = TRUE), length.out = 300))
    p <- predict(m, newdata = nd, se.fit = TRUE)
    data.frame(
      prop_mpa = nd$prop_mpa,
      est_type = .x$est_type[1],
      measure = .x$measure[1],
      measure_clean = .x$measure_clean[1],
      line_fit = (p$fit),
      survey_abbrev = NA
    )
  })
line_dat_survey <- group_by(temp, est_type, survey_abbrev) |>
  group_split() |>
  purrr::map_dfr(function(.x) {
    m <- glmmTMB::glmmTMB(est ~ prop_mpa, dispformula = ~ prop_mpa, data = .x)
    nd <- data.frame(prop_mpa = seq(min(.x$prop_mpa, na.rm = TRUE), max(.x$prop_mpa, na.rm = TRUE), length.out = 300))
    p <- predict(m, newdata = nd, se.fit = TRUE)
    data.frame(
      prop_mpa = nd$prop_mpa,
      est_type = .x$est_type[1],
      measure = .x$measure[1],
      measure_clean = .x$measure_clean[1],
      line_fit = (p$fit),
      survey_abbrev = .x$survey_abbrev[1]
    )
  })
g <- g + geom_line(data = line_dat, mapping = aes(y = line_fit), lwd = 1, colour = "grey20") +
  geom_line(data = line_dat_survey, mapping = aes(y = line_fit), lwd = 1)

g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")
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

fit <- glm(est ~ log(prop_mpa_set + 0.001) * est_type,
  data = filter(met_dat, measure == "slope_re"), family = Gamma(link = "log")
)
arm::display(fit)

d <- filter(met_dat, measure == "cv_perc")
plot(log(d$prop_mpa), log(d$est))
plot(log(d$prop_mpa), log(d$est))

fit <- glm(est ~ log(prop_mpa + 0.001) * est_type,
  data = filter(met_dat, measure == "cv_perc"), family = Gamma(link = "log")
)
arm::display(fit)

get_slope <- function(
    .measure = "mare",
    levels = c("Geostatistical", "Design-based")) {
  d <- filter(met_dat, measure == .measure)
  d$est_type <- factor(d$est_type, levels = levels)
  fit <- glm(est ~ log(prop_mpa) * est_type, data = d, family = Gamma(link = "log"))
  b <- coef(fit)[[2]]
  suppressMessages(ci <- confint(fit))
  data.frame(
    b = b, lwr = ci[2, 1], upr = ci[2, 2], measure = .measure,
    base_level = levels[1], stringsAsFactors = FALSE
  )
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

# CV panel fit ------------------------------

metrics_wide <- readRDS("data-generated/metrics-wide2.rds")
metrics_wide <- filter(metrics_wide, !survey_abbrev %in% c("SYN QCS, SYN HS"))
cv_dat <- metrics_wide |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type %in% c("geostat", "bootstrap"))

get_slope_cv <- function(
    levels = c("geostat", "bootstrap")) {
  d <- cv_dat
  d$est_type <- factor(d$est_type, levels = levels)
  fit <- glm(mare_med ~ log(orig_cv_mean) * est_type, data = d, family = Gamma(link = "log"))
  b <- coef(fit)[[2]]
  suppressMessages(ci <- confint(fit))
  data.frame(
    b = b, lwr = ci[2, 1], upr = ci[2, 2],
    measure = "mare", base_level = levels[1], stringsAsFactors = FALSE
  )
}

slopes_cv <- get_slope_cv()
slopes_cv[[1]] <- sdmTMB:::mround(slopes_cv[[1]], 2)
slopes_cv[[2]] <- sdmTMB:::mround(slopes_cv[[2]], 2)
slopes_cv[[3]] <- sdmTMB:::mround(slopes_cv[[3]], 2)

slopes_cv <- mutate(slopes_cv, text = paste0(b, "\\% (95\\% CI: ", lwr, "--", upr, ")"))
slopes_cv
slopes_cv$measure <- "cv-vs-mare"
slopes <- bind_rows(slopes, slopes_cv)

# save it ------------------

saveRDS(slopes, "data-generated/metrics-slopes-table.rds")
