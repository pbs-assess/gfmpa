library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long2.rds")
metrics_wide <- readRDS("data-generated/metrics-wide2.rds")

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

make_plot <- function(dat, colour_var = survey_abbrev, prop_threshold = 0.15, cv_upper_limit = 100, group_by_survey = TRUE, include_mean_lines = TRUE, group_var = NULL) {

  dat <- filter(dat, prop_mpa >= prop_threshold)

  # get averages...
  temp <- dat |>
    mutate(se = abs((upr - lwr)) / 4) |>
    mutate(variance = se^2)
  if (group_by_survey) {
    temp <- group_by(temp, survey_abbrev, measure, measure_clean, {{ colour_var }})
  } else {
    temp <- group_by(temp, measure, measure_clean, {{ colour_var }})
  }
  # temp <-
    # summarise(mean_est = weighted.mean(est, w = 1 / variance))
   means <- temp |>
    summarise(mean_est = median(est, na.rm = TRUE)) |>
     filter(measure != "slope_re")

  dat <- dat |>
    # hack to cut off upper limit:
    mutate(orig_est = est) |>
    mutate(est = ifelse(grepl("cv", measure) & est >= cv_upper_limit, NA_real_, est)) |>
    mutate(upr = ifelse(grepl("cv", measure) & upr >= cv_upper_limit & lwr <= cv_upper_limit, cv_upper_limit, upr)) |>
    mutate(upr = ifelse(grepl("cv", measure) & upr >= cv_upper_limit & lwr >= cv_upper_limit, NA_real_, upr)) |>
    mutate(lwr = ifelse(grepl("cv", measure) & lwr >= cv_upper_limit, NA_real_, lwr)) |> # put dot on right if removed:
    mutate(est = ifelse(!is.na(orig_est) & is.na(est), cv_upper_limit, est))

  show_prop_mpa <- TRUE
  if (show_prop_mpa) {
    dat <- dat |>
      group_by(species_common_name, survey_abbrev) |>
      mutate(
        spp_lab_plot =
          paste0(
            stringr::str_to_title(species_common_name), " (",
            mround(max(prop_mpa), 2), ")"
          )
      )
  } else {
    dat <- mutate(dat,
      spp_lab_plot =
        stringr::str_to_title(species_common_name)
    )
  }

  dat <- dat |>
    ungroup() |>
    arrange(survey_abbrev, prop_mpa, species_common_name)

  g <- dat |>
    ggplot(aes(
      x = forcats::fct_inorder(paste(survey_abbrev, spp_lab_plot, sep = "-")), # hack!
      y = est,
      ymin = lwr, ymax = upr,
      colour = {{ colour_var }}
    )) +
    geom_hline(yintercept = 0, lty = 1, col = "grey65") +
    scale_x_discrete(labels = function(x) gsub("[a-zA-Z ,]+-", "", x)) # hack!

  g <- g +
    geom_linerange(position = position_dodge(width = 0.5), na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.5), pch = 21, size = 1.8, na.rm = TRUE) +
    geom_point(position = position_dodge(width = 0.5), pch = 20, size = 1.8, alpha = 0.2, na.rm = TRUE)
  # }

  g <- g +
    scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0, 0)) +
    coord_flip() +
    theme(
      strip.text = element_text(colour = "black"),
      legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"),
      strip.placement = "outside",
      axis.title = element_blank(),
    ) +
    facet_grid(survey_abbrev ~ measure_clean,
      scales = "free",
      space = "free_y", switch = "x"
    ) +
    scale_colour_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2")

  if (include_mean_lines)
    g <- g + geom_hline(data = means, mapping = aes(yintercept = mean_est, colour = {{ colour_var }}), lty = 2)
  g + tagger::tag_facets(tag_prefix = "(", position = "tl",
    tag_pool = letters[c(1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12)])
}


.dat <- filter(metrics_long, est_type %in% c("geostat", "bootstrap")) |>
  mutate(est_type = gsub("bootstrap", "Bootstrap", est_type)) |>
  mutate(est_type = gsub("geostat", "Geostatistical", est_type)) |>
  filter(measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk", "Status quo")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS")

temp_mean_lines <- group_by(.dat, survey_abbrev, measure, measure_clean, est_type, type) |>
  summarise(mean_est = median(est, na.rm = TRUE))

g <- .dat |> make_plot(colour_var = type, group_by_survey = FALSE, prop_threshold = 0.15, group_var = survey_abbrev, include_mean_lines = FALSE) +
  coord_flip(ylim = c(0, 1)) +
  geom_hline(yintercept = 0.2, colour = "grey30") +
  facet_grid(
    survey_abbrev ~ est_type,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(colour = "Scenario", fill = "Scenario")
g + geom_hline(data = temp_mean_lines,
  mapping = aes(yintercept = mean_est, colour = type), lty = 2) +
  tagger::tag_facets(tag_prefix = "(", position = "tl",
    tag_pool = letters)
ggsave("figs/abs-cv.pdf", width = 6.5, height = 8.4)
ggsave("figs/abs-cv.png", width = 6.5, height = 8.4)

# restricted_cols <- RColorBrewer::brewer.pal(4, "Set2")
g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  # filter(survey_abbrev != "SYN QCS") |>
  # filter(survey_abbrev != "SYN HS") |>
  mutate(upr = ifelse(measure == "slope_re" & upr > 0.5, 0.5, upr)) |> # FIXME Note
  make_plot(colour_var = survey_abbrev, cv_upper_limit = 100, group_by_survey = FALSE) +
  labs(colour = "Survey", fill = "Survey") +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = restricted_cols) +
  scale_fill_manual(values = restricted_cols)
print(g)
ggsave("figs/metrics-dotplot-main.pdf", width = 7, height = 7.5)
ggsave("figs/metrics-dotplot-main.png", width = 7, height = 7.5)

g <- filter(metrics_long, est_type %in% c("bootstrap")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa_set)) |>
  mutate(prop_mpa = prop_mpa_set) |>  #<
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  # filter(!is.na(prop_pos_set)) |>
  # filter(prop_pos_set > 0.1) |> # TODO note this
  # filter(survey_abbrev != "SYN QCS") |>
  # filter(survey_abbrev != "SYN HS") |>
  # mutate(upr = ifelse(measure == "slope_re" & upr > 0.5, 0.5, upr)) |> # FIXME Note
  mutate(upr = ifelse(measure == "mare" & upr > .83, .83, upr)) |> # FIXME Note
  make_plot(colour_var = survey_abbrev, cv_upper_limit = 1e6, group_by_survey = FALSE) +
  labs(colour = "Survey", fill = "Survey") +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = restricted_cols) +
  scale_fill_manual(values = restricted_cols)
print(g)
ggsave("figs/metrics-dotplot-main-design.pdf", width = 7, height = 7.5)
ggsave("figs/metrics-dotplot-main-design.png", width = 7, height = 7.5)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk", "Restricted")) |>
  mutate(type = case_when(
    type == "Restricted and shrunk" ~ "Restricted and shrunk survey domain",
    type == "Restricted" ~ "Restricted and extrapolated into MPAs",
  )) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS") |>
  filter(survey_abbrev != "SYN HS") |>
  make_plot(colour_var = type, cv_upper_limit = 200, group_by_survey = TRUE) +
  labs(colour = "Scenario", fill = "Scenario") +
  scale_colour_brewer(palette = "Set1")
print(g)
ggsave("figs/metrics-dotplot-extrapolate.pdf", width = 8.5, height = 9)
ggsave("figs/metrics-dotplot-extrapolate.png", width = 8.5, height = 9)

g <- filter(metrics_long, est_type %in% c("geostat", "bootstrap")) |>
  mutate(est_type = case_when(
    est_type == "geostat" ~ "Geostatistical",
    est_type == "bootstrap" ~ "Design-based bootstrap",
  )) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_plot(colour_var = est_type) +
  labs(colour = "Approach", fill = "Approach") +
  scale_colour_brewer(palette = "Set1")
print(g)
ggsave("figs/metrics-dotplot-design-geo.pdf", width = 8.5, height = 9)
ggsave("figs/metrics-dotplot-design-geo.png", width = 8.5, height = 9)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  mutate(upr = ifelse(measure == "slope_re" & upr > 0.8, 0.8, upr)) |> # FIXME Note
  filter(type %in% c("Restricted and shrunk", "Random up-sampled and shrunk 1")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_plot(colour_var = type) +
  labs(colour = "Scenario")
print(g)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk", "Random down-sampled 2")) |>
  filter(!is.na(prop_mpa)) |>
  mutate(upr = if_else(measure == "mare" & upr > 1, 1, upr)) |>
  mutate(upr = if_else(measure == "slope_re" & upr > 1, 1, upr)) |>
  mutate(est = if_else(measure == "slope_re" & est > 1, NA_real_, est)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_plot(colour_var = type) +
  labs(colour = "Scenario")
print(g)

# which models?? ----------

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
m <- bind_rows(list(hbll, syn))
m <- select(m, species_common_name, survey_abbrev, spatiotemporal, family) |>
  distinct() |>
  mutate(species_common_name = stringr::str_to_title(species_common_name)) |>
  mutate(species_common_name = gsub(
    "Rougheye/Blackspotted Rockfish Complex",
    "Rougheye/Blackspotted Rockfish", species_common_name
  )) |>
  mutate(model = paste(stringr::str_to_title(family), paste0("(", spatiotemporal, ")"), sep = "\n")) |>
  mutate(model = gsub("iid", "IID", model)) |>
  mutate(model = gsub("off", "Off", model))
unique(m$model)

m$model <- factor(m$model, levels =
  c(
    "Binomial-Gamma\n(IID, IID)",
    "Binomial-Gamma\n(Off, IID)",
    "Tweedie\n(IID)",
    "Binomial-Gamma\n(Off, Off)"
))

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  left_join(m) |>
  # filter(survey_abbrev != "SYN QCS") |>
  # filter(survey_abbrev != "SYN HS") |>
  mutate(upr = ifelse(measure == "slope_re" & upr > 0.5, 0.5, upr)) |> # FIXME Note
  mutate(survey_abbrev = factor(survey_abbrev,
    levels =
      c("SYN WCHG", "HBLL OUT N", "SYN HS", "SYN QCS", "SYN QCS, SYN HS")
  )) |>
  make_plot(colour_var = model, cv_upper_limit = 100, group_by_survey = FALSE, include_mean_lines = FALSE, prop_threshold = 0.1) +
  labs(colour = "Model", fill = "Model") +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
print(g)
ggsave("figs/metrics-dotplot-by-model.pdf", width = 7.8, height = 9)
ggsave("figs/metrics-dotplot-by-model.png", width = 7.8, height = 9)


# ggplot(metrics, aes(cv_med, geo, colour = survey_abbrev)) +
#   geom_linerange(aes(x = design, ymin = geo, ymax = design), colour = "grey70", alpha = 0.5) +
#   geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
#   geom_point(na.rm = TRUE, pch = 21, fill = "white") +
#   geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
#   scale_colour_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette = "Dark2") +
#   xlab("Design-based") + ylab("Geostatistical model") +
#   coord_fixed() +
#   labs(colour = "Survey", fill = "Survey") +
#   ggrepel::geom_text_repel(aes(label = species_common_name),
#     show.legend = FALSE, size = 2.5, alpha = 0.6) +
#   theme(legend.position = c(0.2,0.8))
# # theme(legend.position = "bottom", legend.direction = "vet")

make_cross_plot <- function(dat1, dat2, xlab = "x", ylab = "y") {
  dat_a <- dat1 |>
    select(species_common_name, survey_abbrev, est, lwr, upr) |>
    rename(est_a = est, lwr_a = lwr, upr_a = upr)
  dat_b <- dat2 |>
    select(species_common_name, survey_abbrev, est, lwr, upr) |>
    rename(est_b = est, lwr_b = lwr, upr_b = upr)
  dat <- left_join(dat_a, dat_b)
  ggplot(dat, aes(est_a, est_b, colour = survey_abbrev)) +
    geom_point(na.rm = TRUE, pch = 21, fill = "white") +
    geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
    xlab(xlab) +
    ylab(ylab) +
    coord_fixed() +
    ggrepel::geom_text_repel(aes(label = species_common_name),
      show.legend = FALSE, size = 2.5, alpha = 0.6
    ) +
    geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
    geom_linerange(aes(x = est_a, ymin = est_a, ymax = est_b), colour = "grey70", alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(colour = "Survey", fill = "Survey") +
    theme(legend.position = c(0.2, 0.8))
}

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Random up-sampled and shrunk 1")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("CV")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Random down-sampled 1")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random down-sampled 1") +
  ggtitle("CV")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Random down-sampled 2")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random down-sampled 2") +
  ggtitle("MARE")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Random down-sampled 2")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random down-sampled 2") +
  ggtitle("slope_re")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Random up-sampled and shrunk 1")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("MARE")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Random up-sampled and shrunk 1")
make_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("Absolute relative error slope")

upper <- 0.9
d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Restricted and shrunk") |>
  filter(survey_abbrev != "SYN QCS, SYN HS")
d2 <- metrics_long |>
  filter(est_type == "bootstrap", measure == "cv", type == "Restricted and shrunk")
g1 <- make_cross_plot(d2, d1, "bootstrap", "geostat") +
  labs(y = "Geostatistical CV", x = "Design-based CV") +
  coord_fixed(xlim = c(0.1, upper), ylim = c(0.1, upper)) +
  ggtitle("Restricted and shrunk")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Status quo") |>
  filter(survey_abbrev != "SYN QCS, SYN HS")
d2 <- metrics_long |>
  filter(est_type == "bootstrap", measure == "cv", type == "Status quo")
g2 <- make_cross_plot(d2, d1, "bootstrap", "geostat") +
  labs(y = "Geostatistical CV", x = "Design-based CV") +
  coord_fixed(xlim = c(0.1, upper), ylim = c(0.1, upper)) +
  ggtitle("Status quo")
cowplot::plot_grid(g2, g1, nrow = 1)
ggsave("figs/geostat-vs-design2.pdf", width = 8, height = 4.5)
ggsave("figs/geostat-vs-design2.png", width = 8, height = 4.5)

metrics_wide |>
  filter(type != "MPA only") |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  filter(!grepl("cochran", est_type)) |>
  filter(!grepl("^Restricted$", type)) |>
  mutate(type = gsub(" [0-9]$", "", type)) |>
  group_by(type, est_type) |>
  filter(est_type == "geostat") |>
  summarise(
    mean_cv = mean(cv_med, na.rm = TRUE),
    mean_cv_change = mean(cv_perc_med, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_mare = mean(mare_med, na.rm = TRUE),
    mean_slope = mean(abs(slope_re_med), na.rm = TRUE)
  ) |>
  select(-mean_cv) |>
  arrange(mean_cv_change) |>
  knitr::kable(digits = 2L)

# metrics_wide |>
#   filter(est_type == "geostat") |>
#   group_by(survey_abbrev, species_common_name) |>
#   summarise()
#

extra <- select(
  metrics_wide, species_common_name, survey_abbrev, type,
  est_type, coverage, prop_mpa
) |>
  distinct() |>
  rename(est = coverage) |>
  mutate(measure = "coverage", measure_clean = "50% CI coverage")

to_plot <- metrics_long |>
  mutate(measure_clean = ifelse(measure == "cv", "CV", measure_clean)) |>
  bind_rows(extra) |>
  filter(!is.na(prop_mpa)) |>
  filter(prop_mpa > 0.1) |>
  filter(type != "MPA only") |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  # filter(survey_abbrev != "SYN QCS") |>
  # filter(survey_abbrev != "SYN HS") |>
  filter(!grepl("cochran", est_type)) |>
  filter(!grepl("^Restricted$", type)) |>
  mutate(type2 = gsub(" [0-9]$", "", type)) |>
  group_by(species_common_name, survey_abbrev, est_type, measure, measure_clean, type2) |>
  summarise(est = mean(est)) |> # compress sims
  mutate(est_diff = est / est[type2 == "Restricted and shrunk"]) |>
  # mutate(type = gsub(" [0-9]$", "", type)) |>

  filter(est_type == "geostat") |>
  ungroup() |>
  mutate(est = if_else(measure =="slope_re", abs(est), est)) # abs() it

th <- ggsidekick::theme_sleek() +
  theme(
    panel.grid = element_line(colour = "grey90"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25))
  )

# to_plot |>
#   filter(type2 != "Restricted and shrunk") |>
#   filter(measure != "cv_perc") |>
#   ggplot(aes(type2, est_diff, colour = survey_abbrev)) +
#   geom_boxplot() +
#   facet_wrap(~measure_clean, scales = "free_x") +
#   scale_y_log10() +
#   geom_hline(yintercept = 1, lty = 2) +
#   coord_flip() +
#   scale_color_manual(values = restricted_cols) +
#   th

to_plot <- to_plot |>
  filter(measure != "cv") |>
  filter(type2 != "Status quo and shrunk") |>
  filter(type2 != "Status quo")

all_conv <- to_plot |>
  group_by(survey_abbrev, species_common_name) |>
  summarise(n = length(unique(type2))) |>
  filter(n == 3L) |>
  select(-n)

pmpa <- select(metrics_wide, survey_abbrev, species_common_name, prop_mpa) |> distinct() |>
  filter(!is.na(prop_mpa))

to_plot <- left_join(to_plot, pmpa)

semi_join(to_plot, all_conv) |>
  group_by(survey_abbrev, measure_clean, type2) |>
  summarise(
    lwr = quantile(est, probs = 0.25),
    upr = quantile(est, probs = 0.75),
    est = median(est),
    mean_prop = mean(prop_mpa, na.rm = TRUE),
  ) |>
  ungroup() |>
  mutate(survey_abbrev = forcats::fct_reorder(survey_abbrev, -mean_prop)) |>
  ggplot(aes(type2, est, colour = type2)) +
  # geom_boxplot() +
  geom_pointrange(aes(ymin = lwr, ymax = upr), pch = 21, size = 0.3) +
  facet_grid(survey_abbrev~measure_clean, scales = "free_x") +
  # facet_grid(~measure_clean, scales = "free_x") +
  # scale_y_log10() +
  coord_flip() +
  # geom_jitter(height = 0, width = 0.3, mapping = aes(colour = survey_abbrev), pch = 21, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Set2")[1:3])) +
  th +
  theme(axis.title.y = element_blank()) +
  labs(y = "Metric value", colour = "Scenario") +
  guides(colour = "none") +
  theme(strip.text.y.right = element_text(size = 7))
ggsave("figs/sampled-dotplot-comparison.pdf", width = 7, height = 3.4)
ggsave("figs/sampled-dotplot-comparison.png", width = 7, height = 3.4)

# plot differently? -------
pal <- as.character(colorBlindness::availableColors())
semi_join(to_plot, all_conv) |>
  group_by(survey_abbrev, measure_clean, type2) |>
  summarise(
    lwr = quantile(est, probs = 0.25),
    upr = quantile(est, probs = 0.75),
    est = median(est),
    mean_prop = mean(prop_mpa, na.rm = TRUE),
  ) |>
  ungroup() |>
  mutate(survey_abbrev = forcats::fct_reorder(survey_abbrev, -mean_prop)) |>
  ggplot(aes(survey_abbrev, est, colour = type2)) +
  # geom_boxplot() +
  geom_pointrange(aes(ymin = lwr, ymax = upr), pch = 21, size = 0.3, position = position_dodge(width = 0.38)) +
  facet_wrap(~measure_clean, scales = "free_x", ncol = 4) +
  # facet_grid(~measure_clean, scales = "free_x") +
  # scale_y_log10() +
  coord_flip() +
  # geom_jitter(height = 0, width = 0.3, mapping = aes(colour = survey_abbrev), pch = 21, alpha = 0.5) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = c("Restricted and shrunk" = "grey30", "Random up-sampled and shrunk" = pal[3], "Random down-sampled" = pal[7])) +
  # scale_color_manual(values = c(RColorBrewer::brewer.pal(4, "Set2")[1:3])) +
  th +
  theme(axis.title.y = element_blank(), legend.position = "top") +
  labs(y = "Metric value", colour = "Scenario") +
  # guides(colour = "none") +
  theme(strip.text.y.right = element_text(size = 7)) +
  tagger::tag_facets(tag_prefix = "(", position = "tl")
ggsave("figs/sampled-dotplot-comparison2.pdf", width = 7, height = 3)
ggsave("figs/sampled-dotplot-comparison2.png", width = 7, height = 3)

# NOTES:
# downsampling helps trend bias a bit
# MARE about same as MPA + shrunk...
# should be random downsample and NOT shrunk


# all_conv <- to_plot |>
#   group_by(survey_abbrev, species_common_name) |>
#   summarise(n = length(unique(type2))) |>
#   filter(n == 3L) |>
#   select(-n)


# FIXME: drop those that don't converge in one or other!!
all_conv <- metrics_long |>
  filter(est_type %in% c("geostat", "bootstrap"), !grepl("sample", type)) |>
  group_by(survey_abbrev, species_common_name) |>
  summarise(n = length(unique(est_type))) |>
  filter(n == 2L) |>
  select(-n)


metrics_long |>
  semi_join(all_conv) |>
  left_join(pmpa) |>
  filter(est_type %in% c("geostat", "bootstrap"), !grepl("sample", type)) |> glimpse() |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  mutate(est = if_else(measure =="slope_re", abs(est), est)) |>  # abs() it
filter(prop_mpa > 0.1) |>
  filter(measure != "cv") |>
  filter(measure != "cv_perc") |>
  ungroup() |>
  group_by(survey_abbrev, est_type, measure_clean) |>
  summarise(
    lwr = quantile(est, probs = 0.25),
    upr = quantile(est, probs = 0.75),
    est = mean(est),
    mean_prop = mean(prop_mpa, na.rm = TRUE),
  )  |>
    # mutate(survey_abbrev = forcats::fct_reorder(survey_abbrev, -mean_prop)) |>
  ggplot(aes(est_type, est, colour = est_type)) +
  coord_flip() +
  geom_pointrange(aes(ymin = lwr, ymax = upr), pch = 21) +
  facet_grid(forcats::fct_reorder(survey_abbrev, -mean_prop)~measure_clean, scales = "free_x") +
  xlab("") + ylab("Metric value") + guides(colour = "none") +
  scale_color_brewer(palette = "Set2")

ggsave("figs/geo-design-mare-trend.pdf", width = 5, height = 5)
ggsave("figs/geo-design-mare-trend.png", width = 5, height = 5)

  # summarise(
  #   mean_cv = mean(cv_med, na.rm = TRUE),
  #   mean_cv_change = mean(cv_perc_med, na.rm = TRUE),
  #   mean_coverage = mean(coverage, na.rm = TRUE),
  #   mean_mare = mean(mare_med, na.rm = TRUE),
  #   mean_slope = mean(abs(slope_re_med), na.rm = TRUE)
  # ) |> knitr::kable(digits = 3)

metrics_long |>
  filter(est_type %in% c("geostat", "bootstrap"), !grepl("sample", type)) |> glimpse() |>
  group_by(survey_abbrev, est_type, measure) |>
  summarise(value = mean(est, na.rm = TRUE)) |>
  filter(measure == "mare") |> knitr::kable(digits = 3)
