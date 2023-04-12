library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long2.rds")
metrics_wide <- readRDS("data-generated/metrics-wide2.rds")

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

make_zee_plot <- function(dat, colour_var = survey_abbrev, prop_threshold = 0.1, cv_upper_limit = 100, group_by_survey = TRUE) {

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
    # summarise(mean_est = weighted.mean(est, w = 1 / variance))
   means <- temp |>  summarise(mean_est = median(est))

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

  g <- g + geom_hline(data = means, mapping = aes(yintercept = mean_est, colour = {{ colour_var }}), lty = 2)
  g
}

metrics_long$survey_abbrev <- factor(metrics_long$survey_abbrev,
  levels =
    c("SYN WCHG", "HBLL OUT N", "SYN HS", "SYN QCS", "SYN QCS, SYN HS")
)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS") |>
  filter(survey_abbrev != "SYN HS") |>
  mutate(upr = ifelse(measure == "slope_re" & upr > 0.5, 0.5, upr)) |> # FIXME Note
  make_zee_plot(colour_var = survey_abbrev, cv_upper_limit = 100, group_by_survey = FALSE) +
  labs(colour = "Survey", fill = "Survey") +
  guides(colour = "none", fill = "none") +
  scale_color_manual(values = restricted_cols) +
  scale_fill_manual(values = restricted_cols)
print(g)

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
  make_zee_plot(colour_var = type, cv_upper_limit = 200, group_by_survey = TRUE) +
  labs(colour = "Scenario", fill = "Scenario")
print(g)

g <- filter(metrics_long, est_type %in% c("geostat", "bootstrap")) |>
  mutate(est_type = case_when(
    est_type == "geostat" ~ "Geostatistical",
    est_type == "bootstrap" ~ "Design-based bootstrap",
  )) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_zee_plot(colour_var = est_type) +
  labs(colour = "Approach", fill = "Approach")
print(g)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  mutate(upr = ifelse(measure == "slope_re" & upr > 0.8, 0.8, upr)) |> # FIXME Note
  filter(type %in% c("Restricted and shrunk", "Random up-sampled and shrunk 1")) |>
  filter(!is.na(prop_mpa)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_zee_plot(colour_var = type) +
  labs(colour = "Scenario")
print(g)

g <- filter(metrics_long, est_type %in% c("geostat")) |>
  filter(!measure %in% "cv") |>
  filter(type %in% c("Restricted and shrunk", "Random down-sampled and shrunk 2")) |>
  filter(!is.na(prop_mpa)) |>
  mutate(upr = if_else(measure == "mare" & upr > 1, 1, upr)) |>
  mutate(upr = if_else(measure == "slope_re" & upr > 1, 1, upr)) |>
  mutate(est = if_else(measure == "slope_re" & est > 1, NA_real_, est)) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  make_zee_plot(colour_var = type) +
  labs(colour = "Scenario")
print(g)

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

make_zee_cross_plot <- function(dat1, dat2, xlab = "x", ylab = "y") {
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
make_zee_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("CV")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "mare", type == "Random up-sampled and shrunk 1")
make_zee_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("MARE")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "geostat", measure == "slope_re", type == "Random up-sampled and shrunk 1")
make_zee_cross_plot(d1, d2, "Restricted and shrunk", "Random up-sampled and shrunk 1") +
  ggtitle("Absolute relative error slope")

d1 <- metrics_long |>
  filter(est_type == "geostat", measure == "cv", type == "Restricted and shrunk")
d2 <- metrics_long |>
  filter(est_type == "bootstrap", measure == "cv", type == "Restricted and shrunk")
make_zee_cross_plot(d1, d2, "geostat", "bootstrap") +
  ggtitle("CV") + labs(x = "Geostatistical", y = "Design-based")

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
    mean_slope = mean(abs_slope_re, na.rm = TRUE)
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
  filter(type != "MPA only") |>
  filter(!is.na(prop_mpa)) |>
  # filter(survey_abbrev != "SYN QCS, SYN HS") |>
  # filter(survey_abbrev != "SYN QCS, SYN HS") |>
  filter(survey_abbrev != "SYN QCS") |>
  filter(survey_abbrev != "SYN HS") |>
  filter(!grepl("cochran", est_type)) |>
  filter(!grepl("^Restricted$", type)) |>
  mutate(type2 = gsub(" [0-9]$", "", type)) |>
  group_by(species_common_name, survey_abbrev, est_type, measure, measure_clean, type2) |>
  summarise(est = mean(est)) |> # compress sims
  mutate(est_diff = est / est[type2 == "Restricted and shrunk"]) |>
  # mutate(type = gsub(" [0-9]$", "", type)) |>

  filter(est_type == "geostat")

th <- ggsidekick::theme_sleek() +
  theme(
    panel.grid = element_line(colour = "grey90"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25))
  )

to_plot |>
  filter(type2 != "Restricted and shrunk") |>
  filter(measure != "cv_perc") |>
  ggplot(aes(type2, est_diff, colour = survey_abbrev)) +
  geom_boxplot() +
  facet_wrap(~measure_clean, scales = "free_x") +
  scale_y_log10() +
  geom_hline(yintercept = 1, lty = 2) +
  coord_flip() +
  scale_color_manual(values = restricted_cols) +
  th

to_plot |>
  filter(measure != "cv") |>
  ggplot(aes(type2, est)) +
  geom_boxplot() +
  facet_wrap(~measure_clean, scales = "free_x") +
  # scale_y_log10() +
  coord_flip() +
  geom_jitter(height = 0, width = 0.3, mapping = aes(colour = survey_abbrev), pch = 21, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  # scale_color_manual(values = restricted_cols) +
  th +
  theme(axis.title.y = element_blank()) +
  labs(y = "Metric value", colour = "Survey")

# NOTES:
# downsampling helps trend bias a bit
# MARE about same as MPA + shrunk...
# should be random downsample and NOT shrunk

# theme(axis.text.x = element_text(angle = 90, hjust = 1))

#
# # ----------------------------------
#
#
#
# make_dotplot <- function(
#     .data, prop_mpa_filter = 0.1, exclude_extrapolation = TRUE,
#     colour_var = survey_abbrev, show_prop_mpa = TRUE, dodge_points = FALSE, standard_plot = TRUE, pal = restricted_cols) {
#   dat <- .data |> dplyr::filter(prop_mpa > prop_mpa_filter)
#
#   if (exclude_extrapolation) {
#     dat <- filter(dat, restr_clean != "Same survey domain")
#   } else { # show both types but exlcude some HUGE ones
#     dat <- dat |>
#       mutate(est = ifelse(grepl("CV", measure) & est >= 500, 500, est)) |>
#       mutate(upr = ifelse(grepl("CV", measure) & upr >= 500, 500, upr)) |>
#       mutate(lwr = ifelse(grepl("CV", measure) & lwr >= 500, 500, lwr))
#   }
#
#   if (standard_plot) {
#     dat <- dat |>
#       mutate(abs_est_avg = abs(est_avg)) |>
#       mutate(survey_abbrev = factor(survey_abbrev,
#         levels = c(
#           "SYN WCHG",
#           "HBLL OUT N",
#           "SYN QCS/HS"
#         )
#       )) |>
#       arrange(survey_abbrev, prop_mpa, species_common_name)
#   }
#
#   if (show_prop_mpa) {
#     dat <- dat |> group_by(species_common_name, survey_abbrev) |>
#       mutate(
#       spp_lab_plot =
#         paste0(
#           stringr::str_to_title(species_common_name), " (",
#           mround(max(prop_mpa), 2), ")"
#         )
#     )
#   } else {
#     dat <- mutate(dat,
#       spp_lab_plot =
#         stringr::str_to_title(species_common_name)
#     )
#   }
#
#   g <- ggplot(dat, aes(
#     x = forcats::fct_inorder(spp_lab_plot),
#     y = est,
#     ymin = lwr, ymax = upr,
#     colour = {{ colour_var }}
#   )) +
#     geom_hline(yintercept = 0, lty = 2, col = "grey60")
#
#   if (exclude_extrapolation && !dodge_points) {
#     g <- g + geom_linerange() +
#       geom_point(pch = 21, size = 1.8) +
#       geom_point(pch = 20, size = 1.8, alpha = 0.2)
#   } else {
#     g <- g +
#       geom_linerange(position = position_dodge(width = 0.4)) +
#       geom_point(position = position_dodge(width = 0.4), pch = 21, size = 1.8) +
#       geom_point(position = position_dodge(width = 0.4), pch = 20, size = 1.8, alpha = 0.2)
#   }
#
#   g <- g +
#     scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0, 0)) +
#     coord_flip() +
#     scale_colour_manual(values = pal) +
#     labs(x = "", y = "", colour = "Survey", fill = "Survey") +
#     guides(colour = "none", fill = "none") +
#     theme(
#       strip.text = element_text(colour = "black"),
#       legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"),
#       strip.placement = "outside",
#       axis.title = element_blank()
#     ) +
#     facet_grid(survey_abbrev ~ measure,
#       scales = "free",
#       space = "free_y", switch = "x"
#     )
#   # g
#   # not sure why but egg didn't work here
#   # devtools::install_github("eliocamp/tagger")
#   g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel", tag_pool = letters[c(1, 4, 7, 2, 5, 8, 3, 6, 9)])
#   # print(g)
#   g
# }
#
# # make_dotplot(metrics_long)
# # ggsave("figs/metric-dotplot.pdf", width = 7.8, height = 7.3)
#
# metrics_long |>
#   arrange(survey_abbrev, prop_mpa) |>
#   make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
#   tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
#     tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
# ggsave("figs/metric-dotplot2.pdf", width = 7.8, height = 8.4)
#
# mm <- metrics_long |> bind_rows(metrics_long_boot)
# mm |>
#   arrange(survey_abbrev, prop_mpa) |>
#   make_dotplot(prop_mpa_filter = 0.1, standard_plot = FALSE,
#     colour_var = est_type, dodge_points = TRUE) +
#     scale_colour_brewer(palette = "Set1", labels = c("Design-based bootstrap", "Geostatistical")) +
#     guides(colour = guide_legend()) +
#     labs(colour = "Scenario") +
#   tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
#     tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
# ggsave("figs/metric-geo-vs-boot.pdf", width = 7.8, height = 8.4)
#
# # make_dotplot(metrics_long, prop_mpa_filter = 0, exclude_extrapolation = FALSE, colour_var = restr_clean) +
# #   scale_colour_brewer(palette = "Set1", labels = c("Extrapolate into MPAs", "Shrink survey domain")) +
# #   guides(colour = guide_legend()) +
# #   labs(colour = "Scenario")
# # ggsave("figs/metric-dotplot-extrapolate.pdf", width = 7.8, height = 11)
#
# # # Compare to worst case scenario with 'as-is-where-is'
# # mla <- readRDS("data-generated-ALL/metrics-long.rds")
# # mla_wide <- readRDS("data-generated-ALL/metrics-wide.rds")
# # m <- select(mla_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
# #   distinct()
# # mla <- mla |> left_join(m, by = join_by(species_common_name, survey_abbrev))
# #
# # mla <- filter(mla, restr_clean == "Shrunk survey domain")
# # # mla <- filter(mla, prop_mpa > 0.1)
# # mla <- mla |> mutate(restr_clean = "Shrunk ALL")
# # mla$survey_abbrev <- gsub("SYN QCS, SYN HS", "SYN QCS/HS", mla$survey_abbrev)
# #
# # metrics_long_keep <- semi_join(metrics_long, select(mla, survey_abbrev, species_common_name))
# # .dat <- bind_rows(metrics_long_keep, mla)
# #
# # .dat <- filter(.dat, !(survey_abbrev == "SYN WCHG" & restr_clean == "Shrunk ALL")) # same!
# #
# # make_dotplot(.dat, prop_mpa_filter = 0, exclude_extrapolation = TRUE, colour_var = restr_clean, show_prop_mpa = TRUE, dodge_points = TRUE) +
# #   scale_colour_brewer(palette = "Set1", labels = c("Restrict within all existing MPAs", "Restrict only within Cat. 1 + 2 MPAs")) +
# #   guides(colour = guide_legend()) +
# #   labs(colour = "Scenario")
# #
# # ggsave("figs/metric-dotplot-all-vs-cat12.pdf", width = 7.8, height = 11)
#
# # design-cochran
# metrics_long_cochran |>
#   arrange(survey_abbrev, prop_mpa) |>
#   make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
#   tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
#     tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
# ggsave("figs/metric-dotplot-cochran.pdf", width = 7.8, height = 7.8)
#
# metrics_long_boot |>
#   arrange(survey_abbrev, prop_mpa) |>
#   make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
#   tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
#     tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
# ggsave("figs/metric-dotplot-boot.pdf", width = 7.8, height = 7.8)
#
# group_by(metrics_long_cochran, survey_abbrev, measure) |>
#   summarise(m = mean(est, na.rm = TRUE))
#
# group_by(metrics_long_boot, survey_abbrev, measure) |>
#   summarise(m = mean(est, na.rm = TRUE))
#
# group_by(metrics_long, survey_abbrev, measure) |>
#   summarise(m = mean(est, na.rm = TRUE))
#
# # ------------------
# # check diff between geo/design:
#
# geo <- filter(metrics_wide, est_type == "geostat", `Restriction type` == "re_shrunk") |>
#   select(survey_abbrev, species_common_name, cv_orig_geo = cv_orig) |>
#   filter(survey_abbrev != "SYN QCS, SYN HS") |>
#   distinct()
# bootstrap <- filter(metrics_wide, est_type == "bootstrap", `Restriction type` == "re_shrunk") |>
#   select(survey_abbrev, species_common_name, cv_orig_boot = cv_orig) |> distinct()
# together <- left_join(geo, bootstrap)
# together <- together[together$cv_orig_geo < 0.6, ]
# # plot(together$cv_orig_geo, together$cv_orig_boot);abline(0, 1)
# mean((together$cv_orig_geo - together$cv_orig_boot) / together$cv_orig_boot * 100, na.rm = TRUE)
# cat("percent lower CV")
#
# ggplot(together, aes(cv_orig_boot, cv_orig_geo, colour = survey_abbrev)) +
#   geom_linerange(aes(x = cv_orig_boot, ymin = cv_orig_geo, ymax = cv_orig_boot), colour = "grey70", alpha = 0.5) +
#   geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
#   geom_point(na.rm = TRUE, pch = 21, fill = "white") +
#   geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
#   scale_colour_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette = "Dark2") +
#   xlab("CV design-based bootstrap") + ylab("CV geostatistical model") +
#   coord_fixed() +
#   labs(colour = "Survey", fill = "Survey")
# ggsave("figs/boot-vs-geo-cv.pdf", width = 5, height = 4)
#
# compare_geo_boot <- function(compare_col, dat = metrics_wide, restr_type = "re_shrunk") {
#
#   geo <- filter(dat, est_type == "geostat", `Restriction type` == restr_type) |>
#     select(survey_abbrev, species_common_name, cv_orig_geo = cv_orig, geo = {{compare_col}}) |>
#     filter(survey_abbrev != "SYN QCS, SYN HS") |>
#     distinct()
#   bootstrap <- filter(dat, est_type == "bootstrap", `Restriction type` == restr_type) |>
#     select(survey_abbrev, species_common_name, design = {{compare_col}}) |> distinct()
#   together <- left_join(geo, bootstrap)
#   together <- together[together$cv_orig_geo < 0.6, ] # FIXME!?
#
#   ggplot(together, aes(design, geo, colour = survey_abbrev)) +
#     geom_linerange(aes(x = design, ymin = geo, ymax = design), colour = "grey70", alpha = 0.5) +
#     geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
#     geom_point(na.rm = TRUE, pch = 21, fill = "white") +
#     geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
#     scale_colour_brewer(palette = "Dark2") +
#     scale_fill_brewer(palette = "Dark2") +
#     xlab("Design-based") + ylab("Geostatistical model") +
#     coord_fixed() +
#     labs(colour = "Survey", fill = "Survey") +
#     ggrepel::geom_text_repel(aes(label = species_common_name),
#       show.legend = FALSE, size = 2.5, alpha = 0.6) +
#     theme(legend.position = c(0.2,0.8))
#     # theme(legend.position = "bottom", legend.direction = "vet")
# }
#
#
# ind <- readRDS("data-generated/index-filtered.rds")
# ind$est_type <- "geostat"
# ind_des <- readRDS("data-generated/stratified-random-design-all.rds")
# ind <- bind_rows(ind, ind_des) |>
#   filter(type == "Restricted and shrunk") |>
#   group_by(survey_abbrev, species_common_name, est_type) |>
#   summarise(mean_cv = mean(cv, na.rm = TRUE), cv_orig = mean(orig_cv, na.rm = TRUE)) |>
#   filter(survey_abbrev != "SYN QCS, SYN HS")
# ind
# ind$`Restriction type` <- "re_shrunk"
# ind$species_common_name <- stringr::str_to_title(ind$species_common_name)
# ind$cv_ratio <- ind$mean_cv / ind$cv_orig
#
# metrics_wide <- mutate(metrics_wide, abs_slope_re = abs(slope_re))
#
# # compare_geo_boot(slope_re)
# g1 <- compare_geo_boot(abs_slope_re) +
#   ggtitle("Absolute RE slope")
# g2 <- compare_geo_boot(mare)+
#   ggtitle("MARE")
# # compare_geo_boot(cv_ratio, dat = ind)+ coord_fixed(ylim = c(0.9, 1.5), xlim = c(0.9, 1.5))
# g3 <- compare_geo_boot(cv_ratio) +
#   coord_fixed(ylim = c(0.9, 1.7), xlim = c(0.9, 1.7))+
#   ggtitle("Mean CV ratio")
# g4 <- compare_geo_boot(cv_orig) +
#   coord_fixed(ylim = c(0.08, 0.65), xlim = c(0.08, 0.65))+
#   ggtitle("CV status quo")
# g5 <- compare_geo_boot(mean_cv, dat = ind)+coord_fixed(ylim = c(0.08, 0.65), xlim = c(0.08, 0.65))+
#   ggtitle("CV restricted and shrunk")
#
# cowplot::plot_grid(g3, g2, g1, ncol = 3)
# ggsave("figs/boot-vs-geo-metrics.pdf", width = 12, height = 4)
#
# cowplot::plot_grid(g4, g5, ncol = 2)
# ggsave("figs/boot-vs-geo-cv.pdf", width = 9, height = 4)
#
# # ggsave("figs/boot-vs-geo-cv.pdf", width = 5, height = 4)
