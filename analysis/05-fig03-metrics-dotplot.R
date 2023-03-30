library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long.rds")
metrics_wide <- readRDS("data-generated/metrics-wide.rds")

m <- select(metrics_wide, survey_abbrev, species_common_name, est_type, cv_orig, prop_mpa) |>
  filter(est_type == "geostat") |>
  distinct()
metrics_long <- distinct(metrics_long)

metrics_long <- metrics_long |> left_join(m) |> ungroup()

metrics_long$species_common_name <- stringr::str_to_title(metrics_long$species_common_name)

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

metrics_long_cochran <- filter(metrics_long, grepl("cochran", est_type))
metrics_long_boot <- filter(metrics_long, grepl("boot", est_type))

metrics_long <- filter(metrics_long, grepl("geo", est_type)) # main
metrics_long$survey_abbrev <- gsub("SYN QCS, SYN HS", "SYN QCS/HS", metrics_long$survey_abbrev)
metrics_long <- filter(metrics_long, survey_abbrev != "SYN QCS/HS")

prop_mpa <- select(metrics_long, survey_abbrev, species_common_name, prop_mpa) |>
  distinct()
# prop_mpa$survey_abbrev <- gsub("SYN QCS/HS", "SYN QCS", prop_mpa$survey_abbrev)
# hs_hack <- filter(prop_mpa, survey_abbrev == "SYN QCS") |> mutate(survey_abbrev = "SYN HS")
# prop_mpa <- bind_rows(prop_mpa, hs_hack) |> distinct()

metrics_long_cochran$prop_mpa <- NULL
metrics_long_boot$prop_mpa <- NULL
metrics_long_cochran <- left_join(metrics_long_cochran, prop_mpa)
metrics_long_boot <- left_join(metrics_long_boot, prop_mpa)

make_dotplot <- function(
    .data, prop_mpa_filter = 0.1, exclude_extrapolation = TRUE,
    colour_var = survey_abbrev, show_prop_mpa = TRUE, dodge_points = FALSE, standard_plot = TRUE, pal = restricted_cols) {
  dat <- .data |> dplyr::filter(prop_mpa > prop_mpa_filter)

  if (exclude_extrapolation) {
    dat <- filter(dat, restr_clean != "Same survey domain")
  } else { # show both types but exlcude some HUGE ones
    dat <- dat |>
      mutate(est = ifelse(grepl("CV", measure) & est >= 500, 500, est)) |>
      mutate(upr = ifelse(grepl("CV", measure) & upr >= 500, 500, upr)) |>
      mutate(lwr = ifelse(grepl("CV", measure) & lwr >= 500, 500, lwr))
  }

  if (standard_plot) {
    dat <- dat |>
      mutate(abs_est_avg = abs(est_avg)) |>
      mutate(survey_abbrev = factor(survey_abbrev,
        levels = c(
          "SYN WCHG",
          "HBLL OUT N",
          "SYN QCS/HS"
        )
      )) |>
      arrange(survey_abbrev, prop_mpa, species_common_name)
  }

  if (show_prop_mpa) {
    dat <- dat |> group_by(species_common_name, survey_abbrev) |>
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

  g <- ggplot(dat, aes(
    x = forcats::fct_inorder(spp_lab_plot),
    y = est,
    ymin = lwr, ymax = upr,
    colour = {{ colour_var }}
  )) +
    geom_hline(yintercept = 0, lty = 2, col = "grey60")

  if (exclude_extrapolation && !dodge_points) {
    g <- g + geom_linerange() +
      geom_point(pch = 21, size = 1.8) +
      geom_point(pch = 20, size = 1.8, alpha = 0.2)
  } else {
    g <- g +
      geom_linerange(position = position_dodge(width = 0.4)) +
      geom_point(position = position_dodge(width = 0.4), pch = 21, size = 1.8) +
      geom_point(position = position_dodge(width = 0.4), pch = 20, size = 1.8, alpha = 0.2)
  }

  g <- g +
    scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0, 0)) +
    coord_flip() +
    scale_colour_manual(values = pal) +
    labs(x = "", y = "", colour = "Survey", fill = "Survey") +
    guides(colour = "none", fill = "none") +
    theme(
      strip.text = element_text(colour = "black"),
      legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"),
      strip.placement = "outside",
      axis.title = element_blank()
    ) +
    facet_grid(survey_abbrev ~ measure,
      scales = "free",
      space = "free_y", switch = "x"
    )
  # g
  # not sure why but egg didn't work here
  # devtools::install_github("eliocamp/tagger")
  g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel", tag_pool = letters[c(1, 4, 7, 2, 5, 8, 3, 6, 9)])
  # print(g)
  g
}

# make_dotplot(metrics_long)
# ggsave("figs/metric-dotplot.pdf", width = 7.8, height = 7.3)

metrics_long |>
  arrange(survey_abbrev, prop_mpa) |>
  make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
    tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
ggsave("figs/metric-dotplot2.pdf", width = 7.8, height = 8.4)

mm <- metrics_long |> bind_rows(metrics_long_boot)
mm |>
  arrange(survey_abbrev, prop_mpa) |>
  make_dotplot(prop_mpa_filter = 0.1, standard_plot = FALSE,
    colour_var = est_type, dodge_points = TRUE) +
    scale_colour_brewer(palette = "Set1", labels = c("Design-based bootstrap", "Geostatistical")) +
    guides(colour = guide_legend()) +
    labs(colour = "Scenario") +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
    tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
ggsave("figs/metric-geo-vs-boot.pdf", width = 7.8, height = 8.4)

# make_dotplot(metrics_long, prop_mpa_filter = 0, exclude_extrapolation = FALSE, colour_var = restr_clean) +
#   scale_colour_brewer(palette = "Set1", labels = c("Extrapolate into MPAs", "Shrink survey domain")) +
#   guides(colour = guide_legend()) +
#   labs(colour = "Scenario")
# ggsave("figs/metric-dotplot-extrapolate.pdf", width = 7.8, height = 11)

# # Compare to worst case scenario with 'as-is-where-is'
# mla <- readRDS("data-generated-ALL/metrics-long.rds")
# mla_wide <- readRDS("data-generated-ALL/metrics-wide.rds")
# m <- select(mla_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
#   distinct()
# mla <- mla |> left_join(m, by = join_by(species_common_name, survey_abbrev))
#
# mla <- filter(mla, restr_clean == "Shrunk survey domain")
# # mla <- filter(mla, prop_mpa > 0.1)
# mla <- mla |> mutate(restr_clean = "Shrunk ALL")
# mla$survey_abbrev <- gsub("SYN QCS, SYN HS", "SYN QCS/HS", mla$survey_abbrev)
#
# metrics_long_keep <- semi_join(metrics_long, select(mla, survey_abbrev, species_common_name))
# .dat <- bind_rows(metrics_long_keep, mla)
#
# .dat <- filter(.dat, !(survey_abbrev == "SYN WCHG" & restr_clean == "Shrunk ALL")) # same!
#
# make_dotplot(.dat, prop_mpa_filter = 0, exclude_extrapolation = TRUE, colour_var = restr_clean, show_prop_mpa = TRUE, dodge_points = TRUE) +
#   scale_colour_brewer(palette = "Set1", labels = c("Restrict within all existing MPAs", "Restrict only within Cat. 1 + 2 MPAs")) +
#   guides(colour = guide_legend()) +
#   labs(colour = "Scenario")
#
# ggsave("figs/metric-dotplot-all-vs-cat12.pdf", width = 7.8, height = 11)

# design-cochran
metrics_long_cochran |>
  arrange(survey_abbrev, prop_mpa) |>
  make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
    tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
ggsave("figs/metric-dotplot-cochran.pdf", width = 7.8, height = 7.8)

metrics_long_boot |>
  arrange(survey_abbrev, prop_mpa) |>
  make_dotplot(standard_plot = FALSE, pal = RColorBrewer::brewer.pal(4, "Dark2")) +
  tagger::tag_facets(tag_prefix = "(", position = "tl", tag = "panel",
    tag_pool = letters[c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)])
ggsave("figs/metric-dotplot-boot.pdf", width = 7.8, height = 7.8)

group_by(metrics_long_cochran, survey_abbrev, measure) |>
  summarise(m = mean(est, na.rm = TRUE))

group_by(metrics_long_boot, survey_abbrev, measure) |>
  summarise(m = mean(est, na.rm = TRUE))

group_by(metrics_long, survey_abbrev, measure) |>
  summarise(m = mean(est, na.rm = TRUE))

# ------------------
# check diff between geo/design:

geo <- filter(metrics_wide, est_type == "geostat", `Restriction type` == "re_shrunk") |>
  select(survey_abbrev, species_common_name, cv_orig_geo = cv_orig) |>
  filter(survey_abbrev != "SYN QCS, SYN HS") |>
  distinct()
bootstrap <- filter(metrics_wide, est_type == "bootstrap", `Restriction type` == "re_shrunk") |>
  select(survey_abbrev, species_common_name, cv_orig_boot = cv_orig) |> distinct()
together <- left_join(geo, bootstrap)
together <- together[together$cv_orig_geo < 0.6, ]
# plot(together$cv_orig_geo, together$cv_orig_boot);abline(0, 1)
mean((together$cv_orig_geo - together$cv_orig_boot) / together$cv_orig_boot * 100, na.rm = TRUE)
cat("percent lower CV")

ggplot(together, aes(cv_orig_boot, cv_orig_geo, colour = survey_abbrev)) +
  geom_linerange(aes(x = cv_orig_boot, ymin = cv_orig_geo, ymax = cv_orig_boot), colour = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
  geom_point(na.rm = TRUE, pch = 21, fill = "white") +
  geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  xlab("CV design-based bootstrap") + ylab("CV geostatistical model") +
  coord_fixed() +
  labs(colour = "Survey", fill = "Survey")
ggsave("figs/boot-vs-geo-cv.pdf", width = 5, height = 4)

compare_geo_boot <- function(compare_col, restr_type = "re_shrunk") {

  geo <- filter(metrics_wide, est_type == "geostat", `Restriction type` == restr_type) |>
    select(survey_abbrev, species_common_name, cv_orig_geo = cv_orig, geo = {{compare_col}}) |>
    filter(survey_abbrev != "SYN QCS, SYN HS") |>
    distinct()
  bootstrap <- filter(metrics_wide, est_type == "bootstrap", `Restriction type` == restr_type) |>
    select(survey_abbrev, species_common_name, design = {{compare_col}}) |> distinct()
  together <- left_join(geo, bootstrap)
  together <- together[together$cv_orig_geo < 0.6, ] # FIXME!?

  ggplot(together, aes(design, geo, colour = survey_abbrev)) +
    geom_linerange(aes(x = design, ymin = geo, ymax = design), colour = "grey70", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 0.5) +
    geom_point(na.rm = TRUE, pch = 21, fill = "white") +
    geom_point(na.rm = TRUE, pch = 21, mapping = aes(fill = survey_abbrev), alpha = 0.18) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    xlab("Design-based") + ylab("Geostatistical model") +
    coord_fixed() +
    labs(colour = "Survey", fill = "Survey") +
    ggrepel::geom_text_repel(aes(label = species_common_name),
      show.legend = FALSE, size = 3)
}

metrics_wide <- mutate(metrics_wide, abs_slope_re = abs(slope_re))

compare_geo_boot(slope_re)
compare_geo_boot(abs_slope_re)
compare_geo_boot(cv_orig)
compare_geo_boot(mare)
compare_geo_boot(cv_ratio) + coord_fixed(ylim = c(0.8, 1.6))


# ggsave("figs/boot-vs-geo-cv.pdf", width = 5, height = 4)
