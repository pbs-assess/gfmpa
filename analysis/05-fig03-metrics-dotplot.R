library(dplyr)
library(ggplot2)
source("analysis/theme.R")
metrics_long <- readRDS("data-generated/metrics-long.rds")
metrics_wide <- readRDS("data-generated/metrics-wide.rds")

m <- select(metrics_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()

metrics_long <- metrics_long |> left_join(m, by = join_by(species_common_name, survey_abbrev))

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

metrics_long$survey_abbrev <- gsub("SYN QCS, SYN HS", "SYN QCS/HS", metrics_long$survey_abbrev)

make_dotplot <- function(
    .data, prop_mpa_filter = 0.1, exclude_extrapolation = TRUE,
    colour_var = survey_abbrev, show_prop_mpa = TRUE, dodge_points = FALSE) {
  dat <- .data |> dplyr::filter(prop_mpa > prop_mpa_filter)

  if (exclude_extrapolation) {
    dat <- filter(dat, restr_clean != "Same survey domain")
  } else { # show both types but exlcude some HUGE ones
    dat <- dat |>
      mutate(est = ifelse(grepl("CV", measure) & est >= 5, 5, est)) |>
      mutate(upr = ifelse(grepl("CV", measure) & upr >= 5, 5, upr)) |>
      mutate(lwr = ifelse(grepl("CV", measure) & lwr >= 5, 5, lwr))
  }

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

  if (show_prop_mpa) {
    dat <- mutate(dat,
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
    scale_colour_manual(values = restricted_cols) +
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

make_dotplot(metrics_long)
ggsave("figs/metric-dotplot.pdf", width = 7.8, height = 7.3)

make_dotplot(metrics_long, prop_mpa_filter = 0, exclude_extrapolation = FALSE, colour_var = restr_clean) +
  scale_colour_brewer(palette = "Set1", labels = c("Extrapolate into MPAs", "Shrink survey domain")) +
  guides(colour = guide_legend()) +
  labs(colour = "Scenario")
ggsave("figs/metric-dotplot-extrapolate.pdf", width = 7.8, height = 11)

# Compare to worst case scenario with 'as-is-where-is'
mla <- readRDS("data-generated-ALL/metrics-long.rds")
mla_wide <- readRDS("data-generated-ALL/metrics-wide.rds")
m <- select(mla_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()
mla <- mla |> left_join(m, by = join_by(species_common_name, survey_abbrev))

mla <- filter(mla, restr_clean == "Shrunk survey domain")
# mla <- filter(mla, prop_mpa > 0.1)
mla <- mla |> mutate(restr_clean = "Shrunk ALL")
mla$survey_abbrev <- gsub("SYN QCS, SYN HS", "SYN QCS/HS", mla$survey_abbrev)

metrics_long_keep <- semi_join(metrics_long, select(mla, survey_abbrev, species_common_name))
.dat <- bind_rows(metrics_long_keep, mla)

.dat <- filter(.dat, !(survey_abbrev == "SYN WCHG" & restr_clean == "Shrunk ALL")) # same!

make_dotplot(.dat, prop_mpa_filter = 0, exclude_extrapolation = TRUE, colour_var = restr_clean, show_prop_mpa = TRUE, dodge_points = TRUE) +
  scale_colour_brewer(palette = "Set1", labels = c("Restrict within all existing MPAs", "Restrict only within Cat. 1 + 2 MPAs")) +
  guides(colour = guide_legend()) +
  labs(colour = "Scenario")

ggsave("figs/metric-dotplot-all-vs-cat12.pdf", width = 7.8, height = 11)
