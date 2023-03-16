library(dplyr)
library(ggplot2)
ggplot2::theme_set(ggsidekick::theme_sleek())
metrics_long <- readRDS("data-generated/metrics-long.rds")
metrics_wide <- readRDS("data-generated/metrics-wide.rds")

m <- select(metrics_wide, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()

metrics_long <- metrics_long |> left_join(m, by = join_by(species_common_name, survey_abbrev))

restricted_cols <- RColorBrewer::brewer.pal(4, "Dark2")[-3][c(2, 3, 1)]

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

g <- metrics_long |>
  filter(prop_mpa > 0.1) |>
  filter(restr_clean == "Shrunk survey domain") |>
  # filter(!survey_abbrev %in% c("SYN HS", "SYN WCHG")) |>
  mutate(abs_est_avg = abs(est_avg)) |>
  mutate(survey_abbrev = factor(survey_abbrev,
    levels = c(
      "SYN WCHG",
      "HBLL OUT N",
      "SYN QCS, SYN HS"
    ))) |>
  arrange(survey_abbrev, prop_mpa, species_common_name) |>
  ggplot(aes(
    forcats::fct_inorder(
      paste0(stringr::str_to_title(species_common_name), " (", mround(prop_mpa, 2), ")")
      ),
    est,
    ymin = lwr, ymax = upr,
    colour = survey_abbrev
    # colour = as.factor(restr_clean)
  )) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  # geom_pointrange(position = position_dodge(width = 0.75), size = 0.35, na.rm = TRUE) +

  geom_linerange() +
  geom_point(pch = 21, size = 1.8) +
  geom_point(pch = 20, size = 1.8, alpha = 0.2) +

  coord_flip() +
  scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0, 0)) +
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
g <- g + tagger::tag_facets(tag_prefix = "(", position = "tl")
# print(g)

ggsave("figs/index-geo-combined-dotplot.pdf", width = 7.8, height = 8)
