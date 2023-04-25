library(dplyr)
library(ggplot2)
source("analysis/theme.R")

# index <- readRDS("data-generated/index-filtered.rds")

# index <- readRDS("data-generated/metrics-long2.rds")

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
index <- bind_rows(hbll, syn)

index <- mutate(index, species_common_name = gsub(
  "rougheye/blackspotted rockfish complex",
  "rougheye/blackspotted rockfish", species_common_name
)) |>
  filter(!(year == 2014 & survey_abbrev == "SYN WCHG"))

metrics <- readRDS("data-generated/metrics-wide2.rds")

prop <- metrics |>
  select(
    prop_mpa, species_common_name,
    survey_abbrev
  ) |>
  distinct() |>
  mutate(species_common_name = tolower(species_common_name))

index$prop_mpa <- NULL
index <- left_join(index, prop)

top_prop <- prop |>
  filter(prop_mpa > 0.15) |>
  # filter(!survey_abbrev %in% c("SYN HS", "SYN QCS")) |>
  filter(!survey_abbrev %in% c("SYN QCS, SYN HS")) |>
  group_by(survey_abbrev) |>
  # top_n(n = 8L, wt = prop_mpa) |>
  arrange(survey_abbrev, -prop_mpa) |>
  mutate(prop_mpa = round(prop_mpa, 2))

top_prop |> as.data.frame()

# top_prop |> readr::write_csv("data-raw/selected-spp.csv")

# fig2_keep <- top_prop |> select(-prop_mpa) |> mutate(species_common_name_lower = tolower(species_common_name))
  # filter(!(species_common_name == "Redbanded Rockfish" & survey_abbrev == "HBLL OUT N"))

dodge_width <- 1

colour_pal <- c("gray50", restricted_cols)

# fig2_keep <- bind_rows(fig2_keep, tibble(species_common_name == ""))

# m <- select(metrics, survey_abbrev, species_common_name, prop_mpa) |>
  # distinct()

# index <- index |> left_join(m, by = join_by(species_common_name, survey_abbrev))

fig2_keep <- readr::read_csv("data-raw/selected-spp.csv") |>
  select(-prop_mpa, -notes) |>
  filter(!is.na(include)) |>
  mutate(species_common_name_lower = tolower(species_common_name))
nrow(fig2_keep)

ind <- index |>
  mutate(species_common_name_lower = tolower(species_common_name)) |>
  # semi_join(top_prop, by = join_by(survey_abbrev, species_common_name_lower)) |>
  semi_join(fig2_keep, by = join_by(survey_abbrev, species_common_name_lower))

labeller_fn <- function(labels, multi_line = TRUE) {
  labels <- lapply(labels, function(x) {
    gsub("[a-zA-Z ,]+-", "", as.character(x))
  })
  if (multi_line) {
    labels
  }
  else {
    collapse_labels_lines(labels)
  }
}

cols <- c(RColorBrewer::brewer.pal(4L, "Set2"), "grey50")
g <- ind |>
  # filter(species_common_name_lower == "english sole", survey_abbrev == "SYN WCHG") |>
  mutate(spp_survey = species_common_name) |>
  mutate(spp_survey = gsub("HBLL OUT N", "HBLL", spp_survey)) |>
  filter(type %in% c("Status quo", "Restricted and shrunk")) |>
  group_by(spp_survey, type) |>
  mutate(
    colour_var =
      ifelse(
        type == "Status quo",
        paste0("Status quo ", survey_abbrev),
        paste0("Restricted")
  )) |>
  mutate(colour_var = factor(colour_var,
    # levels = c(
    #   "Status quo",
    #   "Restricted SYN WCHG",
    #   "Restricted HBLL OUT N",
    #   # "Restricted SYN QCS, SYN HS"
    #   "Restricted SYN HS",
    #   "Restricted SYN QCS"
    # )
    levels = c(
      "Status quo",
      "Status quo SYN WCHG",
      "Status quo HBLL OUT N",
      # "Status quo SYN QCS, SYN HS"
      "Status quo SYN HS",
      "Status quo SYN QCS",
      "Restricted"
    )
  )) |>
  mutate(survey_abbrev = factor(survey_abbrev,
    levels = c(
      "SYN WCHG",
      "HBLL OUT N",
      # "SYN QCS, SYN HS"
      # "SYN QCS, SYN HS"
      "SYN HS",
      "SYN QCS"
    )
  )) |>
  mutate(lwr = lwr / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(upr = upr / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(est = est / exp(mean(log(est), na.rm = TRUE))) |>
  mutate(upr = ifelse(upr > 2 * max(est, na.rm = TRUE), 2 * max(est, na.rm = TRUE), upr)) |>
  mutate(spp_survey = gsub("Rougheye/Blackspotted Rockfish", "Rougheye/Black. Rockfish", spp_survey)) |>
  arrange(survey_abbrev, -prop_mpa, species_common_name) |>
  ggplot(aes(year, est,
    ymin = lwr, ymax = upr,
    colour = colour_var
  )) +
  geom_linerange(position = position_dodge(width = dodge_width)) +
  geom_point(position = position_dodge(width = dodge_width), pch = 21, size = 1.8) +
  geom_point(position = position_dodge(width = dodge_width), pch = 20, size = 1.8, alpha = 0.2) +
  # coord_cartesian(
  #   expand = FALSE,
  #   xlim = range(index$year) + c(-0.5, 0.5),
  #   ylim = c(0, NA)
  # ) +
  labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
  scale_colour_manual(values = cols) +
  # scale_colour_brewer(palette = "Set2") +
  facet_wrap(~ forcats::fct_inorder(stringr::str_to_title(paste(survey_abbrev, spp_survey, sep = "-"))), scales = "free_y", ncol = 5, labeller = labeller_fn) +
  ylab("Relative abundance or biomass") +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  theme(
    strip.text = element_text(colour = "black"),
    legend.position = "top", axis.text.y = element_text(size = 8)
  ) +
  # geom_smooth(method = "gam", se = F, alpha = 0.5, formula = y ~ s(x, k = 8), method.args = list(family = Gamma(link = "log")), lwd = 0.8 )
  geom_smooth(method = "loess", se = F, alpha = 0.5, formula = y ~ x, lwd = 0.8 )
  # scale_y_log10()
g

# g <- g + tagger::tag_facets(tag_prefix = "(", position = list(x = 0.1, y = 0.87), tag = "panel")

ggsave("figs/index-geo-restricted-highlights.pdf", width = 9.5, height = 5.25)
ggsave("figs/index-geo-restricted-highlights.png", width = 9.5, height = 5.25)
