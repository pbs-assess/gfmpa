library(dplyr)
library(ggplot2)
source("analysis/theme.R")

index <- readRDS("data-generated/index-filtered.rds")
metrics <- readRDS("data-generated/metrics-wide.rds")

top_prop <- metrics |>
  select(
    prop_mpa, species_common_name,
    survey_abbrev
  ) |>
  distinct() |>
  group_by(survey_abbrev) |>
  top_n(n = 8, wt = prop_mpa) |>
  arrange(survey_abbrev, -prop_mpa) |>
  mutate(prop_mpa = round(prop_mpa, 2)) |>
  filter(prop_mpa > 0.1)

top_prop |> as.data.frame()

hbll_fig2 <-
  c(
    "shortspine thornyhead",
    "china rockfish",
    "yelloweye rockfish",
    "big skate",
    "copper rockfish"
  )

syn_qcshs_fig2 <-
  c(
    "yellowmouth rockfish",
    "redbanded rockfish"
    # "shortspine thornyhead"
  )


syn_wchg_fig2 <-
  c(
    # "canary rockfish",
    "walleye pollock",
    "rougheye/blackspotted rockfish",
    "petrale sole",
    "bocaccio",
    "north pacific spiny dogfish",
    "pacific ocean perch",
    "redstripe rockfish",
    "silvergray rockfish"
  )

fig2_keep <-
  bind_rows(
    tibble(species_common_name_lower = hbll_fig2, survey_abbrev = "HBLL OUT N"),
    tibble(species_common_name_lower = syn_qcshs_fig2, survey_abbrev = "SYN QCS, SYN HS"),
    tibble(species_common_name_lower = syn_wchg_fig2, survey_abbrev = "SYN WCHG")
  )

dodge_width <- 1

colour_pal <- c("gray50", restricted_cols)

m <- select(metrics, survey_abbrev, species_common_name, cv_orig, prop_mpa) |>
  distinct()

index <- index |> left_join(m, by = join_by(species_common_name, survey_abbrev))

g <- index |>
  mutate(species_common_name_lower = tolower(species_common_name)) |>
  # semi_join(top_prop, by = join_by(survey_abbrev, species_common_name_lower)) |>
  semi_join(fig2_keep, by = join_by(survey_abbrev, species_common_name_lower)) |>
  mutate(spp_survey = species_common_name) |>
  mutate(spp_survey = gsub("HBLL OUT N", "HBLL", spp_survey)) |>
  filter(type %in% c("Status quo", "Restricted and shrunk")) |>
  group_by(spp_survey, type) |>
  mutate(colour_var = ifelse(type == "Status quo", "Status quo",
    paste0("Restricted ", survey_abbrev)
  )) |>
  mutate(colour_var = factor(colour_var,
    levels = c(
      "Status quo",
      "Restricted SYN WCHG",
      "Restricted HBLL OUT N",
      "Restricted SYN QCS, SYN HS"
    )
  )) |>
  mutate(survey_abbrev = factor(survey_abbrev,
    levels = c(
      "SYN WCHG",
      "HBLL OUT N",
      "SYN QCS, SYN HS"
    )
  )) |>
  mutate(lwr = lwr / exp(mean(log(est)))) |>
  mutate(upr = upr / exp(mean(log(est)))) |>
  mutate(est = est / exp(mean(log(est)))) |>
  mutate(upr = ifelse(upr > 2 * max(est), 2 * max(est), upr)) |>
  mutate(spp_survey = gsub("Rougheye/Blackspotted Rockfish", "Rougheye/Black. Rockfish", spp_survey)) |>
  arrange(survey_abbrev, -prop_mpa, species_common_name) |>
  ggplot(aes(year, est,
    ymin = lwr, ymax = upr,
    colour = colour_var
  )) +
  geom_linerange(position = position_dodge(width = dodge_width)) +
  geom_point(position = position_dodge(width = dodge_width), pch = 21, size = 1.8) +
  geom_point(position = position_dodge(width = dodge_width), pch = 20, size = 1.8, alpha = 0.2) +
  coord_cartesian(
    expand = FALSE,
    xlim = range(index$year) + c(-0.5, 0.5),
    ylim = c(0, NA)
  ) +
  labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
  scale_colour_manual(values = colour_pal) +
  facet_wrap(~ forcats::fct_inorder(spp_survey), scales = "free_y", ncol = 5) +
  ylab("Relative abundance or biomass") +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  theme(
    strip.text = element_text(colour = "black"),
    legend.position = "top", axis.text.y = element_text(size = 8)
  )
g

g <- g + tagger::tag_facets(tag_prefix = "(", position = list(x = 0.1, y = 0.87))

ggsave("figs/index-geo-restricted-highlights.pdf", width = 9.5, height = 5.25)
