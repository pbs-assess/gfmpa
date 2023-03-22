library(dplyr)
library(ggplot2)
dall <- readRDS("data-generated-ALL/metrics-wide.rds")
d12 <- readRDS("data-generated/metrics-wide.rds")
source("analysis/theme.R")

dall <- dall |> select(survey_abbrev, species_common_name, restr_clean,
  prop_mpa_all = prop_mpa, mare_all = mare)
d12 <- d12 |> select(survey_abbrev, species_common_name, restr_clean,
  prop_mpa_cat12 = prop_mpa, mare_cat12 = mare)

d <- left_join(dall, d12, by = join_by(survey_abbrev, species_common_name, restr_clean))

d <- mutate(d, mare_inc = mare_all / mare_cat12)

d <- d |>
  mutate(survey_abbrev = factor(survey_abbrev,
    levels = c(
      "SYN WCHG",
      "HBLL OUT N",
      "SYN QCS, SYN HS"))
  )

ggplot(d, aes(prop_mpa_cat12, prop_mpa_all, colour = survey_abbrev, size = mare_inc)) +
  geom_abline(intercept = 0, slope = 1, lty = 2, colour = "grey50") +
  geom_point(pch = 21) +
  geom_point(pch = 19, alpha = 0.3) +
  scale_colour_manual(values = restricted_cols) +
  xlab("Proportion species in MPA\n(NSB category 1 + 2 only)") +
  ylab("Proportion species in MPA\n(NSB category 1 + 2 + 'as-is where-is')") +
  labs(colour = "Survey", size = "Fold increase in MARE") +
  theme(legend.position = "right", legend.direction = "vertical") +
  coord_equal() +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.02))) +
  scale_size_area(max_size = 7, breaks = seq(1, 7, 2))
  # theme(legend.position = c(0.8, 0.2), legend.box.background = element_rect(colour = "grey60"))
ggsave("figs/asis-whereis-1-1.pdf", width = 7, height = 5)


