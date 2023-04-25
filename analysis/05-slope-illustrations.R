library(dplyr)
library(ggplot2)

hbll <- readRDS("data-generated/index-hbll-geo-clean.rds")
syn <- readRDS("data-generated/index-syn-geo-clean.rds")
index <- bind_rows(hbll, syn)

index <- mutate(index, species_common_name = gsub(
  "rougheye/blackspotted rockfish complex",
  "rougheye/blackspotted rockfish", species_common_name
))

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

design <- readRDS("data-generated/stratified-random-design-all.rds") |>
  filter(est_type == "bootstrap") |>
  filter(species_common_name == "english sole", survey_abbrev == "SYN WCHG") |>
  filter(type %in% c("Status quo", "Restricted and shrunk")) |>
  filter(year != 2014) |>
  mutate(est = est * 1000) |>
  mutate(lwr = lwr * 1000) |>
  mutate(upr = upr * 1000)

design$lwr[design$lwr == 0] <- 0.1

d <- index |> filter(species_common_name == "english sole", survey_abbrev == "SYN WCHG") |>
  filter(type %in% c("Status quo", "Restricted and shrunk")) |>
  filter(year != 2014)

# g <- ggplot(d, aes(year, log(est), ymin = log(lwr), ymax = log(upr), colour = type)) + geom_pointrange(position = position_dodge(width = 0.5), pch = 21) +
#   # scale_y_log10() +
#   theme_light() +
#   scale_color_brewer(palette = "Dark2") +
#   geom_smooth(method = "lm", se = FALSE) +
#   # geom_smooth(method = "glm", se = FALSE, method.args = list(family = Gamma(link = "log"))) +
#   ggtitle("Model-based") +
#   theme(legend.position = "bottom")
#
# g2 <- ggplot(design, aes(year, log(est), ymin = log(lwr), ymax = log(upr), colour = type)) + geom_pointrange(position = position_dodge(width = 0.5), pch = 21) +
#   theme_light() +
#   # scale_y_log10() +
#   scale_color_brewer(palette = "Dark2") +
#   # geom_smooth(method = "glm", se = FALSE, method.args = list(family = Gamma(link = "log"))) +
#   geom_smooth(method = "lm", se = FALSE) +
#   ggtitle("Design-based") +
#   theme(legend.position = "bottom")
#
# cowplot::plot_grid(g, g2)

# d <- index |> filter(species_common_name == "redstripe rockfish", survey_abbrev == "SYN WCHG") |>
#   filter(type %in% c("Status quo", "Restricted and shrunk"))
  # filter(year != 2014)

d <- index |>
  filter(type %in% c("Status quo", "Restricted and shrunk")) |>
  filter(!(year == 2014 & survey_abbrev == "SYN WCHG"))
  # mutate(species_common_name = paste(survey_abbrev, species_common_name))

# source("analysis/theme.R")

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

make_ts_plot <- function(dat, k = 8, method = "loess", se = TRUE, center = FALSE) {
  dat$est <- dat$est * 1
  dat$lwr <- dat$lwr * 1
  dat$upr <- dat$upr * 1
  dat$log_est <- NULL

  if (center) {
    dat <- dat |>
      group_by(species_common_name, survey_abbrev, type) |>
      mutate(lwr = lwr / exp(mean(log(est), na.rm = TRUE))) |>
      mutate(upr = upr / exp(mean(log(est), na.rm = TRUE))) |>
      mutate(est = est / exp(mean(log(est), na.rm = TRUE))) |>
      ungroup()
  }

  dat <- mutate(dat,
    facet_panel = paste0(stringr::str_to_title(species_common_name), "\n(", mround(prop_mpa, 2), ")")) |>
    mutate(facet_panel = forcats::fct_reorder(facet_panel, -prop_mpa))

  g <- ggplot(dat, aes(year, log10(est), colour = type, fill = type)) +
    geom_pointrange(mapping = aes(ymin = log10(lwr), ymax = log10(upr)),
      position = position_dodge(width = 0.5), pch = 21, fill = NA) +
    # ggsidekick::theme_sleek() +
    # scale_color_brewer(palette = "Dark2") +
    scale_color_manual(values = c("Status quo" = "grey40", "Restricted and shrunk" = "red")) +
    scale_fill_manual(values = c("Status quo" = "grey40", "Restricted and shrunk" = "red")) +
    # geom_smooth(method = "lm", se = FALSE) +
    # geom_smooth(method = "glm", se = FALSE, method.args = list(family = Gamma(link = "log"))) +
    # geom_smooth(method = "gam", se = FALSE, method.args = list(family = Gamma(link = "log")), formula = y ~ s(x, bs = "cs", k = k)) +
    # geom_smooth(method = "gam", se = TRUE, formula = y ~ s(x, bs = "cs", k = k)) +
    # geom_smooth(se = TRUE, method = "loess", formula = y ~ x, alpha = 0.3) +
    geom_smooth(se = se, method = method, formula = y ~ x, alpha = 0.25) +
    ggtitle(unique(dat$survey_abbrev)) +
    # scale_y_continuous(labels = function(x) round(10^x)) +
    # scale_y_log10() +
    theme(legend.position = "top") +
    facet_wrap(~facet_panel, scales = "free_y") +
    ylab("Log10 Index") + xlab("") +
    labs(colour = "Scenario", fill = "Scenario") +
    ggplot2::theme_set(ggsidekick::theme_sleek()) +
    theme(panel.grid = element_line(colour = "grey85"),
      panel.grid.major = element_line(linewidth = rel(0.5)),
      panel.grid.minor = element_line(linewidth = rel(0.25)))
  g
}

d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.2) |>
   make_ts_plot(method = "lm", se = FALSE)
ggsave("figs/ts-lm.pdf", width = 13, height = 8)

d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.2) |>
  make_ts_plot(method = "loess", se = T)
ggsave("figs/ts-loess.pdf", width = 13, height = 8)

d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.2) |>
   make_ts_plot(method = "lm", se = FALSE, center = TRUE)
ggsave("figs/ts-lm-center.pdf", width = 13, height = 8)

d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.20) |>
  make_ts_plot(method = "loess", se = TRUE, center = TRUE)
ggsave("figs/ts-loess-center.pdf", width = 13, height = 8)


# d |> filter(survey_abbrev == "SYN WCHG") |>
#   filter(prop_mpa > 0.25) |>
#   make_ts_plot(method = "loess")
#
#
# d |> filter(survey_abbrev == "HBLL OUT N") |>
#   filter(prop_mpa > 0.25) |>
#   make_ts_plot(k = 8)

#########

surv_dat <- readRDS("data-generated/dat_to_fit.rds")
surv_dat_ll <- readRDS("data-generated/dat_to_fit_hbll.rds")

out <- surv_dat |>
  group_by(species_common_name, survey_abbrev, year) |>
  summarize(
    prop_set_mpa = sum(density_kgpm2[restricted == TRUE], na.rm = TRUE) / sum(density_kgpm2, na.rm = TRUE), prop_pos = mean(density_kgpm2 > 0)
    # prop_loc_mpa = mean(restricted)
  ) |>
  mutate(prop_pos = mean(prop_pos)) |>
  ungroup()
  # group_split() |>
  # purrr::map_dfr(function(x) {
  #   m <- lm(prop_set_mpa ~ year, data = x)
  #   data.frame(
  #     species_common_name = unique(x$species_common_name),
  #     survey_abbrev = unique(x$survey_abbrev),
  #     slope = coef(m)[[2]])
  # })


ind <- d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.2) |>
  mutate(
    facet_panel = paste0(stringr::str_to_title(species_common_name), "\n(", mround(prop_mpa, 2), ")")) |>
  mutate(facet_panel = forcats::fct_reorder(facet_panel, -prop_mpa)) |>
  select(species_common_name, survey_abbrev, facet_panel) |>
  distinct()

d |> filter(survey_abbrev == "SYN WCHG") |>
  filter(prop_mpa > 0.20) |>
  make_ts_plot(method = "lm", se = F, center = TRUE)

out |>
  filter(year != 2014) |>
  mutate(species_common_name = gsub(
    "rougheye/blackspotted rockfish complex",
    "rougheye/blackspotted rockfish", species_common_name
  )) |>
  right_join(ind) |>
  ggplot(aes(year, prop_set_mpa)) +
  geom_point() +
  facet_wrap(~facet_panel) +
  geom_smooth(method = "gam", se = F, colour = "red", formula = y ~ s(x, k = 8)) +
  ggtitle("SYN WCHG") +
  ggplot2::theme_set(ggsidekick::theme_sleek()) +
  theme(panel.grid = element_line(colour = "grey95"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25)))

ggsave("figs/prop-mpa-hidden.pdf", width = 13, height = 8)

comp_dat <- out |>
  filter(year != 2014) |>
  mutate(species_common_name = gsub(
    "rougheye/blackspotted rockfish complex",
    "rougheye/blackspotted rockfish", species_common_name
  )) |>
  right_join(ind)
comp_dat <- comp_dat |>
  group_by(species_common_name, survey_abbrev) |>
  group_split() |>
  purrr::map_dfr(function(x) {
    m <- lm(prop_set_mpa ~ year, data = x)
    data.frame(
      species_common_name = unique(x$species_common_name),
      survey_abbrev = unique(x$survey_abbrev),
      slope_prop_mpa = coef(m)[[2]] * 10,
      prop_pos = mean(x$prop_pos)
    )
  })

met <- readRDS("data-generated/metrics-wide2.rds") |>
  filter(type %in% "Restricted and shrunk") |>
  filter(est_type == "geostat") |>
  mutate(species_common_name = tolower(species_common_name)) |>
  right_join(comp_dat)

# col <- RColorBrewer::brewer.pal(5, "Blues")[5]

# ggplot(met, aes(prop_mpa, slope_re_med)) +
#   geom_point()

prep_cols_syn <- function(d) {
  d$area_swept1 <- d$doorspread_m * d$tow_length_m
  d$area_swept2 <- d$doorspread_m * d$duration_min * d$speed_mpm
  d$area_swept <- ifelse(!is.na(d$tow_length_m),
    d$area_swept1, d$area_swept2
  )
  d <- dplyr::filter(d, !is.na(area_swept))
  d
}
# prep_cols_hbll <- function(d) {
#   d$offset <- log(d$hook_count)
#   d$response <- d$catch_count
#   d
# }

surv_dat <- prep_cols_syn(surv_dat)


select(surv_dat, survey_abbrev, year, restricted, X, Y) |> distinct() |> group_by(survey_abbrev, year) |> summarise(prop = mean(restricted)) |> ggplot(aes(year, prop, colour = survey_abbrev)) + geom_point() + geom_smooth(method = "lm")

xx <- select(surv_dat, survey_abbrev, year, restricted, X, Y, area_swept) |>
  distinct() |>
  group_by(survey_abbrev, year, restricted) |>
  summarise(area_swept_total = sum(area_swept), .groups = "drop_last") |>
  summarise(prop = area_swept_total[restricted == TRUE] / sum(area_swept_total), .groups = "drop") |>
  filter(survey_abbrev == "SYN WCHG")
m <- lm(prop ~ year, data = xx)
survey_slope <- coef(m)[[2]] * 10


pal <- RColorBrewer::brewer.pal(5, "Dark2")
spp_foc <- c(
  "north pacific spiny dogfish",
  "redstripe rockfish",
  "canary rockfish",
  "petrale sole",
  "widow rockfish"
)

names(pal) <- spp_foc
other_spp <- unique(met$species_common_name)
other_spp <- other_spp[!other_spp %in% spp_foc]
pal2 <- rep("grey50", length(other_spp))
names(pal2) <- other_spp
pal_all <- c(pal, pal2)

g2 <- out |>
  filter(year != 2014) |>
  mutate(species_common_name = gsub(
    "rougheye/blackspotted rockfish complex",
    "rougheye/blackspotted rockfish", species_common_name
  )) |>
  right_join(ind) |>
  filter(species_common_name %in%
      c(
        "north pacific spiny dogfish",
        "redstripe rockfish",
        "canary rockfish",
        "petrale sole",
        "widow rockfish"
        )
    ) |>
  ggplot(aes(year, prop_set_mpa, colour = species_common_name, group = species_common_name)) +
  geom_point() +
  # facet_wrap(~facet_panel, ncol = 1, scales = "free_y") +
  facet_wrap(~stringr::str_to_title(species_common_name), ncol = 1L, scales = "fixed") +
  # geom_smooth(method = "loess", se = T) +
  # geom_smooth(method = "loess", se = F) +
  # geom_smooth(method = "betareg::betareg", se = FALSE) +
  geom_smooth(method = "gam", se = F, formula = y ~ s(x, k = 7)) +
  # geom_smooth(method = "lm", se = F, formula = y ~ x) +
  # ggtitle("SYN WCHG") +
  ggplot2::theme_set(ggsidekick::theme_sleek()) +
  theme(panel.grid = element_line(colour = "grey95"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25))) +
  # scale_colour_brewer(palette = "Dark2") +
  scale_colour_manual(values = pal_all) +
  guides(colour = "none") +
  ylab("Density proportion within MPA") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(0, 1), expand = FALSE, xlim = c(2005, 2023)) +
  xlab("Year")

g1 <- met |>
  filter(species_common_name != "english sole") |>  #FIXME!?
  ggplot(aes(slope_prop_mpa - survey_slope, slope_re_med)) +
  geom_point(aes(colour = species_common_name), pch = 19, alpha= 0.9, size = 2.5) +
  geom_point(pch = 21, col = "grey40", size = 2.5) +
  # scale_colour_viridis_c(end = 0.95, option = "F", direction = -1) +
  scale_colour_manual(values = pal_all) +
  scale_size_area() +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", alpha = 0.2) +
  geom_vline(xintercept = 0, lty = 2, col = "grey60") +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  ggrepel::geom_text_repel(aes(label = stringr::str_to_title(species_common_name)), size = 3, colour = "grey40") +
  xlab("Slope of density proportion within MPA per decade") +
  ylab("Slope of relative error per decade") +
  # labs(size = "Proportion\npositive sets", colour = "Proportion\ndensity in MPA") +
  labs(colour = "Proportion\ndensity in MPA") +
  theme(legend.position = c(0.11, 0.13), legend.box = "horizontal") +
  geom_abline(intercept = 0, slope = -1, col = "grey60", lty = 3) +
  guides(colour = "none") +
  theme(panel.grid = element_line(colour = "grey95"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25)))
  # coord_fixed()

cowplot::plot_grid(g1, g2, ncol = 2, align = "h", axis = "t", rel_widths = c(0.7, 0.3))

ggsave("figs/slopes-wchg2.pdf", width = 7, height = 5)

#
# # with fitted?
#
# glimpse(syn)
#
# slopes2 <- syn |>
#   filter(type %in% c("MPA only", "Status quo"), survey_abbrev == "SYN WCHG") |>
#   group_by(species_common_name, survey_abbrev, year) |>
#   summarise(prop_dens_mpa = sum(est[type == "MPA only"]) / sum(est[type == "Status quo"]), .groups = "drop_last") |>
#   group_split() |>
#   purrr::map_dfr(function(x) {
#     m <- lm(prop_dens_mpa ~ year, data = x)
#     data.frame(
#       species_common_name = unique(x$species_common_name),
#       survey_abbrev = unique(x$survey_abbrev),
#       slope_prop_mpa = coef(m)[[2]] * 10
#     )
#   })
#
# met_geo <- readRDS("data-generated/metrics-wide2.rds") |>
#   filter(type %in% "Restricted and shrunk") |>
#   filter(est_type == "geostat") |>
#   mutate(species_common_name = tolower(species_common_name)) |>
#   right_join(slopes2)
#
# met_geo |>
#   # filter(species_common_name != "english sole") |>  #FIXME!?
#   ggplot(aes(slope_prop_mpa, slope_re_med)) +
#   geom_point(aes(colour = prop_mpa), pch = 19, alpha= 0.9, size = 2.5) +
#   geom_point(aes(colour = prop_mpa), pch = 21, col = "grey40", size = 2.5) +
#   scale_colour_viridis_c(end = 0.95, option = "F", direction = -1) +
#   # scale_color_distiller(palette = "YlOrRd", direction = -1) +
#   scale_size_area() +
#   geom_smooth(method = "lm", se = TRUE, colour = "grey30", alpha = 0.2) +
#   geom_vline(xintercept = 0, lty = 2, col = "grey60") +
#   geom_hline(yintercept = 0, lty = 2, col = "grey60") +
#   ggrepel::geom_text_repel(aes(label = stringr::str_to_title(species_common_name)), size = 3, colour = "grey40") +
#   xlab("Slope of density proportion within MPA per decade") +
#   ylab("Slope of relative error per decade") +
#   # labs(size = "Proportion\npositive sets", colour = "Proportion\ndensity in MPA") +
#   labs(colour = "Proportion\ndensity in MPA") +
#   theme(legend.position = c(0.11, 0.13), legend.box = "horizontal") +
#   geom_abline(intercept = 0, slope = -1, col = "grey60", lty = 3)
#
