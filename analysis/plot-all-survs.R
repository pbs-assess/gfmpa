# make plots that combine surveys

library(tidyverse)
theme_set(ggsidekick::theme_sleek())

include_mpa <- TRUE
include_mpa <- FALSE

y1 <- readRDS(file = "data-generated/index-hbll-geo-clean.rds") %>%
  mutate(est = est/10000, lwr = lwr/10000, upr = upr/10000)
y2 <- readRDS(file = "data-generated/index-syn-geo-clean.rds")
y <- bind_rows(y1,y2)

mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
filter(y, orig_cv <= 1)
index <- filter(y, orig_cv < 1)

index$type_label <- index$type
index[index$type == "Restricted",]$type_label <- "Extrapolated"
index[index$type == "Restricted and shrunk",]$type_label <- "Shrunk"
index$type_label <- factor(index$type_label, levels = c("Status quo", "Extrapolated", "Shrunk", "MPA only"))

index$species_common_name <- stringr::str_to_title(index$species_common_name)

sp_list <- index %>% select(survey_abbrev, species_common_name) %>%
    distinct() %>% group_by(species_common_name) %>%
    summarise(n = n())
sp_list2 <- sp_list %>% filter(n > 1)

# syn_highlights <- c(sp_list2$species_common_name)

syn_highlights <- c(
    "North Pacific Spiny Dogfish",
    # "Big Skate",
    "Sandpaper Skate",#
    "Longnose Skate",#HS
    # "Spotted Ratfish",

    "Pacific Cod",
    "Walleye Pollock",#
    # "Pacific Hake",
    "Sablefish",
    # "Lingcod",
    # "Blackbelly Eelpout",

    # "Canary Rockfish",
    # "Copper Rockfish",
    # "Greenstriped Rockfish",#HS
    "Quillback Rockfish",
    "Redbanded Rockfish",#HS
    # # "Redstripe Rockfish",
    # "Rosethorn Rockfish",
    # "Rougheye/Blackspotted Rockfish Complex",#HBLL
    # "Silvergray Rockfish",
    "Yelloweye Rockfish",
    # "Yellowtail Rockfish",
    # "Shortspine Thornyhead",

    "Arrowtooth Flounder",#
    # "Curlfin Sole",# QCS
    # "Dover Sole",
    # "English Sole",
    "Flathead Sole"
    # "Pacific Halibut",
    # "Rex Sole",
    # "Slender Sole",
    # "Southern Rock Sole"
  )


if (!include_mpa) index <- index %>% filter(type != "MPA only")


  i1 <- index %>% filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    filter((survey_abbrev == "HBLL OUT N"))
  i2 <- index %>% filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    filter(!(survey_abbrev == "HBLL OUT N"))

  spp <- as.data.frame(syn_highlights) %>% rename(species_common_name=syn_highlights)
  i <- left_join(spp, i1) %>% mutate(survey_abbrev = "HBLL OUT N")

  # get two middle colours from plot above
  library(RColorBrewer)
  brewer.pal(4, "Set2")

  if (include_mpa) colour_pal <- c("gray60", "#FC8D62", "#8DA0CB", "gray80")
  if (!include_mpa) colour_pal <- c("gray60", "#FC8D62", "#8DA0CB")
  if (include_mpa) line_pal <- c("dotted", "solid", "solid", "solid")
  if (!include_mpa) line_pal <- c("dotted", "solid", "solid")

  g1 <- i %>%
    filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr,
               colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.6) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Type", fill = "Type", linetype = "Type") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    xlim(2005, 2020) +
    facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
    scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    ggtitle("Index type:   ") +
    ylab("Relative abundance in 1000s (HBLL) or biomass in tonnes (SYN)") +
    theme(plot.title = element_text(hjust = 0.9),
          strip.text.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")

  g2 <- i2 %>% filter(survey_abbrev != "SYN WCHG") %>%
    filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr,
               colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.6) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
    scale_y_continuous(breaks = waiver(), n.breaks = 3) +
    ylab("Relative biomass") +
    ggtitle("") +
    theme(legend.justification=c(0,1), legend.position=c(-0.28,1.085), legend.direction = "horizontal",
      strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0.2))

  g1 + g2 + patchwork::plot_layout(widths = c(1,2))

if (include_mpa) ggsave("figs/index-geo-restricted-highlights.pdf", width = 6.5, height = 8)
if (!include_mpa) ggsave("figs/index-geo-restricted-highlights-noMPA.pdf", width = 6.5, height = 8)

  x <- index %>%
    group_by(species_common_name, survey_abbrev, type) %>%
    mutate(est = est / exp(mean(log(est)))) %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
      re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
    )

  lu <- tibble(
    "Restriction type" = c("re_restr", "re_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )

  x_long <- x %>%
    tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") %>%
    left_join(lu)

  # library(scales)
  # S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  # IS_sqrt <- function(x){x^2*sign(x)}
  # S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

  g <- x_long %>% filter(survey_abbrev != "SYN WCHG") %>%
    filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    ylab("Relative error") + xlab("Year") +
    facet_grid(species_common_name~survey_abbrev, scales = "free") +
    scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    # scale_y_continuous(trans = "S_sqrt", breaks = c(-0.5,-0.1,0, 0.1, 0.5)) +
    # coord_cartesian(ylim = c(-0.35, 0.4)) +
    labs(colour = "") +
    ggtitle("") +
    theme(legend.justification=c(0.5,1), legend.position=c(0.5,1.065), legend.direction = "horizontal",
          # legend.position = "top", legend.justification=c(0.5,1),
          strip.text.y = element_text(size = 7, angle = 0, hjust = 0))
  g

  ggsave("figs/index-geo-restricted-re-highlights.pdf", width = 5.75, height = 8)


# supplements
  g <- index %>% filter(survey_abbrev == "HBLL OUT N") %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    ylab("Relative abundance in 1000s") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5)
  g
  ggsave("figs/index-hbll-geo-restricted.pdf", width = 12, height = 8)

  g <- index %>% filter(survey_abbrev == "SYN QCS") %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    ylab("Relative biomass in tonnes") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5)
  g
  ggsave("figs/index-qcs-geo-restricted.pdf", width = 12, height = 13, limitsize = FALSE)

  g <- index %>% filter(survey_abbrev == "SYN HS") %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    ylab("Relative biomass in tonnes") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5)
  g
  ggsave("figs/index-hs-geo-restricted.pdf", width = 12, height = 13, limitsize = FALSE)

  g <- index %>% filter(survey_abbrev == "SYN WCHG") %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    ylab("Relative biomass in tonnes") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5)
  g
  ggsave("figs/index-wchg-geo-restricted.pdf", width = 12, height = 5, limitsize = FALSE)

  #relative error

  g <- x_long %>% filter(survey_abbrev == "HBLL OUT N") %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    ylab("Relative error") + xlab("Year") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
    # scale_y_continuous(trans = "S_sqrt", breaks = c(-0.5,-0.1,0, 0.1, 0.5)) +
    # coord_cartesian(ylim = c(-0.35, 0.4)) +
    labs(colour = "")
  g
  ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 12, height = 8)

  g <- x_long %>% filter(survey_abbrev == "SYN QCS") %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    ylab("Relative error") + xlab("Year") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
    labs(colour = "")
  g
  ggsave("figs/index-qcs-geo-restricted-re.pdf", width = 12, height = 13)

  g <- x_long %>% filter(survey_abbrev == "SYN HS") %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    ylab("Relative error") + xlab("Year") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
    labs(colour = "")
  g
  ggsave("figs/index-hs-geo-restricted-re.pdf", width = 12, height = 13)

  g <- x_long %>% filter(survey_abbrev == "SYN WCHG") %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    ylab("Relative error") + xlab("Year") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
    labs(colour = "")
  g
  ggsave("figs/index-wchg-geo-restricted-re.pdf", width = 12, height = 5)

  g <- x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    left_join(lu) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est), est, colour = restr_clean)) +
    geom_point(size = 2, position = position_dodge(width = 0)) +
    xlab("") +
    ylab("Median absolute relative error (MARE) of restricted index compared to status quo") +
    labs(colour = "") +
    coord_flip() + scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90")) +
    # guides(colour = guide_legend(nrow = 2L))+
    facet_wrap(~survey_abbrev, ncol = 4)
  g
  ggsave("figs/index-geo-mare-dotplot.pdf", width = 9, height = 8)


  x <- index %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      cv_ratio_restr = cv[type == "Restricted"] /
        cv[type == "Status quo"],
      cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
        cv[type == "Status quo"]
    )
  x_long <- x %>%
    tidyr::pivot_longer(starts_with("cv"),
                        names_to = "Restriction type", values_to = "CV ratio")
  x_long %>%
    group_by(`Restriction type`, survey_abbrev) %>%
    summarise(mean_ratio = mean(`CV ratio`)) %>%
    knitr::kable(digits = 2)
  lu <- tibble(
    "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )
  g <- x_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
    left_join(lu) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est), est, colour = restr_clean, ymin = lwr, ymax = upr)) +
    geom_hline(yintercept = 1, lty = 2, col = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.5)) +
    coord_flip() +
    xlab("") +
    ylab("Ratio of index CV (restricted/status quo)") +
    labs(colour = " ") +
    scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk")) +
    theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"))+
    facet_wrap(~survey_abbrev, ncol = 4, scales = "free_x")
  g
  ggsave("figs/index-geo-cv-ratio-dotplot.pdf", width = 9, height = 9)
