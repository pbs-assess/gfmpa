# make plots that combine surveys

library(tidyverse)
theme_set(ggsidekick::theme_sleek())

include_mpa <- TRUE
# include_mpa <- FALSE

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
index$type_label <- factor(index$type_label, levels = c("Extrapolated", "Shrunk", "Status quo", "MPA only"))

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
    # "Longnose Skate",#HS
    # "Spotted Ratfish",

    "Pacific Cod",
    "Walleye Pollock",#
    # "Pacific Hake",
    "Sablefish",
    # "Lingcod",
    # "Blackbelly Eelpout",

    "Canary Rockfish",
    # "Copper Rockfish",
    "Greenstriped Rockfish",#HS
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
    "Curlfin Sole"# QCS
    # "Dover Sole",
    # "English Sole",
    # "Flathead Sole",
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

  if (include_mpa) colour_pal <- c("#FC8D62", "#8DA0CB", "gray60", "gray80")
  if (!include_mpa) colour_pal <- c("#FC8D62", "#8DA0CB", "gray")
  if (include_mpa) line_pal <- c("solid", "solid", "dotted", "solid")
  if (!include_mpa) line_pal <- c("solid", "solid", "dotdash")

  g1 <- i %>%
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
    theme(legend.justification=c(0,1), legend.position=c(-0.28,1.065), legend.direction = "horizontal",
      strip.text.y = element_text(size = 7, angle = 0, hjust = 0),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0.2))

  g1 + g2 + patchwork::plot_layout(widths = c(1,2))

if (include_mpa) ggsave("figs/index-geo-restricted-highlights.pdf", width = 7, height = 10)

if (!include_mpa) ggsave("figs/index-geo-restricted-highlights-noMPA.pdf", width = 6.5, height = 10)

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

  library(scales)
  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

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

  ggsave("figs/index-geo-restricted-re-highlights.pdf", width = 5.75, height = 9.5)
