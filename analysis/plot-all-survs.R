# make plots that combine results from all species and surveys
# need to run geostat-all.R and then rel-error-slopes.R
# each with globals set for each survey separately
# currently using the "binomial-gamma" model type for both surveys
# for this script the only global is whether to include the mpa index in figure 1

library(tidyverse)
library(egg)
theme_set(ggsidekick::theme_sleek())

# Globals to set ------------------------------
# include_mpa <- TRUE
include_mpa <- FALSE
# ---------------------------------------------


# set plot colours and linetypes ----

# choose two colours for restricted indices that are colourblind friendly
# library(RColorBrewer)
# brewer.pal(4, "Set2")
restricted_cols <- c("#FC8D62", "#8DA0CB")
if (include_mpa) colour_pal <- c("gray50", restricted_cols, "gray80")
if (!include_mpa) colour_pal <- c("gray50", restricted_cols)
if (include_mpa) line_pal <- c("dotted", "solid", "solid", "solid")
if (!include_mpa) line_pal <- c("dotted", "solid", "solid")


# prep index data ----

y1 <- readRDS(file = "data-generated/index-hbll-geo-clean-nbinom2.rds") %>%
  # (x/4km left in predict function from density version)*199 sets to sample each grid cell
  # then /1000 to change counts to per 1000 fish
  # work out to  x 0.2
  mutate(est = est*0.2, lwr = lwr*0.2, upr = upr*0.2)
# y1 <- readRDS(file = "data-generated/index-hbll-geo-clean-binomial-gamma.rds") %>%
#   mutate(est = est/10000, lwr = lwr/10000, upr = upr/10000)
y2 <- readRDS(file = "data-generated/index-syn-geo-clean.rds")
y <- bind_rows(y1,y2)

mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
filter(y, orig_cv <= 1)
index <- filter(y, orig_cv < 1)

index$species_common_name <- stringr::str_to_title(index$species_common_name)
if (!include_mpa) index <- index %>% filter(type != "MPA only")


# choose labels for each index type ----
restricted_labels <- c("Extrapolated", "Shrunk")
index$type_label <- index$type
index[index$type == "Restricted",]$type_label <- "Extrapolated"
index[index$type == "Restricted and shrunk",]$type_label <- "Shrunk"
index$type_label <- factor(index$type_label, levels = c("Status quo", "Extrapolated", "Shrunk", "MPA only"))


# select species to highlight in figures 1 and 2 ----

# I first filtered to get species with indices from at least two surveys
# sp_list <- index %>% select(survey_abbrev, species_common_name) %>%
#     distinct() %>% group_by(species_common_name) %>%
#     summarise(n = n())
# sp_list2 <- sp_list %>% filter(n > 1)
# syn_highlights <- c(sp_list2$species_common_name)

# then I manual selected ones with interesting patterns in the scatterplots
syn_highlights <- c(
    "North Pacific Spiny Dogfish",
    # "Big Skate",
    "Sandpaper Skate",#
    # "Longnose Skate",#HS
    # "Spotted Ratfish",

    # "Pacific Cod",
    "Walleye Pollock",#
    # "Pacific Hake",
    "Sablefish",
    # "Lingcod",
    # "Blackbelly Eelpout",

    # "Canary Rockfish",
    # "Copper Rockfish",
    # "Greenstriped Rockfish",#HS
    # "Quillback Rockfish",
    "Redbanded Rockfish",#HS
    # # "Redstripe Rockfish",
    # "Rosethorn Rockfish",
    # "Rougheye/Blackspotted Rockfish Complex",#HBLL
    # "Silvergray Rockfish", # an interesting one
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

# FIGURE 1: indices ----

# split out HBLL to give it a different Y axis scale than the SYN surveys
  i1 <- index %>% filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    filter((survey_abbrev == "HBLL OUT N"))
  i2 <- index %>% filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    filter(!(survey_abbrev == "HBLL OUT N"))

  spp <- as.data.frame(syn_highlights) %>% rename(species_common_name=syn_highlights)
  i <- left_join(spp, i1) %>% mutate(survey_abbrev = "HBLL OUT N")

  g1 <- i %>%
    filter(species_common_name %in% syn_highlights) %>%
    mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr,
               colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.6) +
    geom_ribbon(alpha = 0.1, colour = NA) +
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
    geom_ribbon(alpha = 0.1, colour = NA) +
    labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
    scale_colour_manual(values = colour_pal) +
    scale_fill_manual(values = colour_pal) +
    scale_linetype_manual(values = line_pal) +
    facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
    scale_y_continuous(breaks = waiver(), n.breaks = 3) +
    ylab("Relative biomass") +
    ggtitle("") +
    theme(legend.justification=c(0,1), legend.position=c(-0.28,1.085), legend.direction = "horizontal",
      strip.text.y = element_text(size = 9, angle = 0, hjust = 0),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0.2))

  g1 + g2 + patchwork::plot_layout(widths = c(1,2))

if (include_mpa) ggsave("figs/index-geo-restricted-highlights.pdf", width = 6.5, height = 8)
if (!include_mpa) ggsave("figs/index-geo-restricted-highlights-noMPA.pdf", width = 7, height = 8)


# FIGURE 2: RE through time ----

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
    geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
    geom_line(size = 0.84, alpha = 0.9) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = restricted_cols, label= ) +
    ylab("Relative error") + xlab("Year") +
    # scale_y_continuous(trans = "S_sqrt", breaks = c(-0.5,-0.1,0, 0.1, 0.5)) +
    # coord_cartesian() + #ylim = c(-0.35, 0.4)
    labs(colour = "") +
    ggtitle("") +
    theme(legend.justification=c(0.5,1), legend.position=c(0.5,1.065), legend.direction = "horizontal",
          # legend.position = "top", legend.justification=c(0.5,1),
          strip.text.y = element_text(size = 9, angle = 0, hjust = 0))

  # with free y-axis scales
  (g <- g +
      scale_y_continuous(breaks = waiver(), n.breaks = 5) +
      facet_grid(species_common_name~survey_abbrev, scales = "free"))
  ggsave("figs/index-geo-restricted-re-highlights.pdf", width = 6.75, height = 8)

  # with fixed y-axis scales
  (g <- g +
    scale_y_continuous(breaks = waiver(), n.breaks = 3) +
    facet_grid(species_common_name~survey_abbrev, scales = "free_x"))
  ggsave("figs/index-geo-restricted-re-highlights-fixed.pdf", width = 6.75, height = 8)



# combine precision, accuracy and bias data from all surveys ----
  cvdata1 <- readRDS("data-generated/hbll-cv-w-lm-slopes-nb.rds")
  cvdata2 <- readRDS("data-generated/syn-cv-w-lm-slopes.rds")

  # glimpse(cvdata1)
  # glimpse(cvdata2)
  cvdata <- bind_rows(cvdata1, cvdata2)

  d <- cvdata %>%
    tidyr::pivot_longer(c("cv_ratio", "mare"), names_to = "Response", values_to = "cv_index") %>%
    ungroup() %>%
    filter(cv_index < 1.6) %>%
    mutate(Response = factor(Response, labels = c("CV Ratio", "MARE")))


# function for scatterplots ----
  plot_scatter <- function(dat, x, y) {
    ggplot(dat, aes_string(x, y,
                           colour = "restr_clean",
                           group = "species_common_name")) +
      geom_line(colour = "gray95") +
      ggrepel::geom_text_repel(data = filter(dat, `Restriction type` == "re_restr"),
                               aes(label = species_common_name),
                               colour = "darkgray",
                               force = 2, direction = "y", max.overlaps = 2,
                               min.segment.length = 10, size = 2) +
      geom_point() +
      theme(legend.position = "none", legend.title = element_blank()) +
      scale_colour_manual(values = restricted_cols, label=restricted_labels)
  }


# FIGURE 3: HBLL - CV ratio and MARE by status quo CV and prop MPA ----

  d2 <- filter(d, survey_abbrev == "HBLL OUT N")

  (g1 <- plot_scatter(d2, "cv_orig", "cv_index") +
      xlab("CV of 'Status quo' index") +
      guides(shape = "none") +
      facet_grid(rows=vars(Response),
                 # cols = vars(survey_abbrev),
                 switch = "y",
                 scales = "free_y") +
      theme(axis.title.y = element_blank(),
            legend.position = c(0.2,0.93),
            strip.placement = "outside"))

  (g2 <- plot_scatter(d2, "prop_mpa", "cv_index") +
      xlab("Biomass proportion inside MPAs") +
      guides(shape = "none") +
      facet_wrap(~Response, strip.position = "top", nrow = 2, scales = "free_y") +
      theme(strip.placement = "outside", strip.text.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()))

  g1 <- tag_facet(g1, tag_pool = c("a","c"), hjust = -0.5, vjust = 2, fontface = 1)
  g2 <- tag_facet(g2, tag_pool = c("b","d"), hjust = -0.5, vjust = 2, fontface = 1)

  g1 + g2 + patchwork::plot_layout(nrow = 1)

  ggsave("figs/explore-hbll-cv.pdf", width = 7, height = 7)

# FIGURE 4: QCS - CV ratio and MARE by status quo CV and prop MPA ----
  d3 <- filter(d, survey_abbrev == "SYN QCS")

  (g1 <- plot_scatter(d3, "cv_orig", "cv_index") +
      xlab("CV of 'Status quo' index") +
      guides(shape = "none") +
      facet_grid(rows=vars(Response),
                 # cols = vars(survey_abbrev),
                 switch = "y",
                 scales = "free_y") +
      theme(axis.title.y = element_blank(),
            # legend.position = c(0.2,0.95),
            strip.placement = "outside"))

  (g2 <- plot_scatter(d3, "prop_mpa", "cv_index") +
      xlab("Biomass proportion inside MPAs") +
      guides(shape = "none") +
      facet_wrap(~Response, strip.position = "top", nrow = 2, scales = "free_y") +
      theme(strip.placement = "outside", strip.text.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"),
            legend.position = c(0.25,0.97),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()))

  g1 <- tag_facet(g1, tag_pool = c("a","c"), hjust = -0.5, vjust = 2, fontface = 1)
  g2 <- tag_facet(g2, tag_pool = c("b","d"), hjust = -0.5, vjust = 2, fontface = 1)

  g1 + g2 + patchwork::plot_layout(nrow = 1)

  ggsave("figs/explore-qcs-cv.pdf", width = 7, height = 7)


# CV ratio and MARE by status quo CV for all surveys at once for appendix ----
    if(length(unique(d$survey_abbrev))> 3){ # makes sure all surveys

      (g0 <- plot_scatter(d, "cv_orig", "cv_index") +
         xlab("CV of 'Status quo' index") +
         guides(shape = "none") +
         facet_grid(rows=vars(Response),
                    cols = vars(survey_abbrev),
                    switch = "y",
                    scales = "free_y") +
         theme(axis.title.y = element_blank(),
               legend.position = c(0.1,0.95),
               strip.placement = "outside"))
      (g0 <- tag_facet(g0, fontface = 1))

      ggsave("figs/explore-all-cv-by-cv.pdf", width = 10, height = 5.5)


# CV ratio and MARE by prop MPA for all surveys at once for appendix ----
    (g <- plot_scatter(d, "prop_mpa", "cv_index") +
       xlab("Biomass proportion inside MPAs") +
       guides(shape = "none") +
       facet_grid(rows=vars(Response),
                  cols = vars(survey_abbrev),
                  switch = "y",
                  scales = "free_y") +
       theme(axis.title.y = element_blank(),
             legend.position = c(0.1,0.95),
             strip.placement = "outside"))
    (g <- tag_facet(g, fontface = 1))

    ggsave("figs/explore-all-cv-by-mpa.pdf", width = 10, height = 5.5)
  }


# FIGURE 5: MARE by CV ratio ----

  cvdata2 <- cvdata
  cvdata2[cvdata2$cv_ratio > 1.5 , ]$cv_ratio <- 1.5

  (g <- plot_scatter(cvdata2, "cv_ratio", "mare") +
      xlab("CV Ratio") +
      ylab("MARE") +
      guides(shape = "none") +
      facet_wrap(~survey_abbrev,
                 ncol = 2,
                 # cols = vars(survey_abbrev),
                 # switch = "y",
                 scales = "free_x") +
      theme(legend.position = c(0.12,0.96),
            strip.placement = "outside"))

  g <- tag_facet(g, hjust = -0.5, vjust = 2, fontface = 1)
  # g
  ggsave("figs/explore-all-mare-by-cv-ratio.pdf", width = 7, height = 7)


# FIGURE 6: slopes ----

  d_mare <- filter(d, Response == "MARE")

  (g <- plot_scatter(d_mare, "slope_mpa", "slope_re/100") +
      facet_wrap(~survey_abbrev,
                 scales = "free") +
      ylab("Change in RE per decade") +
      xlab("Change in proportion of biomass inside MPAs") +
      geom_hline(yintercept = 0, colour = "gray80") +
      geom_vline(xintercept = 0, colour = "gray70") +
      theme(legend.position = c(0.37,0.95)))
  (g <- tag_facet(g, hjust = -0.5, vjust = 2, fontface = 1))

  ggsave("figs/explore-all-slopes.pdf", width = 7, height = 7)

  (g <- plot_scatter(d_mare, "prop_mpa", "abs(slope_re/100)") +
      ylab("Absolute change in RE per decade") +
      xlab("Proportion of biomass inside MPAs") +
      facet_wrap(~survey_abbrev,
                 scales = "free") + theme(legend.position = c(0.08,0.96)))
  ggsave("figs/explore-abs-slope.pdf", width = 7, height = 7)


  # ### benefits of interpolation
  #
  # d <- cvdata %>% group_by(species_common_name, survey_abbrev) %>% mutate(
  #   temp_bias_interp_ratio = abs(slope_re[restr_clean == "Same survey domain"])/abs(slope_re[restr_clean == "Shrunk survey domain"]),
  #   cv_interp_ratio = (cv_mean[restr_clean == "Same survey domain"])/(cv_mean[restr_clean == "Shrunk survey domain"]),
  #   mare_interp_ratio = (mare[restr_clean == "Same survey domain"])/(mare[restr_clean == "Shrunk survey domain"])
  # )
  #
  # median(d$temp_bias_interp_ratio)
  # median(d$cv_interp_ratio)
  # median(d$mare_interp_ratio)


# SUPPLEMENTAL FIGURES ----

# INDEX PLOTS for each survey ----
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
    # CI on exrapolated Yellowmouth index blows up,
    # truncating here to see uncertainty on other indices
    mutate(upr = ifelse(
      species_common_name == "Yellowmouth Rockfish"&upr>10000, 10000, upr
      )) %>%
    ggplot(aes(year, est, ymin = lwr, ymax = upr,
               colour = type_label, fill = type_label, linetype = type_label)) +
    geom_line(lwd = 0.9) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    coord_cartesian() +
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


# RE PLOTS for each survey ----
  g <- x_long %>% filter(survey_abbrev == "HBLL OUT N") %>%
    ggplot(aes(year, re, colour = restr_clean)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_line(size = 0.8, alpha =0.8) +
    # scale_color_brewer(palette = "Set2") +
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
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
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
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
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
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
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
    ylab("Relative error") + xlab("Year") +
    facet_wrap(~species_common_name, scales = "free_y", ncol = 5) +
    labs(colour = "")
  g
  ggsave("figs/index-wchg-geo-restricted-re.pdf", width = 12, height = 5)


# CV ratio dotplot ----
  cv <- index %>%
    group_by(species_common_name, survey_abbrev, year) %>%
    summarise(
      cv_ratio_restr = cv[type == "Restricted"] /
        cv[type == "Status quo"],
      cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
        cv[type == "Status quo"]
    )
  cv_long <- cv %>%
    tidyr::pivot_longer(starts_with("cv"),
                        names_to = "Restriction type", values_to = "CV ratio")
  cv_long %>%
    group_by(`Restriction type`, survey_abbrev) %>%
    summarise(mean_ratio = mean(`CV ratio`)) %>%
    knitr::kable(digits = 2)
  lu_cv <- tibble(
    "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )
  g <- cv_long %>%
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
    left_join(lu_cv) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est), est,
               colour = restr_clean, ymin = lwr, ymax = upr)) +
    geom_hline(yintercept = 1, lty = 2, col = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.75), size = 0.35) +
    coord_flip(ylim = c(0.9, 1.6)) +
    xlab("") +
    ylab("Ratio of index CV (restricted/status quo)") +
    labs(colour = " ") +
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
    theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"))+
    facet_wrap(~survey_abbrev, ncol = 4, scales = "free_x")
  g
  # (g <- tag_facet_outside(g, fontface = 1))

  ggsave("figs/index-geo-cv-ratio-dotplot.pdf", width = 9, height = 9)


# MARE dotplot ----
  g <- x_long %>% # first created for FIG 2
    group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
    summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
    left_join(lu) %>%
    ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est), est,
               colour = restr_clean)) +
    geom_point(size = 2, position = position_dodge(width = 0)) +
    xlab("") +
    ylab("Median absolute relative error (MARE) of restricted index compared to status quo") +
    labs(colour = "") +
    coord_flip() + scale_y_continuous(breaks = waiver(), n.breaks = 4) +
    scale_colour_manual(values = restricted_cols, label=restricted_labels) +
    theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90")) +
    # guides(colour = guide_legend(nrow = 2L))+
    facet_wrap(~survey_abbrev, ncol = 4)
  g
  ggsave("figs/index-geo-mare-dotplot.pdf", width = 9, height = 8)

# RAW data checks ----
  # dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
  # ggplot(dat_to_fit, aes(year, hook_count)) + geom_jitter(alpha = 0.1)
  # ggplot(dat_to_fit, aes(hook_count, area_km2)) + geom_jitter(alpha = 0.1) + facet_wrap(~year)
