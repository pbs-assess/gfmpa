# make plots that combine results from all species and surveys
# need to run geostat-all.R and then rel-error-slopes.R
# each with globals set for each survey separately
# currently using the "binomial-gamma" model type for both surveys
# for this script the only global is whether to include the mpa index in figure 1
# devtools::install_github("eliocamp/tagger")

library(tidyverse)
library(egg)
library(tagger) # egg didn't work for the dotplot https://github.com/eliocamp/tagger
library(patchwork)
theme_set(ggsidekick::theme_sleek())
library(methods)

# Globals to set ------------------------------
# include_mpa <- TRUE
include_mpa <- FALSE
#
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0L)
#   stop("This script is meant to be run from the command line.", call. = FALSE)
# include_mpa <- as.logical(as.integer(args[[1]]))

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

.pal <- viridisLite::viridis(4, begin = 0, end = 0.95)
names(.pal) <- c("SYN HS", "SYN QCS", "HBLL OUT N", "SYN WCHG")

# prep index data ----

y1 <- readRDS(file = "data-generated/index-hbll-geo-clean-binomial_gamma.rds") %>%
  # (x/4km left in predict function from density version)*199 sets to sample each grid cell
  # then /1000 to change counts to per 1000 fish
  # work out to  x 0.2
  mutate(est = est * 0.2, lwr = lwr * 0.2, upr = upr * 0.2)
# y1 <- readRDS(file = "data-generated/index-hbll-geo-clean-binomial-gamma.rds") %>%
#   mutate(est = est/10000, lwr = lwr/10000, upr = upr/10000)
y2 <- readRDS(file = "data-generated/index-syn-geo-clean-binomial_gamma.rds")
y <- bind_rows(y1, y2)

mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
filter(y, orig_cv <= 1)
index <- filter(y, orig_cv < 1)

index$species_common_name <- stringr::str_to_title(index$species_common_name)
index <- mutate(index, species_common_name = gsub("Rougheye/Blackspotted Rockfish Complex", "Rougheye/Blackspotted Rockfish", species_common_name))

if (!include_mpa) index <- index %>% filter(type != "MPA only")

# choose labels for each index type ----
restricted_labels <- c("Extrapolated", "Shrunk")
index$type_label <- index$type
index[index$type == "Restricted", ]$type_label <- "Extrapolated"
index[index$type == "Restricted and shrunk", ]$type_label <- "Shrunk"
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
  # "North Pacific Spiny Dogfish",
  # "Big Skate",
  # "Sandpaper Skate", #
  # "Longnose Skate",#HS
  "Spotted Ratfish",
  # "Pacific Cod",
  "Walleye Pollock", #
  # "Pacific Hake",
  # "Sablefish",
  # "Lingcod",
  # "Blackbelly Eelpout",
  "Bocaccio",
  "Canary Rockfish",
  # "Greenstriped Rockfish",#HS
  # "Redbanded Rockfish", # HS
  # # "Redstripe Rockfish",
  # "Rosethorn Rockfish",
  # "Rougheye/Blackspotted Rockfish Complex",#HBLL
  # "Silvergray Rockfish", # an interesting one
  # "Yellowtail Rockfish",
  # "Shortspine Thornyhead",
  "Widow Rockfish",
  "Arrowtooth Flounder", #
  "Curlfin Sole",# QCS
  # "Dover Sole",
  # "English Sole",
  "Flathead Sole",
  # "Pacific Halibut",
  # "Rex Sole",
  # "Slender Sole",
  "Southern Rock Sole"
)

hbll_highlights <- c(
  "North Pacific Spiny Dogfish",
  "Big Skate",
  # "Sandpaper Skate", #
  "Longnose Skate",#HS
  # "Spotted Ratfish",
  # "Pacific Cod",
  # "Walleye Pollock", #
  # "Pacific Hake",
  # "Sablefish",
  "Lingcod",
  # "Blackbelly Eelpout",
  # "Canary Rockfish",
  "China Rockfish",
  # "Greenstriped Rockfish",#HS
  "Quillback Rockfish",
  # "Redbanded Rockfish", # HS
  # # "Redstripe Rockfish",
  # "Rosethorn Rockfish",
  # "Rougheye/Blackspotted Rockfish Complex",#HBLL
  # "Silvergray Rockfish",
  "Tiger Rockfish" # an interesting one
  # "Yelloweye Rockfish"
  # "Yellowtail Rockfish",
  # "Shortspine Thornyhead"
)

# FIGURE 1: indices ----

# split out HBLL to give it a different Y axis scale than the SYN surveys
i1 <- index %>%
  filter(species_common_name %in% hbll_highlights) %>%
  mutate(species_common_name = factor(species_common_name, levels = hbll_highlights)) %>%
  filter((survey_abbrev == "HBLL OUT N"))
i2 <- index %>%
  filter(species_common_name %in% syn_highlights) %>%
  mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
  filter((survey_abbrev == "SYN QCS"))

# spp <- as.data.frame(syn_highlights) %>% rename(species_common_name = syn_highlights)
# i <- left_join(spp, i1) %>% mutate(survey_abbrev = "HBLL OUT N")
#
# g1 <- i %>%
#   filter(species_common_name %in% syn_highlights) %>%
#   mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
#   ggplot(aes(year, est,
#     ymin = lwr, ymax = upr,
#     colour = type_label, fill = type_label, linetype = type_label
#   )) +
#   geom_line(lwd = 0.6) +
#   geom_ribbon(alpha = 0.1, colour = NA) +
#   labs(x = "Year", colour = "Type", fill = "Type", linetype = "Type") +
#   scale_colour_manual(values = colour_pal) +
#   scale_fill_manual(values = colour_pal) +
#   scale_linetype_manual(values = line_pal) +
#   xlim(2005, 2020) +
#   facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") +
#   scale_y_continuous(breaks = waiver(), n.breaks = 4) +
#   ggtitle("Index type:   ") +
#   ylab("Relative abundance in 1000s (HBLL) or biomass in tonnes (SYN)") +
#   theme(
#     plot.title = element_text(hjust = 0.9),
#     strip.text.y = element_blank(),
#     axis.title.x = element_blank(),
#     legend.position = "none"
#   )
#
# g2 <- i2 %>%
# filter(survey_abbrev != "SYN WCHG") %>%
# filter(species_common_name %in% syn_highlights) %>%
# mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
#   ggplot(aes(year, est,
#     ymin = lwr, ymax = upr,
#     colour = type_label, fill = type_label, linetype = type_label
#   )) +
#   geom_line(lwd = 0.6) +
#   geom_ribbon(alpha = 0.1, colour = NA) +
#   labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
#   scale_colour_manual(values = colour_pal) +
#   scale_fill_manual(values = colour_pal) +
#   scale_linetype_manual(values = line_pal) +
#   facet_grid(species_common_name ~ survey_abbrev, scales = "free_y") +
#   scale_y_continuous(breaks = waiver(), n.breaks = 3) +
#   ylab("Relative biomass") +
#   ggtitle("") +
#   theme(
#     legend.justification = c(0, 1), legend.position = c(-0.28, 1.085), legend.direction = "horizontal",
#     strip.text.y = element_text(size = 9, angle = 0, hjust = 0),
#     axis.title.y = element_blank(),
#     axis.title.x = element_text(hjust = 0.2)
#   )
#
# g1 + g2 + patchwork::plot_layout(widths = c(1, 2))

i3 <- bind_rows(i1,i2)%>% ungroup() %>% mutate(
  which_survey = ifelse(survey_abbrev== "HBLL OUT N", "HBLL", "SYN QCS"),
  # spp_survey = paste0(species_common_name, " (", which_survey, ")"),
  spp_survey = paste0(species_common_name),
  spp_survey1 = forcats::fct_reorder(spp_survey, orig_cv))
  # spp_survey1 = forcats::fct_reorder(spp_survey, year, max))
glimpse(i3)

g <- i3  %>%
  ggplot(aes(year, est,
             ymin = lwr, ymax = upr,
             colour = type_label, fill = type_label, linetype = type_label
  )) +
  geom_line(lwd = 0.6) +
  geom_ribbon(alpha = 0.1, colour = NA) +
  labs(x = "Year", colour = " ", fill = " ", linetype = " ") +
  scale_colour_manual(values = colour_pal) +
  scale_fill_manual(values = colour_pal) +
  scale_linetype_manual(values = line_pal) +
  facet_wrap(~spp_survey, scales = "free_y", ncol = 4) +
  # scale_y_continuous(breaks = waiver(), n.breaks = 3) +
  ylab("Relative abundance in 1000s (HBLL) or biomass in tonnes (SYN)") +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  # ggtitle("Index type:   ") +
  theme(
    # legend.justification = c(0, 1), legend.position = c(0.1, 1.095), legend.direction = "horizontal"
    strip.text = element_text(colour = "black"),
    legend.position = "top",axis.text.y = element_text(size = 8)
        )
g


(g <- g + tagger::tag_facets(tag_prefix = "(", position = list(x = 0.1, y = 0.87)))

if (include_mpa) ggsave("figs/index-geo-restricted-highlights.pdf", width = 6.5, height = 8)
if (!include_mpa) ggsave("figs/index-geo-restricted-highlights-noMPA.pdf", width = 8, height = 7)




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

g <- x_long %>%
  filter(survey_abbrev != "SYN WCHG") %>%
  filter(species_common_name %in% syn_highlights) %>%
  mutate(species_common_name = factor(species_common_name, levels = syn_highlights)) %>%
  ggplot(aes(year, re, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.5) +
  geom_line(size = 0.84, alpha = 0.9) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = restricted_cols) +
  ylab("Relative error") +
  xlab("Year") +
  # scale_y_continuous(trans = "S_sqrt", breaks = c(-0.5,-0.1,0, 0.1, 0.5)) +
  # coord_cartesian() + #ylim = c(-0.35, 0.4)
  labs(colour = "") +
  ggtitle("") +
  theme(
    legend.justification = c(0.5, 1), legend.position = c(0.5, 1.065), legend.direction = "horizontal",
    # legend.position = "top", legend.justification=c(0.5,1),
    strip.text.y = element_text(size = 9, angle = 0, hjust = 0)
  )

# with free y-axis scales
(g <- g +
  scale_y_continuous(breaks = waiver(), n.breaks = 5) +
  facet_grid(species_common_name ~ survey_abbrev, scales = "free"))
ggsave("figs/index-geo-restricted-re-highlights.pdf", width = 6.75, height = 8)

# with fixed y-axis scales
(g <- g +
  scale_y_continuous(breaks = waiver(), n.breaks = 3) +
  facet_grid(species_common_name ~ survey_abbrev, scales = "free_x"))
ggsave("figs/index-geo-restricted-re-highlights-fixed.pdf", width = 6.75, height = 8)



# combine precision, accuracy and bias data from all surveys ----
cvdata1 <- readRDS("data-generated/hbll-cv-w-lm-slopes.rds")
cvdata2 <- readRDS("data-generated/syn-cv-w-lm-slopes.rds")

# glimpse(cvdata1)
# glimpse(cvdata2)
cvdata <- bind_rows(cvdata1, cvdata2)
cvdata <- mutate(cvdata, species_common_name = gsub("Rougheye/Blackspotted Rockfish Complex", "Rougheye/Blackspotted Rockfish", species_common_name))


d <- cvdata %>%
  tidyr::pivot_longer(c("cv_ratio", "mare", "slope_re"), names_to = "Response", values_to = "cv_index") %>%
  ungroup() %>%
  filter(cv_index < 1.6) %>%
  mutate(Response = factor(Response, labels = c("CV Ratio", "MARE", "RE trend")))

# change to proportional change in CV
d[d$Response == "CV Ratio",]$cv_index <- d[d$Response == "CV Ratio",]$cv_index - 1

# make RE trend absolute values
d[d$Response == "RE trend",]$cv_index <- abs(d[d$Response == "RE trend",]$cv_index)

# update labels
d$Response <- factor(d$Response, labels = c("CV Ratio - 1", "MARE", "| RE trend |"))


# function for scatterplots ----
plot_scatter <- function(dat, x, y, col_var = "restr_clean",
                         col_vec = restricted_cols,
                         label_vec = restricted_labels,
                         spp = TRUE) {

  g <- ggplot(dat, aes_string(x, y,
    colour = col_var,
    group = "species_common_name"
  )) +
    geom_line(colour = "gray95") +
    geom_point() +
    theme(legend.position = "none", legend.title = element_blank()) +
    scale_colour_manual(values = col_vec, label = label_vec)
# browser()
  if(spp){
    g <- g + ggrepel::geom_text_repel(
        data = filter(dat, `Restriction type` == "re_restr"),
        aes(label = species_common_name),
        colour = "darkgray",
        force = 2, direction = "y", max.overlaps = 4,
        min.segment.length = 10, size = 2
      )
  }
g
}


# HBLL - CV ratio and MARE by status quo CV and prop MPA ----
#
# d2 <- filter(d, survey_abbrev == "HBLL OUT N")
#
# (g1 <- plot_scatter(d2, "cv_orig", "cv_index") +
#   xlab("CV of 'Status quo' index") +
#   guides(shape = "none") +
#   facet_grid(
#     rows = vars(Response),
#     # cols = vars(survey_abbrev),
#     switch = "y",
#     scales = "free_y"
#   ) +
#   theme(
#     axis.title.y = element_blank(),
#     legend.position = c(0.2, 0.93),
#     strip.placement = "outside"
#   ))
#
# (g2 <- plot_scatter(d2, "prop_mpa", "cv_index") +
#   xlab("Biomass proportion inside MPAs") +
#   guides(shape = "none") +
#   facet_wrap(~Response, strip.position = "top", nrow = 2, scales = "free_y") +
#   theme(
#     strip.placement = "outside", strip.text.x = element_blank(),
#     plot.margin = unit(c(0, 0, 0, 0), "cm"),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   ))
#
# g1 <- tag_facet(g1, tag_pool = c("a", "c"), hjust = -0.5, vjust = 2, fontface = 1)
# g2 <- tag_facet(g2, tag_pool = c("b", "d"), hjust = -0.5, vjust = 2, fontface = 1)
#
# g1 + g2 + patchwork::plot_layout(nrow = 1)
#
# ggsave("figs/explore-hbll-cv.pdf", width = 7, height = 7)
#
# # QCS - CV ratio and MARE by status quo CV and prop MPA ----
# d3 <- filter(d, survey_abbrev == "SYN QCS")
#
# (g1 <- plot_scatter(d3, "cv_orig", "cv_index") +
#   xlab("CV of 'Status quo' index") +
#   guides(shape = "none") +
#   facet_grid(
#     rows = vars(Response),
#     # cols = vars(survey_abbrev),
#     switch = "y",
#     scales = "free_y"
#   ) +
#   theme(
#     axis.title.y = element_blank(),
#     # legend.position = c(0.2,0.95),
#     strip.placement = "outside"
#   ))
#
# (g2 <- plot_scatter(d3, "prop_mpa", "cv_index") +
#   xlab("Biomass proportion inside MPAs") +
#   guides(shape = "none") +
#   facet_wrap(~Response, strip.position = "top", nrow = 2, scales = "free_y") +
#   theme(
#     strip.placement = "outside", strip.text.x = element_blank(),
#     plot.margin = unit(c(0, 0, 0, 0), "cm"),
#     legend.position = c(0.25, 0.97),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   ))
#
# g1 <- tag_facet(g1, tag_pool = c("a", "c"), hjust = -0.5, vjust = 2, fontface = 1)
# g2 <- tag_facet(g2, tag_pool = c("b", "d"), hjust = -0.5, vjust = 2, fontface = 1)
#
# g1 + g2 + patchwork::plot_layout(nrow = 1)
#
# ggsave("figs/explore-qcs-cv.pdf", width = 7, height = 7)
#

##  CV ratio and MARE by status quo CV for all surveys at once ----
if (length(unique(d$survey_abbrev)) > 3) { # makes sure all surveys

  (g0 <- plot_scatter(d, "cv_orig", "cv_index", spp=FALSE) +
    xlab("CV of 'Status quo' index") +
    guides(shape = "none") +
    facet_grid(
      rows = vars(Response),
      cols = vars(survey_abbrev),
      switch = "y",
      scales = "free_y"
    ) +
    theme(
      axis.title.y = element_blank(),
      # legend.position = c(0.1, 0.95),
      strip.placement = "outside"
    ))
  (g0 <- tag_facet(g0, fontface = 1))

  # ggsave("figs/explore-all-cv-by-cv.pdf", width = 6, height = 3) #for two rows
  ggsave("figs/explore-all-cv-by-cv.pdf", width = 6, height = 4.5) # for 3 rows


  # CV ratio and MARE by prop MPA for all surveys at once for appendix ----
  (g <- plot_scatter(d, "prop_mpa", "cv_index", spp=FALSE) +
    xlab("Biomass proportion inside MPAs") +
    guides(shape = "none") +
    facet_grid(
      rows = vars(Response),
      cols = vars(survey_abbrev),
      switch = "y",
      scales = "free_y"
    ) +
    theme(
      axis.title.y = element_blank(),
      # legend.position = c(0.1, 0.95),
      strip.placement = "outside"
    ))
  (g <- tag_facet(g, fontface = 1))

  # ggsave("figs/explore-all-cv-by-mpa.pdf", width = 6, height = 3)
  ggsave("figs/explore-all-cv-by-mpa.pdf", width = 6, height = 4.5) # for 3 rows


  # CV ratio and MARE by prop MPA for all surveys at once----
  d_shrunk <- filter(d, `Restriction type` == "re_shrunk")
  (g <- g <- ggplot(d_shrunk,
                    aes_string("prop_mpa", "cv_index", colour = "survey_abbrev")) +
     geom_point(alpha = 0.8) +
     xlab("Biomass proportion inside MPAs") +
     guides(shape = "none") +
      scale_colour_manual(name = "Survey", values = .pal) +
     facet_grid(
       rows = vars(Response),
       # cols = vars(survey_abbrev),
       switch = "y",
       scales = "free_y"
     ) +
     theme(
       axis.title.y = element_blank(),
       # legend.position = c(0.1, 0.95),
       strip.placement = "outside"
     ))
  (g <- tag_facet(g, fontface = 1))

  # ggsave("figs/explore-all-cv-by-mpa.pdf", width = 6, height = 3)
  ggsave("figs/explore-all-cv-by-mpa2.pdf", width = 4.5, height = 5) # for 3 rows


  (g <- g <- ggplot(d_shrunk,
                    aes_string("cv_orig", "cv_index", colour = "survey_abbrev")) +
      geom_point(alpha = 0.8) +
      xlab("CV of 'Status quo' index") +
      guides(shape = "none") +
      scale_colour_manual(name = "Survey", values = .pal) +
      facet_grid(
        rows = vars(Response),
        # cols = vars(survey_abbrev),
        switch = "y",
        scales = "free_y"
      ) +
      theme(
        axis.title.y = element_blank(),
        # legend.position = c(0.1, 0.95),
        strip.placement = "outside"
      ))
  (g <- tag_facet(g, fontface = 1))
  ggsave("figs/explore-all-cv-by-cv2.pdf", width = 4.5, height = 5) # for 3 rows

}


# MARE by CV ratio ----

cvdata2 <- cvdata
cvdata2[cvdata2$cv_ratio > 1.5, ]$cv_ratio <- 1.5

cvdata2$cv_ratio <- cvdata2$cv_ratio - 1

(g <- plot_scatter(cvdata2, "cv_ratio", "mare") +
  xlab("CV ratio - 1") +
  ylab("MARE") +
  guides(shape = "none") +
  facet_wrap(~survey_abbrev,
    ncol = 2,
    # cols = vars(survey_abbrev),
    # switch = "y",
    scales = "free_x"
  ) +
  theme(
    legend.position = c(0.12, 0.96),
    strip.placement = "outside"
  ))

g <- tag_facet(g, hjust = -0.5, vjust = 2, fontface = 1)
# g
ggsave("figs/explore-all-mare-by-cv-ratio.pdf", width = 7, height = 7)


# FIGURE 6: slopes ----
#
(g <- plot_scatter(cvdata, "slope_mpa", "slope_re") +
  facet_wrap(~survey_abbrev,
    scales = "free"
  ) +
  ylab("Change in RE per decade") +
  xlab("Change in proportion of biomass inside MPAs") +
  geom_hline(yintercept = 0, colour = "gray80") +
  geom_vline(xintercept = 0, colour = "gray70") +
  theme(legend.position = c(0.37, 0.95)))
(g <- tag_facet(g, hjust = -0.5, vjust = 2, fontface = 1))

ggsave("figs/explore-all-slopes.pdf", width = 8, height = 8)

(g <- plot_scatter(cvdata, "prop_mpa", "abs(slope_re)") +
  ylab("Absolute change in RE per decade") +
  xlab("Proportion of biomass inside MPAs") +
  facet_wrap(~survey_abbrev,
    scales = "free"
  ) + theme(legend.position = c(0.08, 0.96)))
ggsave("figs/explore-abs-slope.pdf", width = 7, height = 7)


d_shrunk2 <- filter(cvdata, `Restriction type` == "re_shrunk")
(g <- g <- ggplot(d_shrunk2,
                  aes_string("slope_mpa", "slope_re", colour = "survey_abbrev")) +
    geom_point(alpha = 0.8) +
    ylab("Change in RE per decade") +
    xlab("Change in proportion of biomass inside MPAs") +
    guides(shape = "none") +
    geom_hline(yintercept = 0, colour = "gray80") +
    geom_vline(xintercept = 0, colour = "gray70") +
    scale_colour_manual(name = "Survey", values = .pal)+
    theme(legend.position = c(0.8, 0.2))+
    ggrepel::geom_text_repel(
      aes(label = species_common_name),
      colour = "darkgray",
      force = 2, direction = "y", max.overlaps = 4,
      min.segment.length = 10, size = 2.5
    ))

ggsave("figs/explore-all-slopes2.pdf", width = 6, height = 6)

# ### benefits of interpolation
#
d <- cvdata %>% group_by(species_common_name, survey_abbrev) %>% mutate(
  temp_bias_interp_ratio = abs(slope_re[restr_clean == "Same survey domain"])/abs(slope_re[restr_clean == "Shrunk survey domain"]),
  cv_interp_ratio = (cv_mean[restr_clean == "Same survey domain"])/(cv_mean[restr_clean == "Shrunk survey domain"]),
  mare_interp_ratio = (mare[restr_clean == "Same survey domain"])/(mare[restr_clean == "Shrunk survey domain"])
)

median(d$temp_bias_interp_ratio)
median(d$cv_interp_ratio)
median(d$mare_interp_ratio)


# SUPPLEMENTAL FIGURES ----

# INDEX PLOTS for each survey ----
g <- index %>%
  filter(survey_abbrev == "HBLL OUT N") %>%
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  scale_colour_manual(values = colour_pal) +
  scale_fill_manual(values = colour_pal) +
  scale_linetype_manual(values = line_pal) +
  ylab("Relative abundance in 1000s") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  theme(legend.position = "top")
g
ggsave("figs/index-hbll-geo-restricted.pdf", width = 9, height = 8)

g <- index %>%
  filter(survey_abbrev == "SYN QCS") %>%
  # CI on exrapolated Yellowmouth index blows up,
  # truncating here to see uncertainty on other indices
  # mutate(upr = ifelse(
  #   species_common_name == "Yellowmouth Rockfish" & upr > 10000, 10000, upr
  # )) %>%
  # filter(species_common_name == "Southern Rock Sole") %>%
  # filter(species_common_name == "Dover Sole") %>%
  # filter(type != "Restricted and shrunk") %>%
  # filter(type != "Restricted") %>%
  ggplot(aes(year, est,
    ymin = lwr, ymax = upr,
    colour = type_label, fill = type_label, linetype = type_label
  )) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  # geom_ribbon(alpha = 0.2, fill = NA) +
  coord_cartesian() +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  scale_colour_manual(values = colour_pal) +
  scale_fill_manual(values = colour_pal) +
  scale_linetype_manual(values = line_pal) +
  ylab("Relative biomass in tonnes") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4)+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-qcs-geo-restricted.pdf", width = 9.1, height = 11, limitsize = FALSE)

g <- index %>%
  filter(survey_abbrev == "SYN HS") %>%
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  scale_colour_manual(values = colour_pal) +
  scale_fill_manual(values = colour_pal) +
  scale_linetype_manual(values = line_pal) +
  ylab("Relative biomass in tonnes") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4)+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-hs-geo-restricted.pdf", width = 9, height = 9, limitsize = FALSE)

g <- index %>%
  filter(survey_abbrev == "SYN WCHG") %>%
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type_label, fill = type_label, linetype = type_label)) +
  geom_line(lwd = 0.9) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  labs(x = "Year", colour = "Index type", fill = "Index type", linetype = "Index type") +
  scale_colour_manual(values = colour_pal) +
  scale_fill_manual(values = colour_pal) +
  scale_linetype_manual(values = line_pal) +
  ylab("Relative biomass in tonnes") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4)+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-wchg-geo-restricted.pdf", width = 9.2, height = 6, limitsize = FALSE)


# RE PLOTS for each survey ----
g <- x_long %>%
  filter(survey_abbrev == "HBLL OUT N") %>%
  ggplot(aes(year, re, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(size = 0.8, alpha = 0.8) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  ylab("Relative error") +
  xlab("Year") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  # scale_y_continuous(trans = "S_sqrt", breaks = c(-0.5,-0.1,0, 0.1, 0.5)) +
  # coord_cartesian(ylim = c(-0.35, 0.4)) +
  labs(colour = "")+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 9, height = 8)

g <- x_long %>%
  filter(survey_abbrev == "SYN QCS") %>%
  ggplot(aes(year, re, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(size = 0.8, alpha = 0.8) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  ylab("Relative error") +
  xlab("Year") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  labs(colour = "")+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-qcs-geo-restricted-re.pdf", width = 9.1, height = 11, limitsize = FALSE)

g <- x_long %>%
  filter(survey_abbrev == "SYN HS") %>%
  ggplot(aes(year, re, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(size = 0.8, alpha = 0.8) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  ylab("Relative error") +
  xlab("Year") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  labs(colour = "")+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-hs-geo-restricted-re.pdf", width = 9, height = 9, limitsize = FALSE)

g <- x_long %>%
  filter(survey_abbrev == "SYN WCHG") %>%
  ggplot(aes(year, re, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_line(size = 0.8, alpha = 0.8) +
  # scale_color_brewer(palette = "Set2") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  ylab("Relative error") +
  xlab("Year") +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  labs(colour = "")+
  theme(legend.position = "top",axis.text.y = element_text(size = 7))
g
ggsave("figs/index-wchg-geo-restricted-re.pdf", width = 9.2, height = 6, limitsize = FALSE)


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
    names_to = "Restriction type", values_to = "CV ratio"
  )
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
  ungroup() %>%
  group_by(species_common_name) %>%
  mutate(est_avg = mean(est, na.rm = TRUE)) %>%
  left_join(lu_cv) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est_avg), est,
    colour = restr_clean, ymin = lwr, ymax = upr
  )) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.75), size = 0.35) +
  coord_flip(ylim = c(0.95, 1.8)) +
  scale_y_continuous(trans = "log") +
  xlab("") +
  ylab("Ratio of index CV (restricted/status quo)") +
  labs(colour = " ") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90")) +
  facet_wrap(~survey_abbrev, ncol = 4, scales = "free_x")
g
# (g <- tag_facet_outside(g, fontface = 1))

ggsave("figs/index-geo-cv-ratio-dotplot.pdf", width = 8.5, height = 6)


# CV percent change dotplot ----
cv2 <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(
    cv_change_restr = round((cv[type == "Restricted"]-cv[type == "Status quo"]) /
      cv[type == "Status quo"], 2),
    cv_change_shrunk = round((cv[type == "Restricted and shrunk"]-cv[type == "Status quo"]) /
      cv[type == "Status quo"], 2)
  )
cv_long2 <- cv2 %>%
  tidyr::pivot_longer(starts_with("cv"),
                      names_to = "Restriction type", values_to = "CV change"
  )
cv_long2 %>%
  group_by(`Restriction type`, survey_abbrev) %>%
  summarise(mean_ratio = mean(`CV change`)) %>%
  knitr::kable(digits = 2)
lu_cv2 <- tibble(
  "Restriction type" = c("cv_change_restr", "cv_change_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)
g <- cv_long2 %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV change`), upr = max(`CV change`), est = mean(`CV change`)) %>%
  ungroup() %>%
  group_by(species_common_name) %>%
  mutate(est_avg = mean(est, na.rm = TRUE)) %>%
  left_join(lu_cv2) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est_avg), est,
             colour = restr_clean, ymin = lwr, ymax = upr
  )) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.75), size = 0.35) +
  coord_flip() +
  # scale_y_continuous(trans = "log") +
  xlab("") +
  ylab("% change in CV with restriction") +
  labs(colour = " ") +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90")) +
  facet_wrap(~survey_abbrev, ncol = 4, scales = "free_x")
g
# (g <- tag_facet_outside(g, fontface = 1))

ggsave("figs/index-geo-cv-change-dotplot.pdf", width = 8.5, height = 6)


# MARE dotplot ----
g <- x_long %>% # first created for FIG 2
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
  left_join(lu) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), -est), est,
    colour = restr_clean
  )) +
  geom_point(size = 2, position = position_dodge(width = 0)) +
  xlab("") +
  ylab("Median absolute relative error (MARE) of restricted index compared to status quo") +
  labs(colour = "") +
  coord_flip() +
  scale_y_continuous(breaks = waiver(), n.breaks = 4) +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey90")) +
  # guides(colour = guide_legend(nrow = 2L))+
  facet_wrap(~survey_abbrev, ncol = 4)
g
ggsave("figs/index-geo-mare-dotplot.pdf", width = 9, height = 8)

# RAW data checks ----
# dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
# ggplot(dat_to_fit, aes(year, hook_count)) + geom_jitter(alpha = 0.1)
# ggplot(dat_to_fit, aes(hook_count, area_km2)) + geom_jitter(alpha = 0.1) + facet_wrap(~year)

dd1 <- cv_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
  ungroup() %>%
  group_by(species_common_name) %>%
  mutate(
    # est_avg = mean(est, na.rm = TRUE),
         measure = "CV ratio (lost precision)") %>%
  left_join(lu_cv)

dd1b <- cv_long2 %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = quantile(`CV change`, 0.025),
            upr = quantile(`CV change`, 0.975),
            est = mean(`CV change`)) %>%
  ungroup() %>%
  group_by(species_common_name) %>%
  mutate(
    est_avg = mean(est, na.rm = TRUE),
    measure = "CV ratio - 1 (precision loss)") %>%
  left_join(lu_cv2)

dd2 <-  x_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = quantile(abs(re), 0.025),
            upr = ifelse(quantile(abs(re), 0.975) > 0.5, 0.5, quantile(abs(re), 0.975)),
            est = median(abs(re))) %>%
  mutate(
    est_avg = mean(est, na.rm = TRUE),
    measure = "MARE (accuracy loss)")

dd3 <- d %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(
    lwr = (median(slope_re)-1.98*mean(se_slope_re))/1,
    upr = (median(slope_re)+1.98*mean(se_slope_re))/1,
    est = median((slope_re))/1) %>%
  mutate(
    est_avg = mean(abs(est), na.rm = TRUE)/1,
    measure = "RE trend (bias)") %>%
  select(
    survey_abbrev, species_common_name, `Restriction type`,
    lwr, upr,
    est_avg,
    est, measure
    )

dd <- bind_rows(dd2, dd3)%>%
  left_join(lu) %>% bind_rows(
    # dd1
    dd1b
    )

g <- dd %>% filter(!survey_abbrev %in% c("SYN HS", "SYN WCHG")) %>%
  ggplot(aes(
    forcats::fct_reorder(stringr::str_to_title(species_common_name), -est_avg, mean, na.rm = TRUE),
    est, ymin = lwr, ymax = upr,
             colour = as.factor(restr_clean)
  )) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.75), size = 0.35) +
  # xlab("") +
  # ylab("") +
  coord_flip() +
  scale_y_continuous(breaks = waiver(), n.breaks = 4, expand = c(0,0)) +
  scale_colour_manual(values = restricted_cols, label = restricted_labels) +
  labs(x = "", y = "", colour = "Index type", fill = "Index type", linetype = "Index type") +
  theme(
    strip.text = element_text(colour = "black"),
    legend.position = "top", panel.grid.major.y = element_line(colour = "grey90"),
    strip.placement = "outside") +
  facet_grid(survey_abbrev~measure, scales = "free",
             space="free_y", switch="x")
g
# not sure why but egg didn't work here
# devtools::install_github("eliocamp/tagger")
(g <- g + tagger::tag_facets(tag_prefix = "(", position = list(x = 0.1, y = 0.96)))

ggsave("figs/index-geo-combind-dotplot.pdf", width = 7.8, height = 8)

