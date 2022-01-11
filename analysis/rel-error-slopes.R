# get slopes of relative error
library(dplyr)
library(ggplot2)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

# survey <- "HBLL"
survey <- "SYN"


# for now using delta-gamma?
if (survey == "HBLL") y <- readRDS(file = "data-generated/index-hbll-geo-clean.rds")
if (survey == "SYN") y <- readRDS(file = "data-generated/index-syn-geo-clean.rds")


# keep only species with original cv <= 1
mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
filter(y, orig_cv <= 1)
index <- filter(y, orig_cv < 1)
# index <- filter(index, cv < 2)

index$species_common_name <- stringr::str_to_title(index$species_common_name)

# get cv ratios and prop mpa
cv2 <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(
    prop_mpa = est[type == "MPA only"]/est[type == "Status quo"],
    cv_ratio_restr = cv[type == "Restricted"] /
      cv[type == "Status quo"],
    cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
      cv[type == "Status quo"]
  )

cv_long <- cv2 %>% select(-prop_mpa) %>%
  tidyr::pivot_longer(starts_with("cv"),
                      names_to = "Restriction type", values_to = "CV ratio")

prop_mpa <- cv2 %>% select(species_common_name, survey_abbrev, year, prop_mpa) %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(
    prop_mpa = mean(prop_mpa)
  )

lu <- tibble(
  "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

cvratio <- cv_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
  left_join(lu)



# get raw cvs
lu <- tibble(
  "type" = c("Restricted", "Restricted and shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

cvraw <- index %>% filter(type !="MPA only") %>%
  group_by(species_common_name, survey_abbrev, type) %>%
  summarise(
    cv_orig = mean(orig_cv),
    cv_max = max(cv),
    cv_mean = mean(cv)
  ) %>% filter(type !="Status quo") %>% left_join(., lu)


# bad year in SYN, removed from plot so removing here too
if (survey == "SYN") {
  index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
}

re <- index %>%
  group_by(species_common_name, survey_abbrev, type) %>%
  mutate(est = est / exp(mean(log(est)))) %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(
    prop_mpa = est[type == "MPA only"]/est[type == "Status quo"],
    re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
    re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
  )

lu <- tibble(
  "Restriction type" = c("re_restr", "re_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

re_long <- re %>%
  tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re") %>%
  left_join(lu)

mare <- re_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
  left_join(lu)


# rename variables for other previously calculated indices
cvratio2 <- cvratio %>% rename(cv_ratio = est, cv_lwr =lwr, cv_upr = upr) %>% select(-`Restriction type`)
mare2 <- mare %>% rename(mare = est, mare_lwr = lwr, mare_upr = upr) %>%
  mutate(spp_by_survey = paste(species_common_name, survey_abbrev))



# mean(re_long$year)
re_long <- re_long %>%
  mutate(decade = (year - 2012)/10,
         re100 = re*100,
         spp_by_survey = paste(species_common_name,survey_abbrev)
         )

survs <- unique(re_long$survey_abbrev)

re_long %>% filter(`Restriction type` == "re_restr") %>%
  ggplot(aes(year, re100, colour = species_common_name)) +
  geom_line() + facet_wrap(~survey_abbrev, dir= "v")

index %>% filter(type == "Restricted and shrunk") %>%
  ggplot(aes(year, cv, colour = species_common_name)) +
  geom_line() + facet_wrap(~survey_abbrev, dir= "v")

index %>% filter(type == "Status quo") %>%
  ggplot(aes(year, cv, colour = species_common_name)) +
  geom_line() + facet_wrap(~survey_abbrev, dir= "v")


# # calculate slope from fig 2 relative error plots using random effects
# library(lme4)
#
# by_surv <- list()
# for (j in seq_along(survs)){
#
# restr <- re_long %>% filter(survey_abbrev == survs[j]) %>%
#   filter(`Restriction type` == "re_restr")
# m_restr <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = restr)
# m_restr
#
# shrunk <- re_long %>% filter(survey_abbrev == survs[j]) %>%
#   filter(`Restriction type` == "re_shrunk")
# m_shrunk  <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = shrunk)
# m_shrunk
#
# # These values are a combination of the fixed effects and the variance components
# restr_coefs <- coef(m_restr)$species_common_name %>%
#   tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
# restr_coefs$restr_clean <- "Same survey domain"
# restr_coefs$survey_abbrev <- survs[j]
#
# shrunk_coefs <- coef(m_shrunk)$species_common_name %>%
#   tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
# shrunk_coefs$restr_clean <- "Shrunk survey domain"
# shrunk_coefs$survey_abbrev <- survs[j]
#
# by_surv[[j]] <- bind_rows(restr_coefs, shrunk_coefs)
# }
# coefs <- do.call(rbind, by_surv)
#
# # join everything together
# cvdata <- left_join(cvratio2, mare2) %>% left_join(., coefs) %>%
#   left_join(., cvraw) %>% left_join(., prop_mpa)
#
# if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-re-slopes.rds")
# if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-re-slopes.rds")



# OR use simple linear models

survs <- unique(re_long$survey_abbrev)
spp <- unique(re_long$species_common_name)
sps <- list()


for (i in seq_along(spp)){
  for (j in seq_along(survs)){
  dat <- re_long %>% filter(survey_abbrev == survs[j]) %>%
    filter(species_common_name == spp[i])

  d1 <- filter(dat, `Restriction type` == "re_restr")
  m1  <- tryCatch(lm(re100 ~ 1 + decade, data = d1), error=function(err) NA)

  d2 <- filter(dat, `Restriction type` == "re_shrunk")
  m2  <- tryCatch(lm(re100 ~ 1 + decade, data = d2), error=function(err) NA)

  m3  <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d1), error=function(err) NA)
  m4  <- tryCatch(lm(prop_mpa ~ 1 + decade, data = d2), error=function(err) NA)

  k <- length(spp)*(j-1) + i

  sps[[k]] <- tibble(
    survey_abbrev = c(survs[j], survs[j]),
    species_common_name = c(spp[i],spp[i]),
    `(Intercept)` = c(tryCatch(m1$coefficients[[1]], error=function(err) NA),
                      tryCatch(m2$coefficients[[1]], error=function(err) NA)),
    slope_re = c(tryCatch(m1$coefficients[["decade"]], error=function(err) NA),
              tryCatch(m2$coefficients[["decade"]], error=function(err) NA)),
    slope_mpa = c(tryCatch(m3$coefficients[["decade"]], error=function(err) NA),
              tryCatch(m4$coefficients[["decade"]], error=function(err) NA)),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )
  }
}
coefs <- do.call(rbind, sps)

# join everything together
cvdata <- left_join(cvratio2, mare2) %>% left_join(., coefs) %>%
  left_join(., cvraw) %>% left_join(., prop_mpa)

if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-lm-slopes.rds")
if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-lm-slopes.rds")



# exploratory plots
# cvdata <- readRDS("data-generated/hbll-cv-w-re-slopes.rds")

# combine all surveys
cvdata1 <- readRDS("data-generated/hbll-cv-w-lm-slopes.rds")
cvdata2 <- readRDS("data-generated/syn-cv-w-lm-slopes.rds")

glimpse(cvdata1)
glimpse(cvdata2)

cvdata <- bind_rows(cvdata1, cvdata2)

#
# if (survey == "SYN") {
# plot_scatter <- function(dat, x, y) {
#   ggplot(dat, aes_string(x, y,
#                          colour = "survey_abbrev", # this upends the way colour was used before
#                          shape = "restr_clean",
#                          # colour = "restr_clean",
#                          # shape = "survey_abbrev",
#                          group = "species_common_name")) +
#   # geom_line(colour = "lightgray") +
#   ggrepel::geom_text_repel(aes(label = species_common_name),
#                            colour = "darkgray",
#                            force = 2, direction = "both", max.overlaps = 3,
#                            min.segment.length = 10, size = 2,
#                            data = filter(dat, `Restriction type` == "re_restr")) +
#   geom_point() +
#   scale_color_brewer(palette = "Set2") +
#   # # might be worth combining with HBLL for 4 survey facets?
#   # facet_wrap(~survey_abbrev, dir = "v") +
#   theme(legend.position = "none")
# }
# } else {
  plot_scatter <- function(dat, x, y) {
    ggplot(dat, aes_string(x, y,
                           colour = "restr_clean",
                           group = "species_common_name")) +
      geom_line(colour = "gray95") +
      ggrepel::geom_text_repel(data = filter(dat, `Restriction type` == "re_restr"),
                               aes(label = species_common_name),
                               colour = "darkgray",
                               force = 2, direction = "y", max.overlaps = 5,
                               min.segment.length = 10, size = 2) +
      geom_point() +
      theme(legend.position = "none", legend.title = element_blank()) +
      scale_colour_manual(values = c("#FC8D62", "#8DA0CB"), label=c("Extrapolated", "Shrunk"))
      # scale_color_brewer(palette = "Set2")
  }
# }


#### fig 1
# cv_ratio = precision - shrunk reduces loss of precision
# mare = accuracy of mean - uncertain becomes more uncertain

d <- cvdata %>%
  tidyr::pivot_longer(c("cv_ratio", "mare"), names_to = "Response", values_to = "cv_index") %>%
  ungroup() %>%
  filter(cv_index < 1.6) %>%
  mutate(Response = factor(Response, labels = c("CV Ratio", "MARE")))

# plots of just HBLL
d2 <- filter(d, survey_abbrev == "HBLL OUT N")

(g1 <- plot_scatter(d2, "cv_orig", "cv_index") +
    xlab("CV of 'Status quo' index") +
    guides(shape = "none") +
    facet_grid(rows=vars(Response),
               # cols = vars(survey_abbrev),
               switch = "y",
               scales = "free_y") +
    theme(axis.title.y = element_blank(),
          legend.position = c(0.2,0.95),
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

g1 + g2 + patchwork::plot_layout(nrow = 1)
# if (survey == "HBLL")
ggsave("figs/explore-hbll-cv.pdf", width = 7, height = 7)
# if (survey == "SYN") ggsave("figs/explore-syn-cv.pdf", width = 9, height = 9)
#guides = "collect", design = layout


# plots of just QCS
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
          legend.position = c(0.2,0.95),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))

g1 + g2 + patchwork::plot_layout(nrow = 1)

ggsave("figs/explore-qcs-cv.pdf", width = 7, height = 7)



# plot all surveys at once for appendix
if(length(unique(d$survey_abbrev))> 3){

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

ggsave("figs/explore-all-cv-by-mpa.pdf", width = 10, height = 5.5)

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
ggsave("figs/explore-all-cv-by-cv.pdf", width = 10, height = 5.5)

}


# fig 2
# good - once you get a slope of prop mpa overtime it predicts a given bias
# bottomright = increasing hiding sp, so negative bias over time - fishing?, by chance = local climate velocity
# topleft = increasing species, so positive bias
# - what is ultimately expected to happen for species target for protection in mpas

d_mare <- filter(d, Response == "MARE")
(g <- plot_scatter(d_mare, "slope_mpa", "slope_re/100") +
    facet_wrap(~survey_abbrev,
      # rows = vars(survey_abbrev), # switch = "y",
               scales = "free") +
    ylab("Change in RE per decade") +
    xlab("Change in proportion of biomass inside MPAs") +
    geom_hline(yintercept = 0, colour = "gray80") + geom_vline(xintercept = 0, colour = "gray70"))

ggsave("figs/explore-all-slopes.pdf", width = 7, height = 7)

(g <- plot_scatter(d_mare, "prop_mpa", "abs(slope_re/100)") +
    ylab("Absolute change in RE per decade") +
    xlab("Proportion of biomass inside MPAs") +
    facet_wrap(~survey_abbrev,# rows = vars(survey_abbrev), # switch = "y",
               scales = "free") + theme(legend.position = c(0.08,0.96)))
ggsave("figs/explore-abs-slope.pdf", width = 7, height = 7)


### benefits of interpolation

d <- cvdata %>% group_by(species_common_name, survey_abbrev) %>% mutate(
  temp_bias_interp_ratio = abs(slope_re[restr_clean == "Same survey domain"])/abs(slope_re[restr_clean == "Shrunk survey domain"]),
  cv_interp_ratio = (cv_mean[restr_clean == "Same survey domain"])/(cv_mean[restr_clean == "Shrunk survey domain"]),
  mare_interp_ratio = (mare[restr_clean == "Same survey domain"])/(mare[restr_clean == "Shrunk survey domain"])
)

median(d$temp_bias_interp_ratio)
median(d$cv_interp_ratio)
median(d$mare_interp_ratio)

# (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
# (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
# (g <- ggplot(d) + geom_point(aes_string("prop_mpa", "mare_interp_ratio"))+ geom_hline(yintercept = 1))
#
# (g <- ggplot(d) + geom_point(aes_string("cv_orig", "temp_bias_interp_ratio"))+ geom_hline(yintercept = 1))
# (g <- ggplot(d) + geom_point(aes_string("cv_orig", "cv_interp_ratio"))+ geom_hline(yintercept = 1))
# (g <- ggplot(d) + geom_point(aes_string("cv_orig", "mare_interp_ratio"))+ geom_hline(yintercept = 1))

