# get slopes of relative error
library(dplyr)
library(ggplot2)
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

survey <- "HBLL"
survey <- "SYN"


# for now using delta-gamma?
if (survey == "HBLL") y <- readRDS(file = "data-generated/index-hbll-geo-clean.rds")
if (survey == "SYN") y <- readRDS(file = "data-generated/index-syn-geo-clean.rds")


# keep only species with original cv <= 1
mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
filter(y, orig_cv <= 1)
index <- filter(y, orig_cv < 1)

index$species_common_name <- stringr::str_to_title(index$species_common_name)

# get raw cvs
lu <- tibble(
  "type" = c("Restricted", "Restricted and shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

cvraw <- index %>%
  group_by(species_common_name, survey_abbrev, type) %>%
  summarise(
    cv_orig = mean(orig_cv),
    cv_max = max(cv),
    cv_mean = mean(cv)
  ) %>% filter(type != "Status quo") %>% left_join(., lu)


# get cv ratios
cv2 <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(
    cv_ratio_restr = cv[type == "Restricted"] /
      cv[type == "Status quo"],
    cv_ratio_shrunk = cv[type == "Restricted and shrunk"] /
      cv[type == "Status quo"]
  )

cv_long <- cv2 %>%
  tidyr::pivot_longer(starts_with("cv"),
                      names_to = "Restriction type", values_to = "CV ratio")

lu <- tibble(
  "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)
cvratio <- cv_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
  left_join(lu)


# bad year in SYN, removed from plot so removing here too
if (survey == "SYN") {
  index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
}

re <- index %>%
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



# calculate slope from fig 2 relative error plots using random effects
library(lme4)

# mean(re_long$year)
re_long <- re_long %>% mutate(decade = (year - 2012)/10,
                              re100 = re*100,
                              spp_by_survey = paste(species_common_name, survey_abbrev)
                              )

restr <- filter(re_long, `Restriction type` == "re_restr")


m_restr <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = restr)
m_restr

shrunk <- filter(re_long, `Restriction type` == "re_shrunk")
m_shrunk  <- lmer(re100 ~ 1 + decade + (decade|species_common_name), data = shrunk)
m_shrunk

# These values are a combination of the fixed effects and the variance components
restr_coefs <- coef(m_restr)$species_common_name  %>%
  tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
restr_coefs$restr_clean <- "Same survey domain"

shrunk_coefs <- coef(m_shrunk)$species_common_name %>%
  tibble::rownames_to_column(., "species_common_name") %>% rename(slope = decade)
shrunk_coefs$restr_clean <- "Shrunk survey domain"

coefs <- bind_rows(restr_coefs, shrunk_coefs)

# join everything together
cvdata <- left_join(cvratio2, mare2) %>% left_join(., coefs) %>% left_join(., cvraw)

if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-re-slopes.rds")
if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-re-slopes.rds")



# OR use simple linear models

survs <- unique(re_long$survey_abbrev)
spp <- unique(re_long$species_common_name)
sps <- list()


for (i in seq_along(spp)){
  for (j in seq_along(survs)){
  dat <- re_long %>% filter(survey_abbrev == survs[j]) %>% filter(species_common_name == spp[i])

  d1 <- filter(dat, `Restriction type` == "re_restr")
  m1  <- tryCatch(lm(re100 ~ 1 + decade, data = d1), error=function(err) NA)

  d2 <- filter(dat, `Restriction type` == "re_shrunk")
  m2  <- tryCatch(lm(re100 ~ 1 + decade, data = d2), error=function(err) NA)


  k <- length(spp)*(j-1) + i

  sps[[k]] <- tibble(
    survey_abbrev = c(survs[j], survs[j]),
    species_common_name = c(spp[i],spp[i]),
    `(Intercept)` = c(tryCatch(m1$coefficients[[1]], error=function(err) NA),
                      tryCatch(m2$coefficients[[1]], error=function(err) NA)),
    slope = c(tryCatch(m1$coefficients[["decade"]], error=function(err) NA),
              tryCatch(m2$coefficients[["decade"]], error=function(err) NA)),
    restr_clean = c("Same survey domain", "Shrunk survey domain")
  )
  }
}
coefs <- do.call(rbind, sps)

# join everything together
cvdata <- left_join(cvratio2, mare2) %>% left_join(., coefs) %>% left_join(., cvraw)

if (survey == "HBLL") saveRDS(cvdata, "data-generated/hbll-cv-w-lm-slopes.rds")
if (survey == "SYN") saveRDS(cvdata, "data-generated/syn-cv-w-lm-slopes.rds")



# exploratory plots
leg <- cvdata %>% filter(`Restriction type` == "re_restr")

if (survey == "SYN") {
plot_scatter <- function(dat, x, y) {
  ggplot(dat, aes_string(x, y,
                         colour = "survey_abbrev", # this upends the way colour was used before
                         shape = "restr_clean",
                         # colour = "restr_clean",
                         # shape = "survey_abbrev",
                         group = "species_common_name")) +
  # geom_line(colour = "lightgray") +
  ggrepel::geom_text_repel(aes(label = species_common_name), colour = "darkgray",
                           force = 2, direction = "both", max.overlaps = 3,
                           min.segment.length = 10, size = 2,
                           data = leg) +
  geom_point() +
  scale_color_brewer(palette = "Set2") +
  # # might be worth combining with HBLL for 4 survey facets?
  # facet_wrap(~survey_abbrev, dir = "v") +
  theme(legend.position = "none")
}
} else {
  plot_scatter <- function(dat, x, y) {
    ggplot(dat, aes_string(x, y,
                           shape = "survey_abbrev",
                           colour = "restr_clean",
                           group = "species_common_name")) +
      geom_line(colour = "lightgray") +
      ggrepel::geom_text_repel(aes(label = species_common_name), colour = "darkgray",
                               force = 2, direction = "both", max.overlaps = 3,
                               min.segment.length = 10, size = 2,
                               data = leg) +
      geom_point() +
      scale_color_brewer(palette = "Set2") +
      theme(legend.position = "none")
  }

}

g1 <- plot_scatter(cvdata, "mare", "slope") + theme(legend.position = c(0.4,0.2), legend.title = element_blank())
g2 <- plot_scatter(cvdata, "cv_orig", "slope")
## this has stonger outliers, so keeping mean instead
# (g3 <- plot_scatter(cvdata, "cv_max", "slope"))
g4 <- plot_scatter(cvdata, "cv_mean", "slope")
g5 <- plot_scatter(cvdata, "cv_ratio", "slope")


(g1 + g2 + g4 + g5) + patchwork::plot_layout(guides = "collect")&theme(axis.title.y = element_blank()) #guides = "collect", design = layout

if (survey == "HBLL") ggsave("figs/explore-hbll-slopes.pdf", width = 9, height = 9)
if (survey == "SYN") ggsave("figs/explore-syn-slopes.pdf", width = 8, height = 6)

(g6 <- plot_scatter(cvdata, "`(Intercept)`", "slope")+ theme(legend.position = c(0.8,0.3), legend.title = element_blank()))
if (survey == "HBLL") ggsave("figs/hbll-slope-by-intercept-for-lm.pdf", width = 3, height = 3)
if (survey == "SYN") ggsave("figs/syn-slope-by-intercept-for-lm.pdf", width = 6, height = 6)


(g7 <- plot_scatter(cvdata, "cv_ratio", "cv_mean")+ theme(legend.position = c(0.8,0.3), legend.title = element_blank()))
if (survey == "HBLL") ggsave("figs/hbll-cv-mean-ratio.pdf", width = 5, height = 5)
if (survey == "SYN") ggsave("figs/syn-cv-mean-ratio.pdf", width = 5, height = 5)
