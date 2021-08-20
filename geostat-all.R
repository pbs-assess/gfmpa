library(dplyr)
library(ggplot2)
library(sf)
library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 8L) else plan(multisession, workers = 8L)
options(future.rng.onMisuse = "ignore")
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())

# source("load-data.R")
source("functions.R")

# survey <- "HBLL"
survey <- "SYN"

if (survey == "HBLL") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
  grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
}
if (survey == "SYN") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
  grid <- readRDS("data-generated/syn-grid-w-restr.rds")
}

ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
trawl_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

index_orig <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      fit_geo_model(pred_grid = grid, survey = survey) %>%
      mutate(type = "Status quo")
  }, .progress = TRUE)

index_restr <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      fit_geo_model(pred_grid = grid, survey = survey) %>%
      mutate(type = "Restricted")
  }, .progress = TRUE)

index_shrunk <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
  # purrr::map_dfr(function(.x) {
    out <- .x %>%
      fit_geo_model(pred_grid = filter(grid, !restricted), survey = survey) %>%
      mutate(type = "Restricted and shrunk")
  }, .progress = TRUE)
  # })

index_all <- bind_rows(index_orig, index_restr, index_shrunk)

if (survey == "HBLL") saveRDS(index_all, file = "index-hbll-geo.rds")
if (survey == "SYN") saveRDS(index_all, file = "index-syn-geo.rds")

plan(sequential) # don't crash!

if (survey == "HBLL") index <- readRDS("index-hbll-geo.rds")
if (survey == "SYN") index <- readRDS("index-syn-geo.rds")

index <- filter(index, !is.na(est), !is.na(se)) # didn't fit
index <- index %>% mutate(cv = sqrt(exp(se^2) - 1))
stopifnot(max(index$max_gradient) < 0.01)

# how many didn't fit in orig?
# .u1 <- distinct(index, survey_abbrev, species_common_name, type) %>%
#   group_by(species_common_name, survey_abbrev) %>%
#   summarise(orig = "Status quo" %in% type, restr = grepl("Restricted", type))
# # didn't converge?
# sum(.u1$orig)
# sum(.u1$restr)
# .u1 <- filter(.u1, orig & restr)
# index <- left_join(.u1, index)
y <- index %>%
  group_by(survey_abbrev, species_common_name) %>%
  mutate(orig_cv = max(cv[type == "Status quo"]))

if (survey == "HBLL") saveRDS(y, file = "index-hbll-geo-clean.rds")
if (survey == "SYN") saveRDS(y, file = "index-syn-geo-clean.rds")

mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
index <- filter(y, orig_cv < 1)

g <- ggplot(index, aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_ribbon(alpha = 0.2, colour = NA)

if (survey == "HBLL") ggsave("figs/index-hbll-geo-restricted.pdf", width = 18, height = 18)
if (survey == "SYN") ggsave("figs/index-syn-geo-restricted.pdf", width = 18, height = 18, limitsize = FALSE)

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
  ggplot(aes(`CV ratio`)) +
  facet_wrap(~survey_abbrev) +
  geom_histogram() +
  geom_vline(xintercept = 1, col = "red") +
  facet_wrap(vars(`Restriction type`)) +
  coord_cartesian(expand = FALSE)

if (survey == "HBLL") ggsave("figs/index-geo-hbll-cv-ratios.pdf", width = 8, height = 4)
if (survey == "SYN") ggsave("figs/index-geo-syn-cv-ratios.pdf", width = 8, height = 4)

x_long %>%
  group_by(`Restriction type`) %>%
  summarise(mean_ratio = mean(`CV ratio`))

lu <- tibble(
  "Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)
x_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
  left_join(lu) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  coord_flip() +
  xlab("") +
  ylab("Ratio of index CV\n(restricted/original)") +
  labs(colour = "Restricted survey domain") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
ggsave("figs/index-geo-hbll-cv-ratio-dotplot.pdf", width = 7, height = 7)
ggsave("figs/index-geo-hbll-cv-ratio-dotplot.png", width = 7, height = 7)

# bad year in SYN
if (survey == "SYN") {
  index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))
}

x <- index %>%
  group_by(species_common_name, survey_abbrev, type) %>%
  mutate(est = est / exp(mean(log(est)))) %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(
    re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
    re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
  )

x_long <- x %>%
  tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re")

g <- ggplot(x_long, aes(year, re, colour = `Restriction type`)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_hline(yintercept = 0, lty = 2)

if (survey == "HBLL") ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 10, height = 9)
if (survey == "SYN") ggsave("figs/index-syn-geo-restricted-re.pdf", width = 10, height = 9)

lu <- tibble(
  "Restriction type" = c("re_restr", "re_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain")
)

x_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
  left_join(lu) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  geom_point(position = position_dodge(width = 0)) +
  coord_flip() +
  xlab("") +
  ylab("Median absolute relative error (MARE)\n(restricted compared to original)") +
  labs(colour = "Restricted survey domain") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))

if (survey == "HBLL") ggsave("figs/index-geo-hbll-mare-dotplot.pdf", width = 7, height = 7)
if (survey == "SYN") ggsave("figs/index-geo-syn-mare-dotplot.png", width = 7, height = 7)
