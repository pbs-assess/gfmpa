# TODO:
# - include depth quadratic?
# - what to do with time? s()? AR1? linear?
# - which tests to do after?
# - what time span to check?
# - what ramps of 'recovery' to check?

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 7L) else plan(multisession, workers = 7L)
options(future.rng.onMisuse = "ignore")

library(mgcv)

theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

source("analysis/functions.R")

# Globals to set ------------------------------
survey <- "HBLL"
# survey <- "SYN"
# ---------------------------------------------

if (survey == "HBLL") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds") %>% rename(depth = depth_m)
  grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
  grid$survey_abbrev <- "HBLL OUT N"
  save_file <- paste0("data-generated/sim-mpa-hbll.rds")
}

if (survey == "SYN") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit.rds") %>% rename(depth = depth_m)
  grid <- readRDS("data-generated/syn-grid-w-restr.rds")
  save_file <- paste0("data-generated/sim-mpa-syn.rds")
}

# ## test run on just one spp
# d <- filter(dat_to_fit, species_common_name == "north pacific spiny dogfish")
# s1 <- sim_mpa_surv(surv_dat = d, grid = grid,
#                   survey = "HBLL",
#                   mpa_increase_per_year = log(1.05),
#                   spatial = "on",
#                   spatiotemporal = "IID"
#                   # spatiotemporal = "AR1" #worse
# )

s <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    # purrr::map_dfr(function(.x) {
    out <- .x %>%
      sim_mpa_surv(grid = grid, survey = survey,
                   mpa_increase_per_year = log(1.05),
                   ## these go through as ... for true model underlying sim
                   spatial = "on",
                   spatiotemporal = "IID"
                   )
  }, .progress = TRUE)
saveRDS(s, file = save_file)
plan(sequential)



s <- readRDS(file = save_file) %>%
  mutate(species = as.factor(stringr::str_to_title(species_common_name)))


# 1. visual on geostat index?
s %>% group_by(survey_abbrev, species_common_name, type) %>%
  # remove a couple extreme values in final year to make axes more reasonable
  filter(!(species_common_name %in% c("sandpaper skate","petrale sole") & year == max(s$year))) %>%
  mutate(lwr = lwr / est[1], upr = upr/est[1], est = est/est[1],
         mean_eta = mean_eta / mean_eta[1]) %>%
  ggplot(aes(year, est, ymin = lwr, ymax = ifelse(upr < est * 3, upr, est * 3),
             colour = type, fill = type)) +
  geom_ribbon(alpha = 0.1, colour = NA) + geom_line() +
  # add mean eta from simulated data
  geom_line(aes(year, mean_eta, colour = type), linetype = "dashed") +
  facet_wrap(~species_common_name, ncol = 4, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")

ggsave(paste0("sim-index-", survey, ".png"), width = 8, height = 8)


# 2. fancy geostat BACI?
slopes <- s %>% select(survey_abbrev, species, type, true_slope,
             m = m_slope, lwr = m_slope_lwr, upr = m_slope_upr) %>% distinct()

g1 <- slopes %>% ggplot(aes(group = type, colour = type)) +
  geom_pointrange(aes(m, species,  xmin = lwr, xmax = upr),
                  position = position_jitterdodge(
                    dodge.width =0.5, jitter.height = 0, jitter.width = 0
                    )) +
  geom_vline(aes(xintercept = true_slope, colour = type)) +
  scale_y_discrete(limits = rev(levels(unique(slopes$species)))) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  ggtitle("A. fancy geostat BACI")


# 3. naive GLM BACI on geostat output?
# ignores autocorrelation

a <- s %>% group_by(survey_abbrev, species_common_name) %>% group_split() %>%
  lapply(function(.x){
    try({
    out <- glm(est ~ type*year, family = Gamma(link = "log"), data = .x,  weights = NULL)
    df <- broom::tidy(out)
    df$survey_abbrev <- .x$survey_abbrev[1]
    df$species_common_name <- .x$species_common_name[1]
    df
    })
  }
  )
a <- do.call(rbind, a)

a <- a %>% filter(term %in% c("year", "typeRestricted:year")) %>%
  mutate(type = ifelse(term=="year", "Not restricted", "Restricted"),
         species = as.factor(stringr::str_to_title(species_common_name)),
         true_slope = ifelse(term=="year", 0, s$true_slope[1]),
         m = as.numeric(estimate),
         lwr = m - (as.numeric(std.error)*1.96),
         upr = m + (as.numeric(std.error)*1.96)
         )

g2 <- a %>% ggplot(aes(group = type, colour = type)) +
  geom_pointrange(aes(m, species,  xmin = lwr, xmax = upr),
                  position = position_jitterdodge(
                    dodge.width =0.5, jitter.height = 0, jitter.width = 0
                  )) +
  geom_vline(aes(xintercept = true_slope, colour = type)) +
  scale_y_discrete(limits = rev(levels(unique(slopes$species)))) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  ggtitle("B. naive GLM BACI") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())


# 3.1 weight by 1/SE?
a <- s %>% group_by(survey_abbrev, species_common_name) %>% group_split() %>%
  lapply(function(.x){
    try({
      out <- glm(est ~ type*year, family = Gamma(link = "log"), data = .x,  weights = 1/se)
      df <- broom::tidy(out)
      df$survey_abbrev <- .x$survey_abbrev[1]
      df$species_common_name <- .x$species_common_name[1]
      df
    })
  }
  )
# broom::tidy(a[[1]])
a <- do.call(rbind, a)

a <- a %>% filter(term %in% c("year", "typeRestricted:year")) %>%
  mutate(type = ifelse(term=="year", "Not restricted", "Restricted"),
         species = as.factor(stringr::str_to_title(species_common_name)),
         true_slope = ifelse(term=="year", 0, s$true_slope[1]),
         m = as.numeric(estimate),
         lwr = m - (as.numeric(std.error)*1.96),
         upr = m + (as.numeric(std.error)*1.96)
  )

g3 <- a %>% ggplot(aes(group = type, colour = type)) +
  geom_pointrange(aes(m, species,  xmin = lwr, xmax = upr),
                  position = position_jitterdodge(
                    dodge.width =0.5, jitter.height = 0, jitter.width = 0
                  )) +
  geom_vline(aes(xintercept = true_slope, colour = type)) +
  scale_y_discrete(limits = rev(levels(unique(slopes$species)))) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  ggtitle("C. weight by 1/SE") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())


# 3.2 via lm()?
a <- s %>% group_by(survey_abbrev, species_common_name) %>% group_split() %>%
  lapply(function(.x){
    try({
      out <- lm(log(est) ~ type*year, data = .x)
      df <- broom::tidy(out)
      df$survey_abbrev <- .x$survey_abbrev[1]
      df$species_common_name <- .x$species_common_name[1]
      df
    })
  }
  )
a <- do.call(rbind, a)

a <- a %>% filter(term %in% c("year", "typeRestricted:year")) %>%
  mutate(type = ifelse(term=="year", "Not restricted", "Restricted"),
         species = as.factor(stringr::str_to_title(species_common_name)),
         true_slope = ifelse(term=="year", 0, s$true_slope[1]),
         m = as.numeric(estimate),
         lwr = m - (as.numeric(std.error)*1.96),
         upr = m + (as.numeric(std.error)*1.96)
  )

g4 <- a %>% ggplot(aes(group = type, colour = type)) +
  geom_pointrange(aes(m, species,  xmin = lwr, xmax = upr),
                  position = position_jitterdodge(
                    dodge.width =0.5, jitter.height = 0, jitter.width = 0
                  )) +
  geom_vline(aes(xintercept = true_slope, colour = type)) +
  scale_y_discrete(limits = rev(levels(unique(slopes$species)))) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  ggtitle("D. lm") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())


# 3.3 via gls() with AR1 residuals?
# a <- gls(log(est) ~ type*year, data = s$index, correlation = corAR1(form = ~ 1 | type))
# summary(a)

a <- s %>% group_by(survey_abbrev, species_common_name) %>% group_split() %>%
  lapply(function(.x){
    try({
      out <- gls(log(est) ~ type*year, correlation = corAR1(form = ~ 1 | type), data = .x)
      df <- as.data.frame(summary(out)$coefficients)
      # browser()
      names(df) <- "estimate"
      term <- rownames(df)
      rownames(df) <- NULL
      df <- cbind(term,df)
      df$se <- summary(out)$tTable[,2]
      df$survey_abbrev <- .x$survey_abbrev[1]
      df$species_common_name <- .x$species_common_name[1]
      df
    })
  }
  )

a <- do.call(rbind, a)

a <- a %>% filter(term %in% c("year", "typeRestricted:year")) %>%
  mutate(type = ifelse(term=="year", "Not restricted", "Restricted"),
         species = as.factor(stringr::str_to_title(species_common_name)),
         true_slope = ifelse(term=="year", 0, s$true_slope[1]),
         m = as.numeric(estimate),
         lwr = m - (se*1.96),
         upr = m + (se*1.96)
  )

g5 <- a %>% ggplot(aes(group = type, colour = type)) +
  geom_pointrange(aes(m, species,  xmin = lwr, xmax = upr),
                  position = position_jitterdodge(
                    dodge.width =0.5, jitter.height = 0, jitter.width = 0
                  )) +
  geom_vline(aes(xintercept = true_slope, colour = type)) +
  scale_y_discrete(limits = rev(levels(unique(slopes$species)))) +
  coord_cartesian(xlim = c(-0.25, 0.5)) +
  ggtitle(" E. gls with AR1 residuals") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())


g1 + g2 + g3 + g4 + g5 + patchwork::plot_layout(ncol = 5, guides = 'collect')
ggsave(paste0("BACI-coefs-", survey, ".png"), width = 12, height = 8)
