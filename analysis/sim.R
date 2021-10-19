library(sdmTMB)
library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)
source("analysis/functions.R")

dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
grid$survey_abbrev <- "HBLL OUT N"
ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

d <- filter(dat_to_fit, species_common_name == "north pacific spiny dogfish")
m <- fit_geo_model(d, pred_grid = grid, survey = "HBLL", family = "tweedie",
  return_model = TRUE, silent = FALSE)

fe <- tidy(m)
re <- tidy(m, effects = "ran_pars") %>% distinct()
effects <- bind_rows(fe, re)

####### Fix from here down:

recovery <- seq(5, 7, length.out = ncol(m$tmb_data$X_ij))
plot(unique(m$data$year), exp(recovery))

sim_dat <- sdmTMB_sim2(
  formula = m$formula,
  data = m$data,
  time = "year",
  spde = m$spde,
  family = tweedie(link = "log"),
  range = effects$estimate[effects$term == "range"],
  sigma_E = effects$estimate[effects$term == "sigma_E"],
  tweedie_p = effects$estimate[effects$term == "tweedie_p"],
  phi = effects$estimate[effects$term == "phi"],
  sigma_O = effects$estimate[effects$term == "sigma_O"],
  # seed = 3542,
  B = recovery)
  # B = fe$estimate
)

dat <- m$data
ggplot(dat, aes(X, Y, size = density)) + geom_point() +
  facet_wrap(vars(year))

ggplot(sim_dat, aes(X, Y, size = observed)) + geom_point() +
  facet_wrap(vars(year))

# median(dat$density)
# median(sim_dat$observed)
#
# mean(sim_dat$observed)
# mean(dat$density)
#
# hist(log(dat$density + 1))
# hist(log(sim_dat$observed + 1))

sim_dat$response <- sim_dat$observed
m_sim <- sdmTMB(
  m$formula,
  data = sim_dat,
  time = "year",
  spde = m$spde,
  family = tweedie(link = "log"),
  silent = FALSE
)

m_sim

# now predict only in an MPA and check index!
