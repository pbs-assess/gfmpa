# remotes::install_github("pbs-assess/sdmTMB", ref="62119e0")
# remotes::install_github("pbs-assess/sdmTMB", ref="sim2")

library(sdmTMB)
library(dplyr)
library(ggplot2)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

source("analysis/functions.R")

#### TRY FOR HBLL ####
dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
grid$survey_abbrev <- "HBLL OUT N"
ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
d <- filter(dat_to_fit, species_common_name == "north pacific spiny dogfish")

m <- fit_geo_model(d, pred_grid = grid, survey = "HBLL", family = "tweedie",
  # MPA_trend = TRUE,
  return_model = TRUE, silent = T)


#### TRY FOR SYN ####
dat_to_fit <- readRDS("data-generated/dat_to_fit.rds") %>% filter(survey_abbrev == "SYN QCS")
grid <- readRDS("data-generated/syn-grid-w-restr.rds") %>% filter(survey_abbrev == "SYN QCS")
trawl_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

# d <- filter(dat_to_fit, species_common_name == "north pacific spiny dogfish")
# d <- filter(dat_to_fit, species_common_name == "arrowtooth flounder")
# d <- filter(dat_to_fit, species_common_name == "petrale sole")
d <- filter(dat_to_fit, species_common_name == "redbanded rockfish")

m <- fit_geo_model(d, pred_grid = grid, survey = "SYN", family = "tweedie",
  # MPA_trend = TRUE,
  return_model = TRUE, silent = T)

m

r <- residuals(m)
qq <- qqnorm(r)
plot(qq)
abline(a=0, b=1)
fe <- tidy(m)
re <- tidy(m, effects = "ran_pars") %>% distinct()
effects <- bind_rows(fe, re)

## use model's random field
p <- predict(m, newdata = NULL)

## if as.factor(year)
# recovery <- seq(5, 7, length.out = ncol(m$tmb_data$X_ij))
# plot(unique(m$data$year), exp(recovery))

sim_dat <- sdmTMB_sim2(
  formula = m$formula,
  data = m$data,
  time = "year",
  spde = m$spde,
  family = tweedie(link = "log"),
  tweedie_p = effects$estimate[effects$term == "tweedie_p"],
  phi = effects$estimate[effects$term == "phi"],
  ## reuse modelled random fields
  omega_s = p$omega_s,
  epsilon_st = p$epsilon_st,
  ## if we wanted new sim random fields
  # range = effects$estimate[effects$term == "range"],
  # sigma_E = effects$estimate[effects$term == "sigma_E"],
  # sigma_O = effects$estimate[effects$term == "sigma_O"],
  # seed = 3542,
  B = fe$estimate
)

dat <- m$data
ggplot(dat, aes(X, Y, size = density)) + geom_point() +
  scale_size_continuous(limits= c(0, quantile(dat$density, 0.99)))+
  facet_wrap(vars(year))

ggplot(sim_dat, aes(X, Y, size = observed)) + geom_point() +
  scale_size_continuous(limits= c(0, quantile(dat$density, 0.99)))+
  facet_wrap(vars(year))

# median(dat$density)
# median(sim_dat$observed)
#
# mean(sim_dat$observed)
# mean(dat$density)
#
hist(log(dat$density + 1))
hist(log(sim_dat$observed + 1))

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
