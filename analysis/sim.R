# remotes::install_github("pbs-assess/sdmTMB", ref="sim2")

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(mgcv)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)

source("analysis/functions.R")

#### TRY FOR HBLL ####
dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
grid$survey_abbrev <- "HBLL OUT N"
ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
d <- filter(dat_to_fit, species_common_name == "north pacific spiny dogfish")

m0 <- fit_geo_model(d, pred_grid = grid, survey = "HBLL",
  family = "tweedie",
  # MPA_trend = TRUE,
  return_model = TRUE, silent = FALSE, do_fit = FALSE)

mesh <- make_mesh(m0$data, xy_cols = c("X", "Y"), cutoff = 15)
mesh$mesh$n
dat <- m0$data
m <- sdmTMB(response ~ 1 + s(year, k = 5), data = dat, family = tweedie(), mesh = mesh,
  silent = FALSE, time = "year", spatiotemporal = "IID")
m

s <- sdmTMB::sdmTMB_simulate(
  previous_fit = m,
  simulate_re = TRUE
)

p <- predict(m, newdata = NULL)
ggplot(p, aes(X, Y, colour = omega_s)) + geom_point() +
  scale_colour_gradient2()

ggplot(s, aes(X, Y, colour = omega_s)) + geom_point() +
  scale_colour_gradient2()

ggplot(dat, aes(X, Y, colour = response)) + geom_point() +
  scale_color_viridis_c(trans = "sqrt")

ggplot(s, aes(X, Y, colour = observed)) + geom_point() +
  scale_color_viridis_c(trans = "sqrt")

b1 <- tidy(m)
b2 <- tidy(m, "ran_pars")
b <- bind_rows(b1, b2)

dat$year_covariate <- dat$year - min(dat$year)
mpa_increase_per_year <- 0.1

s <- sdmTMB::sdmTMB_simulate(
  formula = ~ 1 + restricted * year_covariate,
  data = dat,
  mesh = mesh,
  family = tweedie(),
  time = "year",
  # rho = 2 * plogis(m$model$par[['ar1_phi']]) - 1, # TODO!?
  rho = 0,
  sigma_O = b$estimate[b$term == "sigma_O"],
  sigma_E = b$estimate[b$term == "sigma_E"],
  phi = b$estimate[b$term == "phi"],
  tweedie_p = b$estimate[b$term == "tweedie_p"],
  range = b$estimate[b$term == "range"],
  simulate_re = TRUE,
  # (Intercept), restrictedTRUE, year_covariate, restrictedTRUE:year_covariate
  B = c(b$estimate[b$term == "(Intercept)"], 0, 0, mpa_increase_per_year)
)
s <- rename(s, restricted = restrictedTRUE)

ggplot(s, aes(X, Y, colour = omega_s)) + geom_point() +
  scale_colour_gradient2()

ggplot(s, aes(X, Y, colour = epsilon_st)) + geom_point() +
  scale_colour_gradient2() +
  facet_wrap(vars(year))

ggplot(s, aes(X, Y, colour = mu)) + geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~year) +
  geom_point(data = filter(s, restricted == 1L),
    pch = 21, colour = "black", fill = NA)

s %>%
  group_by(restricted, year) %>%
  summarise(m = mean(eta)) %>%
  group_by(restricted) %>%
  mutate(m = m / m[1]) %>%
  ggplot(aes(year, m, colour = as.factor(restricted))) +
  geom_line()

# s$year_covariate <- s$year_covariate - mean(s$year_covariate)
m2 <- sdmTMB(
  observed ~ 1 + restricted * year_covariate,
  data = s,
  mesh = mesh,
  time = "year",
  # spatiotemporal = "AR1",
  spatiotemporal = "IID",
  family = tweedie(),
  silent = FALSE
)

pred_grid <- grid
pred_grid <- filter(pred_grid, year %in% unique(dat$year))
pred_grid$year_covariate <- pred_grid$year - min(pred_grid$year)
p <- predict(m2, newdata = pred_grid, sims = 200L)

p_mpa <- p[pred_grid$restricted, ]
attr(p_mpa, "time") <- "year"

p_out <- p[!pred_grid$restricted, ]
attr(p_out, "time") <- "year"

ind_mpa <- get_index_sims(p_mpa)
ind_out <- get_index_sims(p_out)

ind <- ind_mpa %>% mutate(type = "Restricted") %>%
  bind_rows(mutate(ind_out, type = "Not restricted"))

# How to assess improvement?

# 1. visual on geostat index?
ind %>%
  group_by(type) %>%
  mutate(lwr = lwr / est[1], upr = upr / est[1], est = est / est[1]) %>%
  ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_ribbon(alpha = 0.1, colour = NA) + geom_line() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")

# 2. fancy geostat BACI?
tidy(m2)

# 3. naive GLM BACI on geostat output?
# glm(est ~ year, data = ind_mpa, family = Gamma(link = "log"))
# glm(est ~ year, data = ind_out, family = Gamma(link = "log"))
a <- glm(est ~ type*year, data = ind, family = Gamma(link = "log"))
broom::tidy(a)
# above ignores autocorrelation

# 4. naive GLM BACI on raw observations?
b <- mgcv::gam(observed ~ restricted * year_covariate, data = s, family = tw())
# above ignores spatial correlation and temporal autocorrelation and sampling location
summary(b)

# 5. plot of mean by area?
s %>%
  group_by(restricted, year) %>%
  summarise(m = mean(eta)) %>%
  group_by(restricted) %>%
  mutate(m = m / m[1]) %>%
  ggplot(aes(year, m, colour = as.factor(restricted))) +
  geom_line()
