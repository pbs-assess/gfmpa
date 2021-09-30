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

dat <- m$data
mesh <- make_mesh(m$data, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh$mesh, asp = 1)

# m$tmb_data$X_ij %>% head

loc <- unique(select(dat, X, Y, year))

s <- sdmTMB::sdmTMB_sim(mesh = mesh, x = loc$X, y = loc$Y,
  X = matrix(data = 1, nrow = length(loc$X) * length(unique(loc$year))),
  betas = mean(fe$estimate),
  range = re$estimate[re$term == "range"],
  sigma_O = re$estimate[re$term == "sigma_O"],
  sigma_E = re$estimate[re$term == "sigma_E"],
  family = tweedie(),
  phi = re$estimate[re$term == "phi"],
  thetaf = plogis(m$model$par[["thetaf"]]) + 1,
  time_steps = length(unique(loc$year))
)

yr_lu <- data.frame(year = sort(unique(dat$year)),
  time = seq_len(length(unique(loc$year)))
)
s <- left_join(s, yr_lu)
s <- rename(s, X = x, Y = y)
s2 <- left_join(loc, s)

nrow(s2)
nrow(loc)

ggplot(dat, aes(X, Y, size = density)) + geom_point() +
  facet_wrap(vars(year))

ggplot(s2, aes(X, Y, size = observed)) + geom_point() +
  facet_wrap(vars(year))

mean(s2$observed)
median(s2$observed)
mean(dat$density)
median(dat$density)

