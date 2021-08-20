library(dplyr)
library(ggplot2)
library(sf)
library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 8L) else plan(multisession)
options(future.rng.onMisuse = "ignore")
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())

# source("load-data.R")
source("functions.R")

dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

geo_wrapper <- function(surv_dat, pred_grid, shrink_survey = FALSE) {
  utm_zone9 <- 3156
  coords <- surv_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  coords$X <- coords$X / 1000
  coords$Y <- coords$Y / 1000
  surv_dat <- dplyr::bind_cols(surv_dat, coords)
  surv_dat$density <- surv_dat$density_ppkm2
  mesh <- make_mesh(surv_dat, c("X", "Y"), cutoff = 15)
  # mesh$mesh$n

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )
  m <- try({
    sdmTMB(density ~ 0 + as.factor(year),
      data = surv_dat,
      family = sdmTMB::tweedie(),
      time = "year",
      spde = mesh
    )
  })
  # file_name <- paste0("data-generated/",
  #   gsub("\\/", "-", null_df$species_common_name[1]), "-",
  #   null_df$species_common_name[1], ".rds"
  # )
  if (class(m)[[1]] == "try-error") {
    return(null_df)
  }
  if (max(m$gradients) > 0.01) {
    m <- try({
      run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 0L)
    })
  }
  # saveRDS(m, file = file_name)
  if (class(m)[[1]] == "try-error") {
    return(null_df)
  }

  set.seed(1)
  pred <- try({
    predict(m, newdata = pred_grid, xy_cols = c("X", "Y"), sims = 250L)
  })
  if (class(pred)[[1]] == "try-error") {
    return(null_df)
  }
  ind <- get_index_sims(pred, area = rep(4, nrow(pred)))
  ind$species_common_name <- surv_dat$species_common_name[1]
  ind$survey_abbrev <- surv_dat$survey_abbrev[1]
  ind$max_gradient <- max(m$gradients)
  ind
}

hbll_grid <- gfplot::hbll_n_grid$grid
utm_zone9 <- 3156
coords <- hbll_grid %>%
  sf::st_as_sf(crs = 4326, coords = c("X", "Y")) %>%
  sf::st_transform(utm_zone9)

shrink_a_survey <- function(grid_dat, restriction_dat, plot = FALSE) {
  orig <- grid_dat
  grid_dat <- grid_dat %>% st_transform(sf::st_crs(restriction_dat))
  intersected <- sf::st_intersects(grid_dat, restriction_dat)
  remain <- which(lengths(intersected) == 0)
  lost <- which(lengths(intersected) > 0)
  if (plot) {
    .d <- as.data.frame(st_coordinates(orig))
    remain_df <- .d[remain,]
    lost_df <- .d[lost,]
    plot(remain_df$X, remain_df$Y, col = "black")
    points(lost_df$X, lost_df$Y, col = "red")
  }
  grid_dat$restricted <- lengths(intersected) > 0
  grid_dat <- as.data.frame(grid_dat) %>% select(-geometry)
  grid_dat$longitude <- orig$longitude
  grid_dat$latitude <- orig$latitude
  grid_dat %>% as_tibble()
}

coords_restr <- shrink_a_survey(coords, ll_removed, plot = FALSE)

coords <- coords %>%
  sf::st_coordinates() %>% as.data.frame()
coords$X <- coords$X / 1000
coords$Y <- coords$Y / 1000
coords$restricted <- coords_restr$restricted

hbll_grid <- expand_prediction_grid(coords, years = unique(dat_to_fit$year)) %>% as_tibble()

index_hbll_orig <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = hbll_grid) %>%
      mutate(type = "Status quo")
  }, .progress = TRUE)

index_hbll_restr <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = hbll_grid) %>%
      mutate(type = "Restricted")
  }, .progress = TRUE)

index_hbll_shrunk <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = filter(hbll_grid, restricted)) %>%
      mutate(type = "Restricted")
  }, .progress = TRUE)

index_hbll <- bind_rows(index_hbll_orig, index_hbll_restr)

saveRDS(index_hbll, file = "index-hbll-geo.rds")

plan(sequential)

index_hbll <- readRDS("index-hbll-geo.rds")

index <- filter(index_hbll, !is.na(est), !is.na(se))

max(index$max_gradient)

# remove ones with huge CVs originally; wouldn't use

# how many didn't fit in orig?

.u1 <- distinct(index, survey_abbrev, species_common_name, type) %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(orig = "Status quo" %in% type, restr = "Restricted" %in% type)

# didn't converge?
sum(.u1$orig)
sum(.u1$restr)

.u1 <- filter(.u1, orig & restr)

index <- left_join(.u1, index)

y <- index %>% group_by(survey_abbrev, species_common_name) %>%
  mutate(orig_se = max(se[type == "Status quo"]))

saveRDS(y, file = "index-hbll-geo-clean.rds")

# mean(y$orig_se < 1)
#
# index <- filter(y, orig_se < 1)

g <- ggplot(index, aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_ribbon(alpha = 0.2, colour = NA)

ggsave("figs/index-hbll-geo-restricted.pdf", width = 18, height = 18, limitsize = FALSE)

x <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(se_ratio = se[type == "Restricted"] / se[type == "Status quo"])

x %>%
  ggplot(aes(se_ratio)) + facet_wrap(~survey_abbrev) + geom_histogram() +
  geom_vline(xintercept = 1, col = "red") +
  coord_cartesian(expand = FALSE)

ggsave("figs/index-geo-hbll-se-ratios.pdf", width = 8, height = 4, limitsize = FALSE)

x %>% group_by(survey_abbrev) %>%
  summarise(mean_ratio = mean(se_ratio))

plan(sequential)

# # ??
# index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))

x <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(re = (est[type == "Restricted"] - est[type == "Status quo"]) /
      est[type == "Status quo"])

g <- ggplot(x, aes(year, re)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_hline(yintercept = 0, lty = 2)

ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 18, height = 18, limitsize = FALSE)

ggplot(x, aes(year, re, group = species_common_name)) +
  geom_line(alpha = 0.5) +
  facet_grid(~survey_abbrev, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = c(-0.5, 0.5))

x %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(mare = median(abs(re))) %>%
  ungroup() %>%
  arrange(-mare) %>%
  top_n(25) %>%
  as.data.frame()

x %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(mare = median(abs(re))) %>%
  ggplot(aes(mare)) + geom_histogram(binwidth = 0.05) + facet_grid(~survey_abbrev) +
  coord_cartesian(xlim = c(0, 2))
