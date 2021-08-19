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

dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")

geo_wrapper <- function(trawl_dat) {
  utm_zone9 <- 3156
  coords <- trawl_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  coords$X <- coords$X / 1000
  coords$Y <- coords$Y / 1000
  trawl_dat <- dplyr::bind_cols(trawl_dat, coords)
  trawl_dat$density <- trawl_dat$density_kgpm2 * 1000
  mesh <- make_mesh(trawl_dat, c("X", "Y"), cutoff = 15)
  # mesh$mesh$n

  null_df <- data.frame(
    species_common_name = trawl_dat$species_common_name[1],
    survey_abbrev = trawl_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )
  m <- try({
    sdmTMB(density ~ 0 + as.factor(year),
      data = trawl_dat,
      family = tweedie(),
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
  .grid <- gfplot::synoptic_grid %>%
    dplyr::filter(survey == unique(trawl_dat$survey_abbrev)) %>%
    expand_prediction_grid(years = unique(trawl_dat$year))
  set.seed(1)
  pred <- try({
    predict(m, newdata = .grid, xy_cols = c("X", "Y"), sims = 250L)
  })
  if (class(pred)[[1]] == "try-error") {
    return(null_df)
  }
  ind <- get_index_sims(pred, area = rep(4, nrow(pred)))
  ind$species_common_name <- trawl_dat$species_common_name[1]
  ind$survey_abbrev <- trawl_dat$survey_abbrev[1]
  ind$max_gradient <- max(m$gradients)
  ind
}

index_syn_geo <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    cat(.x$species_common_name[1], .x$survey_abbrev[1], "\n")
    out <- .x %>%
      geo_wrapper() %>%
      mutate(type = "Status quo")
    out_restr <- .x %>%
      filter(!restricted) %>%
      geo_wrapper() %>%
      mutate(type = "Restricted")
    bind_rows(out, out_restr)
  }, .progress = TRUE)

saveRDS(index_syn_geo, file = "index-syn-geo3.rds")

plan(sequential)

index_syn_geo <- readRDS("index-syn-geo3.rds")

index <- filter(index_syn_geo, !is.na(est), !is.na(se))

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

saveRDS(y, file = "index-syn-geo-clean.rds")

mean(y$orig_se < 1)

index <- filter(y, orig_se < 1)

g <- ggplot(index, aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
  geom_ribbon(alpha = 0.2, colour = NA)

ggsave("figs/index-syn-geo-restricted.pdf", width = 8, height = 60, limitsize = FALSE)

x <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(se_ratio = se[type == "Restricted"] / se[type == "Status quo"])

x %>%
  ggplot(aes(se_ratio)) + facet_wrap(~survey_abbrev) + geom_histogram() +
  geom_vline(xintercept = 1, col = "red") +
  coord_cartesian(expand = FALSE)

ggsave("figs/index-geo-syn-se-ratios.pdf", width = 8, height = 4, limitsize = FALSE)

x %>% group_by(survey_abbrev) %>%
  summarise(mean_ratio = mean(se_ratio))

plan(sequential)

# ??
index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))

x <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(re = (est[type == "Restricted"] - est[type == "Status quo"]) /
      est[type == "Status quo"])

g <- ggplot(x, aes(year, re)) +
  geom_line() +
  facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2)

ggsave("figs/index-syn-geo-restricted-re.pdf", width = 8, height = 65, limitsize = FALSE)

ggplot(x, aes(year, re, group = species_common_name)) +
  geom_line(alpha = 0.5) +
  facet_grid(~survey_abbrev, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2) +
  coord_cartesian(ylim = c(-0.5, 0.5))

x %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(mare = median(abs(re))) %>%
  ungroup() %>%
  arrange(-mare) %>%
  top_n(20) %>%
  as.data.frame()

x %>% group_by(species_common_name, survey_abbrev) %>%
  summarise(mare = median(abs(re))) %>%
  ggplot(aes(mare)) + geom_histogram(binwidth = 0.05) + facet_grid(~survey_abbrev) +
  coord_cartesian(xlim = c(0, 2))
