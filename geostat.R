library(dplyr)
library(ggplot2)
library(sf)
library(future)
plan(multisession)
options(future.rng.onMisuse = "ignore")
library(sdmTMB)

source("load-data.R")
source("functions.R")

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
  m <- try({
    sdmTMB(density ~ 0 + as.factor(year),
      data = trawl_dat,
      family = tweedie(),
      time = "year",
      spde = mesh
    )
  })
  if (class(m) == "try-error") {
    return(NULL)
  }
  if (max(m$gradients) > 0.01) {
    m <- try({run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 0L)})
  }
  if (class(m) == "try-error") {
    return(NULL)
  }
  .grid <- gfplot::synoptic_grid %>%
    dplyr::filter(survey == unique(trawl_dat$survey_abbrev)) %>%
    expand_prediction_grid(years = unique(trawl_dat$year))
  set.seed(1)
  pred <- try({
    predict(m, newdata = .grid, xy_cols = c("X", "Y"), sims = 200L)
  })
  if (class(pred) == "try-error") {
    return(NULL)
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
      geo_wrapper() %>%
      mutate(type = "Restricted")
    bind_rows(out, out_restr)
  }, .progress = TRUE)

saveRDS(index_syn, file = "index-syn-geo.rds")
