library(ggplot2)
library(sf)
library(dplyr)

get_stratum_areas <- function(dat) {
  codes <- sort(unique(dat$GROUPING_CO))
  purrr::map_dfr(codes, function(i) {
    a <- dat[dat$GROUPING_CO == i,] |> st_area() |> sum()
    tibble(grouping_code = i, area = as.numeric(a) * 1e-6)
  })
}

d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/", layer = "WCHG_Active_Blocks")
get_stratum_areas(d)

syn <- readRDS("data-raw/syn-survey-data.rds")
x <- filter(syn, species_science_name == "bathyraja abyssicola") |>
  select(survey_abbrev, grouping_code, area_km2) |>
  distinct()
x

x <- filter(syn, species_science_name == "bathyraja abyssicola") |>
  select(survey_abbrev, grouping_code, area_km2) |>
  distinct() |>
  group_by(survey_abbrev) |>
  summarise(a = sum(area_km2))
x

d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/", layer = "QCS_Active_Blocks")
a <- get_stratum_areas(d)
a
sum(a$area)

d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/", layer = "HS_Active_Blocks")
a <- get_stratum_areas(d)
a
sum(a$area)

hbll <- readRDS("data-raw/hbll-survey-data.rds")
unique(hbll$grouping_code)
glimpse(hbll)

# centers_sf <- purrr::map2(x$longitude, x$latitude, ~st_point(c(.x, .y))) %>%
#   st_sfc(crs = st_crs(d)) %>%
#   st_sf(x[,-(1:2)], .)

x <- hbll[1,]
pt <- st_point(c(x$longitude, x$latitude))
st_crs(pt)
st_crs(d)
pt <- st_sfc(pt, crs = st_crs(d))
st_crs(pt)

plot(pt);axis(1);axis(2)
plot(d, max.plot = 1);axis(1);axis(2)
points(pt)

sf::st_intersects(pt, d)

d[as.numeric(sf::st_within(centers_sf, d)),]


d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/", layer = "HBLLOut_Active_Blocks")
d <- d[d$SURVEY_SERI == 22, ]
a <- get_stratum_areas(d)
a
sum(a$area)

hbll <- readRDS("data-raw/hbll-survey-data.rds")

x <- filter(hbll, species_science_name == "bathyraja abyssicola", year == 2006) |>
  select(area_km2) |>
  distinct()
sum(x$area_km2)

# assign these grouping codes:
hbll <- readRDS("data-raw/hbll-survey-data.rds")
bl <- readRDS("~/Downloads/22blocks.rds")
hbll <- left_join(hbll, bl)

bl2 <- select(d, BLOCK_DESIG, GROUPING_CO) |>
  as.data.frame() |> select(1:2) |>
  rename(
    block_designation = BLOCK_DESIG, grouping_code = GROUPING_CO)
# hbll <- left_join(hbll, bl2)
sum(is.na(bl2$grouping_code))

locs <- select(hbll, longitude, latitude) |> distinct()

d2 <- st_transform(d, crs = 32609)
plot(d2, max.plot = 1);axis(1);axis(2)

.crs <- st_crs(d2)

plot(d, max.plot = 1);axis(1);axis(2)

pt <- sf::st_point(c(-130.8217, 53.82167))
pt <- sf::st_sfc(pt, crs = 4326)
pt <- sf::st_transform(pt, crs = st_crs(d))
st_within(pt, d)
plot(pt);axis(1);axis(2)
d <- select(d, -7)


d_points <- data.frame(
  long = -130.595, lat  = 54.15) %>%
  # long = -132.4783 , lat  = 53.01) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = st_crs(d))

ggplot(d) + geom_sf(aes(colour = as.factor(GROUPING_CO))) + geom_sf(data = d_points, pch = 21, colour = "red")

st_within(d_points, d)


gg <- gfplot::hbll_grid$grid
head(gg)

library(future)
plan(multisession)
locs_found <- furrr::future_map_dbl(seq_len(nrow(locs)), function(.x) {
# locs_found <- purrr::map_dbl(1:600, function(.x) {
  cat(.x, "\n")
  x <- locs[.x, ]
  pt <- sf::st_point(c(x$longitude, x$latitude))
  pt <- sf::st_sfc(pt, crs = 4326)
  pt <- st_transform(pt, crs = st_crs(d))
  # a <- st_intersects(d, pt)
  # if (sum(lengths(a) > 0) >= 1) {
    # out <- d[lengths(a) > 0, "GROUPING_CO"][[1]]
  # } else {
    # out <- NA
  # }
  # out
  as.numeric(d[as.numeric(sf::st_within(pt, d)), "BLOCK_DESIG"][[1]])
})
plan(sequential)

which(is.na(locs_found))

d_points <- locs[which(is.na(locs_found)), ] |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(d))

library(ggplot2)
ggplot(d) + geom_sf() + geom_sf(data = d_points, pch = 21, colour = "red")

locs$grouping_code2 <- locs_found
hbll2 <- left_join(hbll, locs)

sum(is.na(hbll2$grouping_code2)) / sum(!is.na(hbll2$grouping_code2)) * 100


glimpse(hbll2)












library(dplyr)
library(ggplot2)

# calculate design-based biomass estimate from output of get_survey_sets()
calc_bio <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>%
    group_by(year, survey_abbrev, area_km2, grouping_code) %>%
    summarise(density = mean(density_kgpm2 * 1e6), .groups = "drop_last") %>%
    group_by(year) %>%
    summarise(
      survey_abbrev = dat$survey_abbrev[1],
      species_common_name = dat$species_common_name[1],
      est = sum(density * area_km2), .groups = "drop")
}

surv_dat <- readRDS("data-generated/dat_to_fit.rds")
# surv_dat <- readRDS("data-generated/dat_to_fit_hbll.rds")

ind_stat_quo <- group_by(surv_dat, species_common_name, survey_abbrev) |>
  group_split() |>
  purrr::map_dfr(calc_bio)

ind_restr_expand <- dplyr::filter(surv_dat, !restricted) |>
  group_by(species_common_name, survey_abbrev) |>
  group_split() |>
  purrr::map_dfr(calc_bio)

