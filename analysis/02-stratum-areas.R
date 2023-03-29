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

get_stratum_areas_restr <- function(dat) {
  codes <- sort(unique(dat$GROUPING_CO))
  purrr::map_dfr(codes, function(i) {
    a <- dat[dat$GROUPING_CO == i & !dat$restricted,] |> st_area() |> sum()
    tibble(grouping_code = i, area = as.numeric(a) * 1e-6)
  })
}

ensure_multipolygons <- function(X) {
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  sf::st_write(X, tmp1)
  gdalUtilities::ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- sf::st_read(tmp2)
  sf::st_sf(sf::st_drop_geometry(X), geom = sf::st_geometry(Y))
}

shrink_domain <- function(grid_dat, restriction_dat) {
  grid_dat$SELECTION_I <- NULL
  grid_dat <- grid_dat |> sf::st_transform(sf::st_crs(restriction_dat))
  intersected <- sf::st_intersects(grid_dat, restriction_dat)
  grid_dat$restricted <- lengths(intersected) > 0
  grid_dat
}

nsb <- sf::read_sf("data-raw/MPATT_P2_Nov25_limited_attributes.gdb/", type = 7)
restriction_dat <- dplyr::filter(nsb, Category_Detailed %in% c("Category 1", "Category 2"))
restriction_dat <- ensure_multipolygons(restriction_dat)

parse_survey <- function(f, restr_dat) {
  d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/",
    layer = f)
  d <- shrink_domain(d, restr_dat)
  filter(d) |>
    ggplot(aes(colour = restricted)) + geom_sf()
  o1 <- get_stratum_areas(d)
  o2 <- get_stratum_areas_restr(d) |> rename(area_nsb = area)
  left_join(o1, o2)
}

d1 <- parse_survey("QCS_Active_Blocks", restriction_dat)
d2 <- parse_survey("HS_Active_Blocks", restriction_dat)
d3 <- parse_survey("WCHG_Active_Blocks", restriction_dat)
d4 <- parse_survey("HBLLOut_Active_Blocks", restriction_dat)

d1$survey_abbrev <- "SYN QCS"
d2$survey_abbrev <- "SYN HS"
d3$survey_abbrev <- "SYN WCHG"
d4$survey_abbrev <- "HBLL OUT N"
nsb_areas <- bind_rows(list(d1, d2, d3, d4))
saveRDS(nsb_areas, "data-generated/stratum-areas.rds")

hbll <- readRDS("data-raw/hbll-survey-data.rds")
unique(hbll$grouping_code)
# different!?

d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/", layer = "HBLLOut_Active_Blocks")
d <- d[d$SURVEY_SERI == 22, ]

# assign these grouping codes:
bl <- readRDS("data-raw/22blocks.rds")
hbll <- left_join(hbll, bl)
bl2 <- select(d, BLOCK_DESIG, GROUPING_CO) |>
  as.data.frame() |> select(1:2) |>
  rename(
    block_designation = BLOCK_DESIG, grouping_code = GROUPING_CO)
hbll$grouping_code <- NULL # wrong ones!
hbll <- left_join(hbll, bl2) # fixed
sum(is.na(hbll$grouping_code)) / nrow(hbll) * 100

stopifnot(sum(is.na(hbll$grouping_code)) == 0L)
