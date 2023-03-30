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

test <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/",
  layer = "HBLLOut_Active_Blocks")
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
    c(xmin = -133.5, ymin = 51, xmax = -127.5, ymax = 55))))
bc_coast <- st_transform(bc_coast, crs = st_crs(test))

parse_survey <- function(f, restr_dat) {
  d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/",
    layer = f)
  d$SELECTION_I <- NULL
  if (f == "HBLLOut_Active_Blocks") {
    d <- filter(d, SURVEY_SERI %in% 22)
    d <- rename(d, grouping_code2 = GROUPING_CO) |>
      left_join(a) |>  # above; global
      rename(GROUPING_CO = grouping_code1) |>
      select(-grouping_code2)
  }
  d <- shrink_domain(d, restr_dat)
  o1 <- get_stratum_areas(d)
  o2 <- get_stratum_areas_restr(d) |> rename(area_nsb = area)
  left_join(o1, o2)
}

d1 <- parse_survey("QCS_Active_Blocks", restriction_dat)
d2 <- parse_survey("HS_Active_Blocks", restriction_dat)
d3 <- parse_survey("WCHG_Active_Blocks", restriction_dat)

# HBLL needs to swap straum codes first:
# d <- gfdata::run_sql("GFBioSQL", "SELECT * FROM GROUPING")
# saveRDS(grouping22, "data-raw/grouping22.rds")
grouping22 <- readRDS("data-raw/grouping22.rds")
grouping22 <- filter(grouping22, SURVEY_SERIES_ID %in% c(22, 23))
# 23!?
a1 <- filter(grouping22, grepl("HBLL OUT North", GROUPING_DESC)) |>
  select(grouping_code1 = GROUPING_CODE, GROUPING_DEPTH_ID)
a2 <- filter(grouping22, !grepl("HBLL", GROUPING_DESC)) |>
  select(grouping_code2 = GROUPING_CODE, GROUPING_DEPTH_ID, SURVEY_SERIES_ID)
a <- left_join(a1, a2)

d4 <- parse_survey("HBLLOut_Active_Blocks", restriction_dat)

d1$survey_abbrev <- "SYN QCS"
d2$survey_abbrev <- "SYN HS"
d3$survey_abbrev <- "SYN WCHG"
d4$survey_abbrev <- "HBLL OUT N"
nsb_areas <- bind_rows(list(d1, d2, d3, d4))
saveRDS(nsb_areas, "data-generated/stratum-areas.rds")

plot_shrunk_grid <- function(f, restr_dat) {
  d <- sf::read_sf("data-raw/Synoptic_Surveys_and_HBLLOutside_Active_Blocks/",
    layer = f)
  d$SELECTION_I <- NULL
  if (f == "HBLLOut_Active_Blocks") {
    d <- filter(d, SURVEY_SERI %in% 22)
    d <- rename(d, grouping_code2 = GROUPING_CO) |>
      left_join(a) |>  # above; global
      rename(GROUPING_CO = grouping_code1) |>
      select(-grouping_code2)
  }
  pal <- c(as.character(colorBlindness::availableColors())[-1], c("grey60"))
  d <- shrink_domain(d, restr_dat)
  filter(d) |>
    ggplot(aes(fill = as.factor(GROUPING_CO))) + geom_sf(colour = NA, linewidth = 0) +
    # scale_colour_manual(values = c("white", "grey20")) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    # ggthemes::scale_colour_colorblind() +
    # ggthemes::scale_fill_colorblind() +
    geom_sf(data = filter(d, restricted), colour = "black", fill = "#00000010", linewidth = 0.05) +
    theme_light() +
    geom_sf(data = bc_coast, inherit.aes = FALSE) +
    labs(fill = "Stratum ID") +
    coord_sf(
      xlim = c(505961.3 - 20000, 842564.7 + 20000),
      ylim = c(758369.3 - 90000, 1092017.2 - 20000)
    )
}

g1 <- plot_shrunk_grid("WCHG_Active_Blocks", restriction_dat)
g2 <- plot_shrunk_grid("HBLLOut_Active_Blocks", restriction_dat)
g3 <- plot_shrunk_grid("HS_Active_Blocks", restriction_dat)
g4 <- plot_shrunk_grid("QCS_Active_Blocks", restriction_dat)
g <- cowplot::plot_grid(g1, g2, g3, g4)
ggsave("figs/grid-strata-restricted.pdf", width = 9, height = 7)
