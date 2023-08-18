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

shrink_domain <- function(grid_dat, restriction_dat) {
  grid_dat$SELECTION_I <- NULL
  grid_dat <- grid_dat |> sf::st_transform(sf::st_crs(restriction_dat))
  intersected <- sf::st_intersects(grid_dat, restriction_dat)
  grid_dat$restricted <- lengths(intersected) > 0
  grid_dat
}

trawl_removed <- readRDS("data-generated/Cat1_2_GH_Apr2023.rds") |>
  sf::st_transform(crs = 3156)
ll_removed <- readRDS("data-generated/LL_Cat1_2_GH_Apr2023.rds") |>
  sf::st_transform(crs = 3156)

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
    c(xmin = -133.5, ymin = 51, xmax = -127.5, ymax = 55))))
bc_coast <- st_transform(bc_coast, crs = st_crs(trawl_removed))

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

d1 <- parse_survey("QCS_Active_Blocks", trawl_removed)
d2 <- parse_survey("HS_Active_Blocks", trawl_removed)
d3 <- parse_survey("WCHG_Active_Blocks", trawl_removed)

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

d4 <- parse_survey("HBLLOut_Active_Blocks", ll_removed)

d1$survey_abbrev <- "SYN QCS"
d2$survey_abbrev <- "SYN HS"
d3$survey_abbrev <- "SYN WCHG"
d4$survey_abbrev <- "HBLL OUT N"
nsb_areas <- bind_rows(list(d1, d2, d3, d4))
saveRDS(nsb_areas, "data-generated/stratum-areas.rds")

plot_shrunk_grid <- function(f, restr_dat, return_data = FALSE) {
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

  g <- d |>
    ggplot(aes(fill = as.factor(GROUPING_CO))) +
    geom_sf(colour = NA, linewidth = 0) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    geom_sf(data = filter(d, restricted), colour = "#00000050", linewidth = 0.03) +
    theme_light() +
    geom_sf(data = bc_coast, inherit.aes = FALSE) +
    labs(fill = "Stratum ID") +
  coord_sf(
    xlim = c(164985.2 + 20000, 274985.2 + 300000),
    ylim = c(5851334.2 - 210000, 6047334.2)
  )
  if (return_data) d else g
}

g1 <- plot_shrunk_grid("WCHG_Active_Blocks", trawl_removed) + ggtitle("SYN WCHG")
g2 <- plot_shrunk_grid("HBLLOut_Active_Blocks", ll_removed) + ggtitle("HBLL OUT N")
g3 <- plot_shrunk_grid("HS_Active_Blocks", trawl_removed) + ggtitle("SYN HS")
g4 <- plot_shrunk_grid("QCS_Active_Blocks", trawl_removed) + ggtitle("SYN QCS")

g <- cowplot::plot_grid(g1, g2, g3, g4)
ggsave("figs/grid-strata-restricted.pdf", width = 9, height = 6.7, plot = g)
ggsave("figs/grid-strata-restricted.png", width = 9, height = 6.7, plot = g)

shrink_survey_grid <- function(f, restr_dat, label = "") {
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
  x <- shrink_domain(d, restr_dat)

  out <- x |> sf::st_centroid() |>
    sf::st_coordinates() |>
    as.data.frame() |>
    select(X, Y) |>
    mutate(X = X/1000, Y=Y/1000) |>
    mutate(restricted = x$restricted, survey_abbrev = label,
      grouping_code = x$GROUPING_CO, block = x$BLOCK_DESIG)
  row.names(out) <- NULL
  out
}

d1 <- shrink_survey_grid("WCHG_Active_Blocks", trawl_removed, label='SYN WCHG')
d2 <- shrink_survey_grid("HS_Active_Blocks", trawl_removed, label='SYN HS')
d3 <- shrink_survey_grid("QCS_Active_Blocks", trawl_removed, label='SYN QCS')
d4 <- shrink_survey_grid("HBLLOut_Active_Blocks", ll_removed, label='HBLL OUT N')

d <- bind_rows(list(d1,d2,d3,d4))
glimpse(d)

ggplot(d, aes(X,Y,color=restricted)) + geom_point()
ggplot(d, aes(X,Y,color=as.factor(grouping_code))) + geom_point()
ggplot(d, aes(X,Y,color=block)) + geom_point()

saveRDS(d, "data-generated/grids-strata-restricted.rds")
