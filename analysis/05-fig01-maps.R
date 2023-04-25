# mapping MPAs
library(dplyr)
library(ggplot2)
library(sf)
ggplot2::theme_set(ggsidekick::theme_sleek())

#### TRAWL ####

dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
trawl <- readRDS("data-generated/syn-grid-w-restr.rds") %>% filter(year == 2018) %>% select(-year)

df <- trawl %>%
  dplyr::mutate(X = X*1000, Y = Y*1000) %>%
  gfplot:::utm2ll(., utm_zone = 9)

coast <- gfplot:::load_coastline(
  range(df$X) + c(-1, 1),
  range(df$Y) + c(-1, 1),
  utm_zone = 9
)

dat <- dat_to_fit

# ggplot(trawl, aes(X, Y)) +
#   geom_tile(aes(fill = restricted), width=2, height=2) +
#   geom_polygon(
#     data = coast, aes(x = X, y = Y, group = PID),
#     fill = "grey87", col = "grey70", lwd = 0.2
#   )+
#   scale_fill_brewer(palette = "Set1", direction = -1) +
#   coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
#   gfplot::theme_pbs() + xlab("UTM") + ylab("UTM")

dat <- dat %>% filter(species_common_name == "arrowtooth flounder") %>%
  mutate(year_pair = case_when(
    year %in% c(2001, 2002) ~ "2001-2002",
    year %in% c(2003, 2004) ~ "2003-2004",
    year %in% c(2005, 2006) ~ "2005-2006",
    year %in% c(2007, 2008) ~ "2007-2008",
    year %in% c(2009, 2010) ~ "2009-2010",
    year %in% c(2011, 2012) ~ "2011-2012",
    year %in% c(2013, 2014) ~ "2013-2014",
    year %in% c(2015, 2016) ~ "2015-2016",
    year %in% c(2017, 2018) ~ "2017-2018",
    year %in% c(2019, 2020) ~ "2019-2020",  # no WCVI
    year %in% c(2021) ~ "2021"
  ),
    year_true = year
  )

ggplot(trawl) +
  geom_tile(aes(X, Y, fill = restricted), alpha = 0.45, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y, colour = restricted), size = 0.1) +
  scale_fill_brewer("Restricted", palette = "Set1", direction = -1) +
  scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.15,0.15)) +
  xlab("Easting (km)") + ylab("Northing (km)")
ggsave("figs/trawl-grid.png", width = 6, height = 6)

g <- ggplot(trawl) +
  geom_tile(aes(X, Y, fill = restricted), alpha = 0.2, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y, colour = restricted), size = 0.1, shape = 20, stroke = 1) +
  scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  facet_wrap(~year_pair) +
  scale_fill_brewer("Restricted", palette = "Set1", direction = -1) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  # xlab("Easting (km)") + ylab("Northing (km)")+
  gfplot::theme_pbs() + theme(legend.position=c(0.07,0.75),
    strip.text.x = element_text(size = 12),
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()
    )
ggsave("figs/trawl-grid-by-year.png", width = 12, height = 12)

#### HBLL ####
dat_to_fit_hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
hbll <- readRDS("data-generated/hbll-n-grid-w-restr.rds") %>% filter(year == 2019) %>% select(-year)

df_hbll <- hbll %>%
  dplyr::mutate(X = X*1000, Y = Y*1000) %>%
  gfplot:::utm2ll(., utm_zone = 9)

dat_hbll <- dat_to_fit_hbll

dat_hbll <- dat_hbll %>% filter(species_common_name == "arrowtooth flounder") %>%
  mutate(year_pair = case_when(
        year %in% c(2001, 2002) ~ "2001-2002",
        year %in% c(2003, 2004) ~ "2003-2004",
        year %in% c(2005, 2006) ~ "2005-2006",
        year %in% c(2007, 2008) ~ "2007-2008",
        year %in% c(2009, 2010) ~ "2009-2010",
        year %in% c(2011, 2012) ~ "2011-2012",
        year %in% c(2013, 2014) ~ "2013-2014",
        year %in% c(2015, 2016) ~ "2015-2016",
        year %in% c(2017, 2018) ~ "2017-2018",
        year %in% c(2019, 2020) ~ "2019-2020",  # no WCVI
        year %in% c(2021) ~ "2021"
      ),
      year_true = year
  )

ggplot(hbll) +
  geom_tile(aes(X, Y, fill = restricted), alpha = 0.5, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat_hbll, aes(X, Y, colour = restricted), size = 0.15) +
  scale_fill_brewer("Restricted", palette = "Set1", direction = -1) +
  scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  coord_fixed(xlim = c(212, 533), ylim = c(5733, 6057)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.15,0.15)) +
  xlab("Easting (km)") + ylab("Northing (km)")
ggsave("figs/hbll-grid.png", width = 6, height = 6)

g <- ggplot(hbll) +
  geom_tile(aes(X, Y, fill = restricted), alpha = 0.2, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat_hbll, aes(X, Y, colour = restricted), size = 0.2, shape = 20) +
  scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  scale_fill_brewer("Restricted", palette = "Set1", direction = -1) +
  coord_fixed(xlim = c(210, 535), ylim = c(5740, 6050)) +
  facet_wrap(~year_pair) +
  # xlab("Easting (km)") + ylab("Northing (km)")+
  gfplot::theme_pbs() + theme(legend.position=c(0.45,0.2),
    # strip.text.x = element_text(size = 12),
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()
  )
ggsave("figs/hbll-grid-by-year.png", width = 7, height = 7)


# both? ------------

all <- trawl %>%
  filter(survey_abbrev != "SYN WCVI") %>%
  bind_rows(mutate(hbll, survey_abbrev = "HBLL OUT N"))

all_rest <- filter(all, restricted)
dat_rest <- filter(dat, restricted)
dat_hbll_rest <- filter(dat_hbll, restricted)

dat_hbll <- mutate(dat_hbll, survey_abbrev = "HBLL OUT N")

#
# all <- mutate(all, survey_abbrev = ifelse(survey_abbrev %in% c("SYN HS", "SYN QCS"), "SYN QCS, SYN HS", survey_abbrev))
# dat <- mutate(dat, survey_abbrev = ifelse(survey_abbrev %in% c("SYN HS", "SYN QCS"), "SYN QCS, SYN HS", survey_abbrev))
# dat_hbll <- mutate(dat_hbll, survey_abbrev = ifelse(survey_abbrev %in% c("SYN HS", "SYN QCS"), "SYN QCS, SYN HS", survey_abbrev))
# dat_rest <- mutate(dat_rest, survey_abbrev = ifelse(survey_abbrev %in% c("SYN HS", "SYN QCS"), "SYN QCS, SYN HS", survey_abbrev))
# dat_hbll_rest <- mutate(dat_hbll_rest, survey_abbrev = ifelse(survey_abbrev %in% c("SYN HS", "SYN QCS"), "SYN QCS, SYN HS", survey_abbrev))

# r_col <- "#B2182B"
r_col <- "black"

# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 3, 5, 4)]
# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(1:4)]
# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 3, 5, 4)]
# .pal <- RColorBrewer::brewer.pal(5, "Dark2")[c(2, 3, 5, 4)]
# .pal <- RColorBrewer::brewer.pal(4, "Paired")

# .pal <- viridisLite::mako(4, begin = 0.3, end = 0.95)
# .pal <- viridisLite::plasma(4, begin = 0.2, end = 0.85)
# .pal <- c(
#   # "#648FFF",
#           "#785EF0",
#           "#DC267F",
#           "#FE6100", "#FFB000")
.pal <- viridisLite::viridis(4, begin = 0, end = 0.95)
# names(.pal) <- c("SYN HS", "SYN QCS", "HBLL OUT N", "SYN WCHG")
# names(.pal) <- c("SYN QCS",  "HBLL OUT N", "SYN WCHG","SYN HS")
source("analysis/theme.R")
.pal <- c(restricted_cols, restricted_cols[3])

.pal <- RColorBrewer::brewer.pal(4, "Set2")[c(2, 3, 1, 4)]

names(.pal) <-
  c(
    "SYN WCHG",
    "HBLL OUT N",
    "SYN QCS",
    "SYN HS")


pt_size <- 0.15
all %>%
    ggplot() +
    geom_tile(aes(X, Y, fill = survey_abbrev), alpha = 0.2, width=2, height=2) +
    geom_tile(aes(X, Y), alpha = 1, width=2, height=2, data = all_rest, fill = "black") +
    # geom_tile(aes(X, Y, fill = restricted), alpha = 0.4, width=2, height=2) +
    geom_polygon(
      data = coast, aes(x = X, y = Y, group = PID),
      fill = "grey70", col = "grey70", lwd = 0
    ) +
    geom_point(data = dat, aes(X, Y, colour = survey_abbrev), size = pt_size) +
    geom_point(data = dat_hbll, aes(X, Y, colour = survey_abbrev), size = pt_size) +
    # geom_point(data = dat_rest, aes(X, Y), size = pt_size * 1, col = "white") +
    geom_point(data = dat_rest, aes(X, Y, colour = survey_abbrev), size = pt_size * 0.5) +
    # geom_point(data = dat_hbll_rest, aes(X, Y), size = pt_size * 1, col = "white") +
    geom_point(data = dat_hbll_rest, aes(X, Y, colour = survey_abbrev), size = pt_size * 0.5) +
    # scale_fill_brewer("Restricted", palette = "Set3", direction = -1) +
    scale_fill_manual("Survey", values = .pal) +
    # scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
    scale_colour_manual("Survey", values = .pal) +
    coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
    theme(legend.position=c(0.15,0.15)) +
    xlab("Easting (km)") + ylab("Northing (km)")

# ggsave("figs/restricted-grid.png", width = 6, height = 6, dpi = 180)
# ggsave("figs/restricted-grid.pdf", width = 6, height = 6)

# pretty version?

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada")
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
    c(xmin = -134, ymin = 50, xmax = -127, ymax = 55))))
utm_zone9 <- 3156
bc_coast_proj <- sf::st_transform(bc_coast, crs = utm_zone9)

.padding <- 1000
.xlim <- range(all$X) * 1000 + c(-.padding, .padding)
.ylim <- range(all$Y) * 1000 + c(-.padding, .padding)

make_map <- function(.dat) {
  ggplot() +
    geom_tile(data = all, mapping = aes(X * 1000, Y * 1000, fill = survey_abbrev), alpha = 0.3, width = 2000, height = 2000) +
    geom_tile(data = .dat, mapping = aes(X * 1000, Y * 1000), alpha = 1, width = 2000, height = 2000, fill = "black") +
    geom_point(data = dat, aes(X * 1000, Y * 1000, colour = survey_abbrev), size = pt_size) +
    geom_point(data = dat_hbll, aes(X * 1000, Y * 1000, colour = survey_abbrev), size = pt_size) +
    geom_point(data = dat_rest, aes(X * 1000, Y * 1000, colour = survey_abbrev), size = pt_size * 0.5, alpha = 0.6) +
    geom_point(data = dat_hbll_rest, aes(X * 1000, Y * 1000, colour = survey_abbrev), size = pt_size * 0.5, alpha = 0.6) +
    geom_tile(data = .dat, mapping = aes(X * 1000, Y * 1000), width = 2000, height = 2000, fill = "#00000050") +
    theme_light() +
    geom_sf(data = bc_coast_proj, colour = "grey40", fill = "grey80", linewidth = 0.25) +
    coord_sf(xlim = .xlim, ylim = .ylim, crs = 3156) +
    labs(fill = "Predicted\ndensity") +
    labs(x = "Longitude", y = "Latitude") +
    scale_fill_manual("Survey", values = .pal) +
    scale_colour_manual("Survey", values = .pal) +
    theme(legend.position = c(0.14, 0.15)) +
    annotate(geom = "text", x = max(.xlim) - 100000, y = max(.ylim) - 50000, label = "British Columbia", size = 4, colour = "grey10", hjust = 0, vjust = 1) +
    annotate(geom = "text", x = max(.xlim) - 310000, y = max(.ylim) - 85000, label = "Haida Gwaii", size = 4, colour = "grey10", hjust = 0.5, vjust = 0.5) +
    ggspatial::annotation_scale(
      location = "bl",
      pad_x = unit(1.5, "in"), pad_y = unit(0.2, "in"),
      bar_cols = c("grey80", "white"), line_width = 0.5
    ) +
    ggspatial::annotation_north_arrow(
      location = "tr", which_north = "true",
      pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm"),
      style = ggspatial::north_arrow_nautical(
        fill = c("grey40", "white"),
        line_col = "grey20"
      ))
}
g <- make_map(all_rest)
ggsave("figs/restricted-grid.pdf", width = 5.4, height = 5.4)

# nsb ----------------
nsb <- sf::read_sf("data-raw/mpatt_survey_overlaps.gdb", type = 7, layer = "MPATT_Q1_full_march2023_Cat1_Cat2_GH_singlepart")
# unique(x$SurveyOverlap)
nsb <- sf::st_transform(nsb, crs = 3156)

lu <- tribble(
  ~SurveyOverlap,  ~SO,
  "Category 1", "Category 1",
  "Category 2", "Category 2",
  "Gwaii Haanas Site- Multi Use - RCA overlap", "Gwaii Haanas Multi Use",
  "Gwaii Haanas Site- Multi Use", "Gwaii Haanas Multi Use",
  "Gwaii Haanas Site- Strict Protection", "Gwaii Haanas Strict",
  "Gwaii Haanas Site- Strict Protection - RCA overlap", "Gwaii Haanas Strict+RCA"
)
nsb <- left_join(nsb, lu)

pal <- c(as.character(colorBlindness::availableColors())[-1], c("grey60"))[c(1, 4, 3, 5, 2)]


.dat <- all_rest

pt_size <- 0.25
g2 <- ggplot() +
  geom_sf(data = nsb, linewidth = 0.55, mapping = aes(colour = SO, fill = SO)) +
  geom_sf(data = bc_coast_proj, colour = "grey40", fill = "grey80", linewidth = 0.25) +
  theme_light() +
  geom_tile(data = all, mapping = aes(X * 1000, Y * 1000), width = 2000, height = 2000, fill = NA, colour = "#00000020", linewidth = 0.2) +
  coord_sf(xlim = .xlim, ylim = .ylim, crs = 3156) +
  theme() +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  # geom_point(data = dat, aes(X * 1000, Y * 1000), size = pt_size, pch = 21, alpha = 0.2, fill = NA, colour = "grey10") +
  # geom_point(data = dat_hbll, aes(X * 1000, Y * 1000), size = pt_size, pch = 21, alpha = 0.2, fill = NA, colour = "grey10") +
  # geom_point(data = dat_rest, aes(X * 1000, Y * 1000), size = pt_size * 0.5, alpha = 0.6) +
  # geom_point(data = dat_hbll_rest, aes(X * 1000, Y * 1000), size = pt_size * 0.5, alpha = 0.6)
  labs(fill = "Zone category", colour = "Zone category") +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = c(0.21, 0.17)) +
  annotate(geom = "text", x = max(.xlim) - 100000, y = max(.ylim) - 50000, label = "British Columbia", size = 4, colour = "grey10", hjust = 0, vjust = 1) +
  annotate(geom = "text", x = max(.xlim) - 310000, y = max(.ylim) - 85000, label = "Haida Gwaii", size = 4, colour = "grey10", hjust = 0.5, vjust = 0.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    ))

g0 <- cowplot::plot_grid(g2,g)
ggsave("figs/fig1.png", width = 10.5, height = 5.4, dpi = 150, plot = g0)
ggsave("figs/fig1.pdf", width = 10.5, height = 5.4, plot = g0)



###########################################

# asis
trawl_asis <- readRDS("data-generated-ALL/syn-grid-w-restr.rds") %>% filter(year == 2018) %>% select(-year) |> filter(restricted)
hbll_asis <- readRDS("data-generated-ALL/hbll-n-grid-w-restr.rds") %>% filter(year == 2019) %>% select(-year)|> filter(restricted)
all_asis <- bind_rows(trawl_asis, hbll_asis)
asis <- all_asis |> anti_join(select(all_rest, X, Y))

g <- make_map(asis)
ggsave("figs/asis-grid.pdf", width = 5.4, height = 5.4)

g <- all %>%
  ggplot() +
  geom_tile(aes(X, Y, fill = survey_abbrev), alpha = 0.45, width=2, height=2) +
  geom_tile(aes(X, Y), alpha = 0.75, width=2, height=2, data = asis, fill = r_col) +
  # geom_tile(aes(X, Y, fill = restricted), alpha = 0.45, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y, colour = survey_abbrev), size = 0.1) +
  geom_point(data = dat_hbll, aes(X, Y, colour = survey_abbrev), size = 0.1) +
  # geom_point(data = dat_rest, aes(X, Y), size = 0.1, col = r_col) +
  # geom_point(data = dat_hbll_rest, aes(X, Y), size = 0.1, col = r_col) +
  # scale_fill_brewer("Restricted", palette = "Set3", direction = -1) +
  scale_fill_manual("Survey", values = .pal) +
  # scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  scale_colour_manual("Survey", values = .pal) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  theme(legend.position=c(0.15,0.15)) +
  xlab("Easting (km)") + ylab("Northing (km)")

ggsave("figs/asis-grid.png", width = 6, height = 6, dpi = 180)
ggsave("figs/asis-grid.pdf", width = 6, height = 6)

# without the proposed restrictions
# add existing restrictions

utm_zone9 <- 3156

Gwaii_Haanas <- sf::read_sf("data-raw/Gwaii Haanas Land-Sea-People plan FINAL ZONING_Nov 2018")
Gwaii_Haanas <- sf::st_cast(Gwaii_Haanas, "MULTIPOLYGON") %>% sf::st_transform(utm_zone9)

Hecate <- sf::read_sf("data-raw/Hecate") %>%
  sf::st_transform(utm_zone9)

RCA <- sf::read_sf("data-raw/RCA2019") %>%
  sf::st_transform(utm_zone9)
# unique(RCA$NAME)


(g <-all %>% mutate(x = X, y = Y, X = X*1000, Y = Y*1000) %>%
  sf::st_as_sf(coords = c("X", "Y"), crs = utm_zone9) %>%
    ggplot() +
    geom_tile(aes(x*1000, y*1000, fill = survey_abbrev), alpha = 0.45, width=2*1000, height=2*1000) +
    geom_polygon(
      data = coast, aes(x = X*1000, y = Y*1000, group = PID),
      fill = "grey67", col = "grey60", lwd = 0.2
    ) +

    geom_point(data = dat, aes(X*1000, Y*1000, colour = survey_abbrev), size = 0.1) +
    geom_point(data = dat_hbll, aes(X*1000, Y*1000, colour = survey_abbrev), size = 0.1) +
    geom_sf(
      data = filter(Gwaii_Haanas, Zone_Type == "Strict Protection (IUCN II)"), aes(group = Name),
      fill = NA, lwd = 0.3, col = "red"
    ) +
    geom_sf(
      data = filter(Hecate, Zone_Type %in% c("AMZ", #"CPZ" # AMZ = buffer, CPZ = main protected area?
                                             )), aes(group = Shape_Area),
      fill = NA, lwd = 0.3, col = "red"
    ) +
    geom_sf(
      data = filter(RCA, !NAME %in% c(# "Scott Islands",
        "Triangle Island")),
      aes(group = NAME),
      fill = NA, lwd = 0.3, col = "red"
    ) +
    scale_fill_manual("Survey", values = .pal) +
    scale_colour_manual("Survey", values = .pal) +
    coord_sf(xlim = c(180*1000, 590*1000), ylim = c(5640*1000, 6050*1000)) +
    gfplot::theme_pbs() + theme(legend.position=c(0.15,0.15)) +
    xlab("Easting (km)") + ylab("Northing (km)"))

ggsave("figs/map-existing-restrictions.png", width = 6, height = 6)


# gg <- g + facet_wrap(~year)
# ggsave("figs/restricted-grid-year.png", width = 8, height = 8)

all %>%
  ggplot() +
  geom_tile(aes(X, Y, fill = restricted), alpha = 0.45, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y, colour = restricted), size = 0.1) +
  geom_point(data = dat_hbll, aes(X, Y, colour = restricted), size = 0.1) +
  scale_fill_brewer("Restricted", palette = "Set1", direction = -1) +
  scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.15,0.15)) +
  xlab("Easting (km)") + ylab("Northing (km)") +
  facet_wrap(~survey_abbrev)

# ggsave("figs/restricted-grid2.png", width = 6, height = 6)
ggsave("figs/restricted-grid3.png", width = 8, height = 8)

all %>% group_by(survey_abbrev) %>%
  summarise(rest = mean(restricted)) %>%
  knitr::kable()

# about equal to CV effect?
# Q: are these 'good' cells in general for fish?
