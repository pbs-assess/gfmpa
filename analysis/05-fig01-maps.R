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

ggplot(hbll) +
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

# r_col <- "#B2182B"
r_col <- "black"

# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 3, 5, 4)]
# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(1:4)]
# .pal <- RColorBrewer::brewer.pal(5, "Set1")[c(2, 3, 5, 4)]
.pal <- RColorBrewer::brewer.pal(5, "Dark2")[c(2, 3, 5, 4)]
# .pal <- RColorBrewer::brewer.pal(4, "Paired")

# .pal <- viridisLite::mako(4, begin = 0.3, end = 0.95)
# .pal <- viridisLite::plasma(4, begin = 0.2, end = 0.85)
# .pal <- c(
#   # "#648FFF",
#           "#785EF0",
#           "#DC267F",
#           "#FE6100", "#FFB000")
.pal <- viridisLite::viridis(4, begin = 0, end = 0.95)
names(.pal) <- c("SYN HS", "SYN QCS", "HBLL OUT N", "SYN WCHG")
# names(.pal) <- c("SYN QCS",  "HBLL OUT N", "SYN WCHG","SYN HS")

(g <- all %>%
  ggplot() +
  geom_tile(aes(X, Y, fill = survey_abbrev), alpha = 0.45, width=2, height=2) +
  geom_tile(aes(X, Y), alpha = 0.45, width=2, height=2, data = all_rest, fill = r_col) +
  # geom_tile(aes(X, Y, fill = restricted), alpha = 0.45, width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y, colour = survey_abbrev), size = 0.1) +
  geom_point(data = dat_hbll, aes(X, Y, colour = survey_abbrev), size = 0.1) +
  geom_point(data = dat_rest, aes(X, Y), size = 0.1, col = r_col) +
  geom_point(data = dat_hbll_rest, aes(X, Y), size = 0.1, col = r_col) +
  # scale_fill_brewer("Restricted", palette = "Set3", direction = -1) +
  scale_fill_manual("Survey", values = .pal) +
  # scale_colour_manual("Restricted", values = c("#2166AC", "#B2182B")) +
  scale_colour_manual("Survey", values = .pal) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.15,0.15)) +
  xlab("Easting (km)") + ylab("Northing (km)"))

ggsave("figs/restricted-grid.png", width = 6, height = 6)

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
