# mapping MPAs
library(dplyr)
library(ggplot2)
library(sf)


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
    year %in% c(2019, 2020) ~ "2019-2020"  # no WCVI
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


ggplot(trawl) +
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
        year %in% c(2019, 2020) ~ "2019-2020"  # no WCVI
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

