# mapping MPAs
library(dplyr)
library(sf)

coast <- gfplot:::load_coastline(
  range(df$X) + c(-1, 1),
  range(df$Y) + c(-1, 1),
  utm_zone = 9
)

dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
trawl <- readRDS("data-generated/syn-grid-w-restr.rds") %>% select(-year)

df <- trawl %>%
  dplyr::mutate(X = X*1000, Y = Y*1000) %>%
  gfplot:::utm2ll(., utm_zone = 9)

dat <- dat_to_fit %>%
  dplyr::mutate(X = longitude, Y = latitude) %>%
  gfplot:::ll2utm(., utm_zone = 9)

# ggplot(trawl, aes(X, Y)) +
#   geom_tile(aes(fill = restricted), width=2, height=2) +
#   geom_polygon(
#     data = coast, aes(x = X, y = Y, group = PID),
#     fill = "grey87", col = "grey70", lwd = 0.2
#   )+
#   scale_fill_brewer(palette = "Set1", direction = -1) +
#   coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
#   gfplot::theme_pbs() + xlab("UTM") + ylab("UTM")

dat <- dat %>% filter(species_common_name == "arrowtooth flounder")
ggplot(trawl) +
  geom_tile(aes(X, Y, fill = restricted), width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat, aes(X, Y), size = 0.1) +
  # facet_wrap(~year) +
  scale_fill_brewer("MPA", palette = "Set1", direction = -1) +
  coord_fixed(xlim = c(180, 590), ylim = c(5640, 6050)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.9,0.9)) +
  xlab("UTM") + ylab("UTM")
ggsave("figs/trawl-grid.pdf", width = 6, height = 6)


dat_to_fit_hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
hbll <- readRDS("data-generated/hbll-n-grid-w-restr.rds")

df_hbll <- hbll %>%
  dplyr::mutate(X = X*1000, Y = Y*1000) %>%
  gfplot:::utm2ll(., utm_zone = 9)

dat_hbll <- dat_to_fit_hbll %>%
  dplyr::mutate(X = longitude, Y = latitude) %>%
  gfplot:::ll2utm(., utm_zone = 9)


dat_hbll <- dat_hbll %>% filter(species_common_name == "arrowtooth flounder")
ggplot(hbll) +
  geom_tile(aes(X, Y, fill = restricted), width=2, height=2) +
  geom_polygon(
    data = coast, aes(x = X, y = Y, group = PID),
    fill = "grey87", col = "grey70", lwd = 0.2
  ) +
  geom_point(data = dat_hbll, aes(X, Y), size = 0.15) +
  # facet_wrap(~year) +
  scale_fill_brewer("MPA", palette = "Set1", direction = -1) +
  coord_fixed(xlim = c(210, 535), ylim = c(5740, 6050)) +
  gfplot::theme_pbs() + theme(legend.position=c(0.9,0.9)) +
  xlab("UTM") + ylab("UTM")

ggsave("figs/hbll-grid.pdf", width = 6, height = 6)

