library(dplyr)
library(ggplot2)
ggplot2::theme_set(ggsidekick::theme_sleek())

survey <- "SYN WCHG"
grid <- readRDS("data-generated/grids-strata-restricted.rds")
grid <- filter(grid, survey_abbrev %in% survey)
gr_full <- dplyr::filter(grid, survey_abbrev %in% survey)
gr_mpa <- dplyr::filter(grid, survey_abbrev %in% survey, restricted == TRUE)
gr_remaining <- dplyr::filter(grid, !restricted, survey_abbrev %in% survey)

dat <- readRDS("data-generated/dat_to_fit.rds")
dat <- filter(dat, species_common_name %in% "pacific ocean perch",
  survey_abbrev %in% survey, year == 2018)
dat

# with data

pal <- c("grey50", "red")
names(pal) <- c("FALSE", "TRUE")
g_full <- ggplot(gr_full, aes(X, Y, fill = restricted)) + geom_tile(width = 2, height = 2) +
  coord_fixed() +
  labs(fill = "MPA\n(restriced survey)", x = "", y = "")  +
  scale_fill_manual(values = pal) +
  guides(fill = "none") +
  theme_void() +
  geom_point(data = dat, aes(size = density_kgpm2), pch = 20) +
  scale_size_area(max_size = 13, guide = "none")
g_full
ggsave("figs/grid-wchg-full-eg-dat.png", width = 4, height = 4)

g_rem <- ggplot(gr_remaining, aes(X, Y, fill = restricted)) + geom_tile(width = 2, height = 2) +
  coord_fixed() +
  labs(fill = "MPA\n(restriced survey)", x = "", y = "")  +
  scale_fill_manual(values = pal) +
  guides(fill = "none") +
  theme_void()  +
  geom_point(data = filter(dat, !restricted), aes(size = density_kgpm2), pch = 20) +
  scale_size_area(max_size = 13, guide = "none")
g_rem
ggsave("figs/grid-wchg-remain-eg-dat.png", width = 4, height = 4)

# without

pal <- c("grey50", "red")
names(pal) <- c("FALSE", "TRUE")
g_full <- ggplot(gr_full, aes(X, Y, fill = restricted)) + geom_tile(width = 2, height = 2) +
  coord_fixed() +
  labs(fill = "MPA\n(restriced survey)", x = "", y = "")  +
  scale_fill_manual(values = pal) +
  guides(fill = "none") +
  theme_void()
ggsave("figs/grid-wchg-full-eg.png", width = 4, height = 4)

g_rem <- ggplot(gr_remaining, aes(X, Y, fill = restricted)) + geom_tile(width = 2, height = 2) +
  coord_fixed() +
  labs(fill = "MPA\n(restriced survey)", x = "", y = "")  +
  scale_fill_manual(values = pal) +
  guides(fill = "none") +
  theme_void()
ggsave("figs/grid-wchg-remain-eg.png", width = 4, height = 4)

