library(dplyr)
library(ggplot2)
library(sf)
library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 8L) else plan(multisession)
options(future.rng.onMisuse = "ignore")
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())

# source("load-data.R")
source("functions.R")

dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

geo_wrapper <- function(surv_dat, pred_grid, shrink_survey = FALSE) {
  utm_zone9 <- 3156
  coords <- surv_dat %>%
    st_as_sf(crs = 4326, coords = c("longitude", "latitude")) %>%
    st_transform(utm_zone9) %>%
    sf::st_coordinates() %>%
    as.data.frame()
  coords$X <- coords$X / 1000
  coords$Y <- coords$Y / 1000
  surv_dat <- dplyr::bind_cols(surv_dat, coords)
  surv_dat$density <- surv_dat$density_ppkm2
  mesh <- make_mesh(surv_dat, c("X", "Y"), cutoff = 10)
  # plot(mesh)
  # mesh$mesh$n

  null_df <- data.frame(
    species_common_name = surv_dat$species_common_name[1],
    survey_abbrev = surv_dat$survey_abbrev[1],
    stringsAsFactors = FALSE
  )
  m <- try({
    sdmTMB(density ~ 0 + as.factor(year),
      data = surv_dat,
      family = sdmTMB::tweedie(),
      time = "year",
      spde = mesh
    )
  })
  # file_name <- paste0("data-generated/",
  #   gsub("\\/", "-", null_df$species_common_name[1]), "-",
  #   null_df$species_common_name[1], ".rds"
  # )
  if (class(m)[[1]] == "try-error") {
    return(null_df)
  }
  if (max(m$gradients) > 0.01) {
    m <- try({
      run_extra_optimization(m, newton_loops = 1L, nlminb_loops = 0L)
    })
  }
  # saveRDS(m, file = file_name)
  if (class(m)[[1]] == "try-error") {
    return(null_df)
  }

  set.seed(1)
  pred <- try({
    predict(m, newdata = pred_grid, xy_cols = c("X", "Y"), sims = 250L)
  })
  if (class(pred)[[1]] == "try-error") {
    return(null_df)
  }
  ind <- get_index_sims(pred, area = rep(4, nrow(pred)))
  ind$species_common_name <- surv_dat$species_common_name[1]
  ind$survey_abbrev <- surv_dat$survey_abbrev[1]
  ind$max_gradient <- max(m$gradients)
  ind
}

hbll_grid <- gfplot::hbll_n_grid$grid
utm_zone9 <- 3156
coords <- hbll_grid %>%
  sf::st_as_sf(crs = 4326, coords = c("X", "Y")) %>%
  sf::st_transform(utm_zone9)
coords_restr <- shrink_a_survey(coords, ll_removed, plot = FALSE)
coords <- coords %>%
  sf::st_coordinates() %>%
  as.data.frame()
coords$X <- coords$X / 1000
coords$Y <- coords$Y / 1000
coords$restricted <- coords_restr$restricted
hbll_grid <- coords %>%
  expand_prediction_grid(years = unique(dat_to_fit$year)) %>%
  as_tibble()

index_hbll_orig <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
  # purrr::map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = hbll_grid) %>%
      mutate(type = "Status quo")
  }, .progress = TRUE)

# library(progressr)
# handlers(global = TRUE)
# handlers("progress")
# p <- progressor(steps =
#     distinct(select(dat_to_fit, survey_abbrev, species_common_name)) %>% nrow())

index_hbll_restr <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = hbll_grid) %>%
      mutate(type = "Restricted")
  }, .progress = TRUE)

index_hbll_shrunk <- dat_to_fit %>%
  filter(!restricted) %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    out <- .x %>%
      geo_wrapper(pred_grid = filter(hbll_grid, !restricted)) %>%
      mutate(type = "Restricted and shrunk")
  }, .progress = TRUE)

index_hbll <- bind_rows(index_hbll_orig, index_hbll_restr, index_hbll_shrunk)
saveRDS(index_hbll, file = "index-hbll-geo.rds")

plan(sequential)
index_hbll <- readRDS("index-hbll-geo.rds")
index <- filter(index_hbll, !is.na(est), !is.na(se))

max(index$max_gradient)

index <- index %>% mutate(cv = sqrt(exp(se^2) - 1))

# how many didn't fit in orig?

.u1 <- distinct(index, survey_abbrev, species_common_name, type) %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(orig = "Status quo" %in% type, restr = "Restricted" %in% type)

# didn't converge?
sum(.u1$orig)
sum(.u1$restr)

.u1 <- filter(.u1, orig & restr)
index <- left_join(.u1, index)

y <- index %>% group_by(survey_abbrev, species_common_name) %>%
  mutate(orig_cv = max(cv[type == "Status quo"]))

saveRDS(y, file = "index-hbll-geo-clean.rds")

mean(y$orig_cv < 1)
filter(y, orig_cv > 1)
index <- filter(y, orig_cv < 1)

g <- ggplot(index, aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_ribbon(alpha = 0.2, colour = NA)

ggsave("figs/index-hbll-geo-restricted.pdf", width = 18, height = 18, limitsize = FALSE)

x <- index %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(cv_ratio_restr = cv[type == "Restricted"] / cv[type == "Status quo"],
    cv_ratio_shrunk = cv[type == "Restricted and shrunk"] / cv[type == "Status quo"])

x %>%
  tidyr::pivot_longer(starts_with("cv"), names_to = "Restriction type", values_to = "CV ratio") %>%
  ggplot(aes(`CV ratio`)) + facet_wrap(~survey_abbrev) + geom_histogram() +
  geom_vline(xintercept = 1, col = "red") +
  facet_wrap(vars(`Restriction type`)) +
  coord_cartesian(expand = FALSE)

ggsave("figs/index-geo-hbll-se-ratios.pdf", width = 8, height = 4, limitsize = FALSE)

x_long <- x %>%
  tidyr::pivot_longer(starts_with("cv"), names_to = "Restriction type", values_to = "CV ratio")

x_long %>% group_by(`Restriction type`) %>%
  summarise(mean_ratio = mean(`CV ratio`))

lu <- tibble("Restriction type" = c("cv_ratio_restr", "cv_ratio_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain"))
x_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(`CV ratio`), upr = max(`CV ratio`), est = mean(`CV ratio`)) %>%
  left_join(lu) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_pointrange(position = position_dodge(width = 0.5)) + coord_flip() +
  xlab("") + ylab("Ratio of index CV\n(restricted/original)") + labs(colour = "Restricted survey domain") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
ggsave("figs/index-geo-hbll-cv-ratio-dotplot.pdf", width = 7, height = 7)
ggsave("figs/index-geo-hbll-cv-ratio-dotplot.png", width = 7, height = 7)

plan(sequential)

# # ??
# index <- filter(index, !(species_common_name == "deepsea sole" & survey_abbrev == "SYN WCHG" & year == 2020))

x <- index %>%
  group_by(species_common_name, survey_abbrev, type) %>%
  mutate(est = est / exp(mean(log(est)))) %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(re_restr = (est[type == "Restricted"] - est[type == "Status quo"]) / est[type == "Status quo"],
  re_shrunk = (est[type == "Restricted and shrunk"] - est[type == "Status quo"]) / est[type == "Status quo"]
  )

x_long <- x %>%
  tidyr::pivot_longer(starts_with("re"), names_to = "Restriction type", values_to = "re")

g <- ggplot(x_long, aes(year, re, colour = `Restriction type`)) +
  geom_line() +
  facet_wrap(~species_common_name, scales = "free_y", ncol = 4) +
  geom_hline(yintercept = 0, lty = 2)

ggsave("figs/index-hbll-geo-restricted-re.pdf", width = 10, height = 9)

lu <- tibble("Restriction type" = c("re_restr", "re_shrunk"),
  restr_clean = c("Same survey domain", "Shrunk survey domain"))

x_long %>%
  group_by(survey_abbrev, species_common_name, `Restriction type`) %>%
  summarise(lwr = min(re), upr = max(re), est = median(abs(re))) %>%
  left_join(lu) %>%
  ggplot(aes(forcats::fct_reorder(stringr::str_to_title(species_common_name), est), est, colour = restr_clean)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  geom_point(position = position_dodge(width = 0)) + coord_flip() +
  xlab("") + ylab("Median absolute relative error (MARE)\n(restricted compared to original)") + labs(colour = "Restricted survey domain") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top") +
  theme(panel.grid.major.y = element_line(colour = "grey90"))
ggsave("figs/index-geo-hbll-mare-dotplot.pdf", width = 7, height = 7)
ggsave("figs/index-geo-hbll-mare-dotplot.png", width = 7, height = 7)

# ggplot(x, aes(year, re, group = species_common_name)) +
#   geom_line(alpha = 0.5) +
#   facet_grid(~survey_abbrev, scales = "free_y") +
#   geom_hline(yintercept = 0, lty = 2) +
#   coord_cartesian(ylim = c(-0.5, 0.5))
#
# x %>% group_by(species_common_name, survey_abbrev) %>%
#   summarise(mare = median(abs(re))) %>%
#   ungroup() %>%
#   arrange(-mare) %>%
#   top_n(25) %>%
#   as.data.frame()
#
# x %>% group_by(species_common_name, survey_abbrev) %>%
#   summarise(mare = median(abs(re))) %>%
#   ggplot(aes(mare)) + geom_histogram(binwidth = 0.05) + facet_grid(~survey_abbrev) +
#   coord_cartesian(xlim = c(0, 2))
