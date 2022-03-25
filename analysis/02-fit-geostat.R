library(dplyr)
library(ggplot2)
library(sf)
library(future)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
if (!is_rstudio && is_unix) plan(multicore, workers = 5L) else plan(multisession, workers = 5L)
options(future.rng.onMisuse = "ignore")
library(sdmTMB)
theme_set(ggsidekick::theme_sleek())
options(dplyr.summarise.inform = FALSE)
dir.create("figs", showWarnings = FALSE)

# source("analysis/load-data.R")
source("analysis/functions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L)
  stop("This script is meant to be run from the command line.", call. = FALSE)
survey <- args[[1]]
fam <- args[[2]]

# survey <- "HBLL"
# fam <- "binomial_gamma"

if (fam == "tweedie") {
  family = tweedie()
} else if (fam == "nbinom2") {
  family = nbinom2()
} else if (fam == "binomial_gamma") {
  family = binomial_gamma()
} else {
  stop("Family not found.")
}

cat("survey: ", survey, "\n")
cat("family: ", family$family, "\n")

# # Globals to set ------------------------------
# # survey <- "HBLL"
# survey <- "SYN"
# # family <- binomial_gamma()
# family <- tweedie()
# # family <- nbinom2
# ---------------------------------------------

if (survey == "HBLL") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit_hbll.rds")
  grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
  grid$survey_abbrev <- "HBLL OUT N"

  highlights <- c(
    "arrowtooth flounder",
    "big skate",
    "canary rockfish",
    "china rockfish",
    "copper rockfish",
    "greenstriped rockfish",
    "lingcod",
    "longnose skate",
    "north pacific spiny dogfish",
    "pacific cod",
    "quillback rockfish",
    "redbanded rockfish",
    "rosethorn rockfish",
    "rougheye/blackspotted rockfish complex",
    "sandpaper skate",
    "shortspine thornyhead",
    "silvergray rockfish",
    "southern rock sole",
    "spotted ratfish",
    "tiger rockfish",
    "yelloweye rockfish"
    ) %>% tolower() %>% sort() %>% unique()
  dat_to_fit <- filter(dat_to_fit, species_common_name %in% highlights)
}
if (survey == "SYN") {
  dat_to_fit <- readRDS("data-generated/dat_to_fit.rds")
  grid <- readRDS("data-generated/syn-grid-w-restr.rds")

  # "Copper"
  syn_highlights <- c(
    "Redstripe Rockfish",
    "Big Skate",
    "Walleye Pollock",
    "Pacific Ocean Perch",
    "Yellowmouth Rockfish",
    "Bocaccio",
    "Petrale Sole",
    "Widow Rockfish",
    "Canary Rockfish",
    "Yellowtail Rockfish",
    "Rex Sole",
    "English Sole",
    "rougheye/blackspotted rockfish complex",
    "Dover Sole",
    "Shortspine Thornyhead",
    "Spotted Ratfish",
    "southern rock sole",
    "Pacific Cod",
    "Lingcod",
    "North Pacific Spiny Dogfish",
    "Shortraker Rockfish",
    "Arrowtooth Flounder",
    "Longnose Skate",
    "Redbanded Rockfish",
    "Silvergray Rockfish",
    "Curlfin Sole",
    "Slender Sole",
    "Flathead Sole"
  ) %>% tolower() %>% sort() %>% unique()
  dat_to_fit <- filter(dat_to_fit, species_common_name %in% syn_highlights)
}

# dat_to_fit <- filter(dat_to_fit, survey_abbrev == "SYN QCS")

ll_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")
trawl_removed <- readRDS("data-generated/hu_co_demersalfishing_bottomlongline_d_X.rds")

if (survey == "SYN") {
  save_file <- paste0("data-generated/index-syn-geo-", fam, ".rds")
}
if (survey == "HBLL") {
  save_file <- paste0("data-generated/index-hbll-geo-", fam, ".rds")
}

if (!file.exists(save_file)) {
  index_orig <- dat_to_fit %>%
    group_by(survey_abbrev, species_common_name) %>%
    group_split() %>%
    furrr::future_map_dfr(function(.x) {
      # purrr::map_dfr(function(.x) {
      out <- .x %>%
        fit_geo_model(pred_grid = grid, survey = survey, family = family) %>%
        mutate(type = "Status quo")
    }, .progress = TRUE)
  # })

  index_restr <- dat_to_fit %>%
    filter(!restricted) %>%
    group_by(survey_abbrev, species_common_name) %>%
    group_split() %>%
    furrr::future_map_dfr(function(.x) {
      out <- .x %>%
        fit_geo_model(pred_grid = grid, survey = survey, family = family) %>%
        mutate(type = "Restricted")
    }, .progress = TRUE) # %>% filter(region != "mpa")

  index_shrunk <- dat_to_fit %>%
    filter(!restricted) %>%
    group_by(survey_abbrev, species_common_name) %>%
    group_split() %>%
    furrr::future_map_dfr(function(.x) {
      out <- .x %>%
        fit_geo_model(pred_grid = filter(grid, !restricted), family = family, survey = survey) %>%
        mutate(type = "Restricted and shrunk")
    }, .progress = TRUE) # %>% filter(region != "mpa")

  index_all <- bind_rows(index_orig, index_restr, index_shrunk)
  saveRDS(index_all, file = save_file)
}
plan(sequential) # don't crash!

# gg <- bind_rows(tw, dg) %>% filter(survey_abbrev != "SYN WCHG") %>% ggplot(aes(year, est, ymin = lwr, ymax = upr, colour = type, fill = type)) +
#   geom_line(lwd = 0.9) +
#   geom_ribbon(alpha = 0.2, colour = NA) +
#   labs(x = "Year", colour = "Type", fill = "Type") +  facet_grid(species_common_name~survey_abbrev, scales = "free_y")

index <- readRDS(save_file) %>% filter(!(region == "mpa" & type %in% c("Restricted", "Restricted and shrunk")))

index <- filter(index, !is.na(est), !is.na(se)) # didn't fit
index[index$region == "mpa", ]$type <- "MPA only"
index <- index %>% mutate(cv = sqrt(exp(se^2) - 1))
# stopifnot(max(index$max_gradient) < 0.01)

# how many didn't fit in orig?
.u1 <- distinct(index, survey_abbrev, species_common_name, type) %>%
  group_by(species_common_name, survey_abbrev) %>%
  summarise(
    orig_converged = "Status quo" %in% type,
    restr_converged = "Restricted" %in% type,
    shrunk_converged = "Restricted and shrunk" %in% type
  ) %>%
  ungroup()

table(index$type)

filter(.u1, orig_converged & !restr_converged)
filter(.u1, !orig_converged & restr_converged)
.u1 <- filter(.u1, restr_converged & orig_converged)
index <- left_join(.u1, index)
y <- index %>%
  group_by(survey_abbrev, species_common_name) %>%
  mutate(orig_cv = max(cv[type == "Status quo"]))

if (survey == "HBLL") saveRDS(y, file = paste0("data-generated/index-hbll-geo-clean-", fam, ".rds"))
if (survey == "SYN") saveRDS(y, file = "data-generated/index-syn-geo-clean-", fam, ".rds")
