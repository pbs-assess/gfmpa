library(dplyr)
library(future)
library(sdmTMB)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
# if (!is_rstudio && is_unix) plan(multicore, workers = 6L) else plan(multisession, workers = 6L)
options(future.rng.onMisuse = "ignore")
options(dplyr.summarise.inform = FALSE)
dir.create("figs", showWarnings = FALSE)
if (is_unix) options(sdmTMB.cores = 4L)

# plan(sequential, split = TRUE)
# plan(sequential) # don't crash!

# source("analysis/load-data.R")
source("analysis/functions.R")

# survey <- "SYN"
fam <- "binomial_gamma"

for (survey in c("SYN", "HBLL")) {
  if (fam == "tweedie") {
    family = tweedie()
  } else if (fam == "nbinom2") {
    family = nbinom2()
  } else if (fam == "binomial_gamma") {
    family = sdmTMB::delta_gamma()
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
    ##  hooks <- readRDS("data-raw/HBLLOUTN-hooks.rds") %>%
    ##    mutate(count_animals = count_target_species + count_non_target_species) %>%
    ##    select(year, fishing_event_id, count_animals,
    ##           count_bait_only, count_empty_hooks, count_bent_broken) %>%
    ##    mutate(total_hooks = count_animals + count_bait_only + count_empty_hooks - count_bent_broken) %>%
    ##    mutate(count_bait_only2 = replace(count_bait_only, which(count_bait_only == 0), 1))
    ##
    ##  hookmeans <- filter(hooks, total_hooks > 350) %>%
    ##    summarise(baited  = mean(count_bait_only),
    ##              prop_baited = mean(count_bait_only / total_hooks))
    ##  # hist(hooks[hooks$total_hooks > 350,]$count_bait_only, breaks = 30)
    ##
    ##  dat_to_fit <- left_join(dat_to_fit, hooks)  %>%
    ##    mutate(missing_hooks = hook_count - total_hooks) %>%
    ##    mutate(prop_bait_hooks = ifelse(total_hooks > 350, count_bait_only2 / total_hooks, hookmeans$prop_baited)) %>%
    ##    mutate(hook_adjust_factor = -log(prop_bait_hooks) / (1 - prop_bait_hooks),
    ##           adjusted_hooks = hook_count/hook_adjust_factor
    ##           )

    # d <- filter(dat_to_fit, species_common_name == "pacific cod")
    # hist(d$adjusted_hooks)

    grid <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
    grid$survey_abbrev <- "HBLL OUT N"
    grid$area <- 4

    highlights <- c(
      "arrowtooth flounder",
      "big skate",
      "canary rockfish",
      "china rockfish", # 13% positive
      "copper rockfish", # 11% positive
      "lingcod",
      "longnose skate",
      "north pacific spiny dogfish",
      "pacific cod",
      "quillback rockfish",
      "redbanded rockfish",
      "rosethorn rockfish",
      "rougheye/blackspotted rockfish complex",
      "sandpaper skate", # 7.6 % positive
      "shortspine thornyhead", # only 6.7% of samples are positive, 8% pos restricted
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
    grid$area <- 4

    # "Copper"
    ### indicates species added that Gabe asked about
    syn_highlights <- c(
      "Big Skate",
      "Longnose Skate",
      "sandpaper skate", ###
      "North Pacific Spiny Dogfish",
      "Spotted Ratfish",
      "kelp greenling",
      "Lingcod",
      "Pacific Cod",
      "Walleye Pollock",
      "Bocaccio",
      "Canary Rockfish",
      "darkblotched rockfish", ###
      "greenstriped rockfish", ###
      "Pacific Ocean Perch",
      "Redstripe Rockfish",
      "Redbanded Rockfish",
      "rougheye/blackspotted rockfish complex",
      "sharpchin rockfish", ###
      # "Shortraker Rockfish", # not fitting well! almost no data when restricted (9% pos)
      "Shortspine Thornyhead",
      "Silvergray Rockfish",
      "Widow Rockfish",
      "Yellowmouth Rockfish",
      "Yellowtail Rockfish",
      "Arrowtooth Flounder",
      "Butter Sole", ###
      "Curlfin Sole",
      "Dover Sole",
      "English Sole",
      "Flathead Sole",
      "Petrale Sole",
      "Rex Sole",
      "Slender Sole",
      "southern rock sole"
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

  # dat_to_fit <- dplyr::filter(dat_to_fit, species_common_name == "yellowtail rockfish", survey_abbrev == "SYN QCS")


  if (!file.exists(save_file)) {

    index_restr <- dat_to_fit %>%
      filter(!restricted) %>%
      group_by(survey_abbrev, species_common_name) %>%
      group_split() %>%
      # furrr::future_map_dfr(function(.x) {
        purrr::map_dfr(function(.x) {
        out <- .x %>%
          fit_geo_model(pred_grid = grid, survey = survey, family = family,
            mpa_dat_removed = TRUE) %>%
          mutate(type = "Restricted")
      # }, .progress = TRUE) # %>% filter(region != "mpa")
    })

    index_shrunk <- dat_to_fit %>%
      filter(!restricted) %>%
      group_by(survey_abbrev, species_common_name) %>%
      group_split() %>%
      # furrr::future_map_dfr(function(.x) {
        purrr::map_dfr(function(.x) {
        out <- .x %>%
          fit_geo_model(pred_grid = filter(grid, !restricted), family = family,
            survey = survey, mpa_dat_removed = TRUE, shrunk = TRUE) %>%
          mutate(type = "Restricted and shrunk")
      # }, .progress = TRUE) # %>% filter(region != "mpa")
    })

    index_orig <- dat_to_fit %>%
      group_by(survey_abbrev, species_common_name) %>%
      group_split() %>%
      # furrr::future_map_dfr(function(.x) {
      purrr::map_dfr(function(.x) {
        out <- .x %>%
          fit_geo_model(pred_grid = grid, survey = survey, family = family) %>%
          mutate(type = "Status quo")
        # }, .progress = TRUE)
      })

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
  if (survey == "SYN") saveRDS(y, file = paste0("data-generated/index-syn-geo-clean-", fam, ".rds"))

}
