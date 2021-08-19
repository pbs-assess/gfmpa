library(dplyr)
library(ggplot2)
library(sf)
library(future)
plan(multisession)
options(future.rng.onMisuse = "ignore")

source("load-data.R")
source("functions.R")

index_syn <- dat_to_fit %>%
  group_by(survey_abbrev, species_common_name) %>%
  group_split() %>%
  furrr::future_map_dfr(function(.x) {
    # cat(.x$species_common_name[1], "\n")
    set.seed(1)
    out <- split(.x, .x$survey_abbrev) %>%
      purrr::map_dfr(function(.y)
        boot_wrapper(.y, reps = 200L), .id = "survey_abbrev") %>%
      mutate(type = "Status quo")
    set.seed(1)
    out_restr <- .x %>%
      filter(!restricted) %>%
      split(.$survey_abbrev) %>%
      purrr::map_dfr(function(.y)
        boot_wrapper(.y, reps = 200L), .id = "survey_abbrev") %>%
      mutate(type = "Restricted")
    combined <- bind_rows(out, out_restr)
    combined$species_common_name <- .x$species_common_name[1]
    combined
  }, .progress = TRUE)

saveRDS(index_syn, file = "index-syn-x-empty.rds")

g <- ggplot(index_syn, aes(year, biomass, ymin = lwr, ymax = upr, colour = type, fill = type)) +
  geom_line() +
  facet_grid(species_common_name~survey_abbrev, scales = "free_y") +
  geom_ribbon(alpha = 0.2, colour = NA)

ggsave("index-syn-restricted-x-empty.pdf", width = 8, height = 85, limitsize = FALSE)

# check:
# ii <- readRDS("/Volumes/Extreme-SSD/src/gfsynopsis-2021/report/data-cache/pacific-cod.rds")$survey_index
#
# ii %>% filter(survey_abbrev %in% c("SYN QCS")) %>%
#   ggplot(aes(year, biomass, ymin = lowerci, ymax = upperci)) +
#   geom_line() +
#   geom_ribbon(alpha = 0.2)

no_na <- index_syn %>%
  group_by(survey_abbrev, species_common_name) %>%
  summarise(nas = sum(is.na(cv) > 0), zeros = sum(median_boot == 0)) %>%
  filter(nas == FALSE & zeros == 0)
nrow(no_na)

nrow(index_syn)
left_join(no_na, index_syn) %>% nrow()

x <- left_join(no_na, index_syn) %>%
  group_by(species_common_name, survey_abbrev, year) %>%
  summarise(cv_ratio = cv[type == "Restricted"] / cv[type == "Status quo"])

x %>%
  ggplot(aes(cv_ratio)) + facet_wrap(~survey_abbrev) + geom_histogram()

x %>% group_by(survey_abbrev) %>%
  summarise(mean_ratio = mean(cv_ratio))

