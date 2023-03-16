library(dplyr)

# prep index data ----

y1 <- readRDS(file = "data-generated/index-hbll-geo-clean.rds") %>%
  # 0.0024384 * 0.009144 * 2 is area swept in km2 per hook
  # = est * (1/10000)/(0.0024384 * 0.009144 * 2) to convert from # fish / hook to 1000 fish / km2
  mutate(est = est * 2.242481, lwr = lwr * 2.242481, upr = upr * 2.242481)
# y1 <- readRDS(file = "data-generated/index-hbll-geo-clean-binomial-gamma.rds") %>%
#   mutate(est = est/10000, lwr = lwr/10000, upr = upr/10000)
y2 <- readRDS(file = "data-generated/index-syn-geo-clean.rds") %>%
  # current offset in strange units so /10000 gets us to tonnes/km2
  mutate(est = est / 10000, lwr = lwr / 10000, upr = upr / 10000)
y <- bind_rows(y1, y2)

y <- y %>% dplyr::filter(type != "MPA only")
y <- y %>% dplyr::filter(type != "MPA only restricted")

ocv <- y |>
  filter(type == "Status quo") |>
  group_by(species_common_name, survey_abbrev) |>
  summarise(orig_cv_mean = mean(cv)) |>
  select(orig_cv_mean, species_common_name, survey_abbrev) |>
  ungroup() |>
  distinct()

y <- left_join(y, ocv)

mean(y$orig_cv_mean < 1)
filter(y, orig_cv_mean > 1)
filter(y, orig_cv_mean <= 1)
index <- filter(y, orig_cv_mean < 1)

index$species_common_name <- stringr::str_to_title(index$species_common_name)
index <- filter(index, !species_common_name == "Shortraker Rockfish")

# remove Shortspine as only HBLL spp below 10% occurance that converges
# index <- filter(index, !(survey_abbrev == "HBLL OUT N" & species_common_name == "Shortspine Thornyhead"))
index <- mutate(index, species_common_name = gsub("Rougheye/Blackspotted Rockfish Complex", "Rougheye/Blackspotted Rockfish", species_common_name))

# choose labels for each index type ----
index$type_label <- index$type
index[index$type == "Restricted", ]$type_label <- "Extrapolated"
index[index$type == "Restricted and shrunk", ]$type_label <- "Shrunk"
index$type_label <- factor(index$type_label, levels = c("Status quo", "Extrapolated", "Shrunk", "MPA only"))

saveRDS(index, file = "data-generated/index-filtered.rds")
