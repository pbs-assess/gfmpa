library(dplyr)
library(ggplot2)
# library(future)
# plan(multisession)

survey_data <- readRDS("data-generated/dat_to_fit.rds")
hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
survey_data <- bind_rows(survey_data, hbll)

down_sample <- function(x, seed = 1) {
  set.seed(seed * x$year[1])
  number_restricted <- sum(x$restricted)
  n_per_stratum <- round(number_restricted / length(unique(x$grouping_code)))
  # available_rows <- seq_len(nrow(x))
  # sampled <- sample(available_rows, size = number_restricted, replace = FALSE)
  # new <- x[-sampled, , drop = FALSE] # random removal of same N
  remove <- group_by(x, grouping_code) |> slice_sample(n = n_per_stratum) |> ungroup()
  new <- anti_join(x, remove, by = join_by(species_common_name, survey_abbrev, year, X, Y, fishing_event_id))
  new |> mutate(downsample_seed = seed)
}

z1 <- group_by(survey_data, species_common_name, survey_abbrev, year) %>%
  group_split() %>%
  purrr::map_dfr(down_sample, seed = 1)

z2 <- group_by(survey_data, species_common_name, survey_abbrev, year) %>%
  group_split() %>%
  purrr::map_dfr(down_sample, seed = 2)

z3 <- group_by(survey_data, species_common_name, survey_abbrev, year) %>%
  group_split() %>%
  purrr::map_dfr(down_sample, seed = 3)

z4 <- group_by(survey_data, species_common_name, survey_abbrev, year) %>%
  group_split() %>%
  purrr::map_dfr(down_sample, seed = 4)

z5 <- group_by(survey_data, species_common_name, survey_abbrev, year) %>%
  group_split() %>%
  purrr::map_dfr(down_sample, seed = 5)

# plan(sequential)

s <- "SYN WCHG"
s <- "HBLL OUT N"
sp <- "arrowtooth flounder"

z2 |>
  filter(survey_abbrev == s, species_common_name == sp) |>
  ggplot(aes(longitude, latitude, size = if (s == "HBLL OUT N") density_ppkm2 else density_kgpm2)) +
  geom_jitter(width = 0.05, height = 0.05, pch = 21) +
  facet_wrap(~year) +
  scale_size_area() +
  geom_point(data = filter(survey_data, survey_abbrev == s, species_common_name == sp, restricted == TRUE), colour = "red", pch = 21)

nrow(z2) / nrow(survey_data)

bind_rows(list(z1, z2, z3, z4, z5)) |>
  saveRDS("data-generated/downsampled-fitting-data.rds")

# ggsave("filename.pdf", width = 6, height = 6)
