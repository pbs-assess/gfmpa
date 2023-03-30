library(dplyr)
library(ggplot2)
library(future)
plan(multisession)

survey_data <- readRDS("data-generated/dat_to_fit.rds")
hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
survey_data <- bind_rows(survey_data, hbll)

up_sample <- function(x, seed = 1) {
  set.seed(seed)
  number_restricted <- sum(x$restricted)
  remaining_samples <- filter(x, !x$restricted)
  available_rows <- seq_len(nrow(remaining_samples))
  sampled <- sample(available_rows, size = number_restricted, replace = TRUE)
  new <- remaining_samples[sampled, , drop = FALSE]
  new$upsampled <- TRUE
  old <- remaining_samples
  old$upsampled <- FALSE
  bind_rows(new, old) |> mutate(upsample_seed = seed)
}

z1 <- group_by(survey_data, grouping_code, species_common_name, year) %>%
  group_split() %>%
  furrr::future_map_dfr(up_sample, seed = 1, .options = furrr::furrr_options(seed = TRUE))

z2 <- group_by(survey_data, grouping_code, species_common_name, year) %>%
  group_split() %>%
  furrr::future_map_dfr(up_sample, seed = 2, .options = furrr::furrr_options(seed = TRUE))

z3 <- group_by(survey_data, grouping_code, species_common_name, year) %>%
  group_split() %>%
  furrr::future_map_dfr(up_sample, seed = 3, .options = furrr::furrr_options(seed = TRUE))

plan(sequential)

s <- "SYN WCHG"
s <- "HBLL OUT N"
sp <- "arrowtooth flounder"

z2 |>
  filter(survey_abbrev == s, species_common_name == sp) |>
  arrange(upsampled) |>
  ggplot(aes(longitude, latitude, color = upsampled, size = if (s == "HBLL OUT N") density_ppkm2 else density_kgpm2)) +
  geom_jitter(width = 0.05, height = 0.05, pch = 21) +
  facet_wrap(~year) +
  scale_size_area() +
  geom_point(data = filter(survey_data, survey_abbrev == s, species_common_name == sp, restricted == TRUE), colour = "black", pch = 21)

# bind_rows(list(z1, z2, z3)) |>
bind_rows(list(z1)) |>
  saveRDS("data-generated/upsampled-fitting-data.rds")
