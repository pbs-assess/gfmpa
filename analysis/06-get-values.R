# get stats for manuscript
library(dplyr)

mround <- function(x, digits) {
  sprintf(paste0("%.", digits, "f"), round(x, digits))
}

write_tex <- function(x, macro, ...) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") %>%
    readr::write_lines("analysis/values.tex", append = TRUE)
}

suppressWarnings(file.remove("analysis/values.tex"))
paste0("% lost on each survey") |>
  readr::write_lines("analysis/values.tex", append = TRUE)

.d <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHBLL")

grid <- readRDS("data-generated/syn-grid-w-restr.rds")
.d <- filter(grid, survey_abbrev %in% c("SYN QCS", "SYN HS"))
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostQCSHS")

.d <- filter(grid, survey_abbrev == "SYN WCHG")
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostWCHG")

.d <- filter(grid, survey_abbrev == "SYN QCS")
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostQCS")

.d <- filter(grid, survey_abbrev == "SYN HS")
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHS")

# with as-is where-is
paste0("\n% lost on each survey with as-is where-is removed too") |>
  readr::write_lines("analysis/values.tex", append = TRUE)

.d <- readRDS("data-generated-ALL/hbll-n-grid-w-restr.rds")
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHBLLall")

grid <- readRDS("data-generated-ALL/syn-grid-w-restr.rds")
.d <- filter(grid, survey_abbrev %in% c("SYN QCS", "SYN HS"))
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted))
, 2), "lostQCSHSall")

.d <- filter(grid, survey_abbrev %in% c("SYN QCS"))
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostQCSall")

.d <- filter(grid, survey_abbrev %in% c("SYN HS"))
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHSall")


.d <- filter(grid, survey_abbrev %in% c("SYN WCHG"))
write_tex(mround((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostWCHGall")


paste0("\n% N values") |>
  readr::write_lines("analysis/values.tex", append = TRUE)

metrics <- readRDS("data-generated/metrics-wide2.rds")

nsp <- filter(metrics, est_type == "geostat") |>
  filter(!is.na(mare_med)) |> pull(species_common_name) |> unique() |>
  length()
write_tex(nsp, "nSpp")

nsp <- filter(metrics, est_type == "geostat", survey_abbrev == "HBLL OUT N") |>
  filter(!is.na(mare_med)) |> pull(species_common_name) |> unique() |>
  length()
write_tex(nsp, "hbllNSpp")

nsp <- filter(metrics, est_type == "geostat", grepl("SYN", survey_abbrev)) |>
  filter(!is.na(mare_med)) |> pull(species_common_name) |> unique() |>
  length()
write_tex(nsp, "synNSpp")

# filter(index, survey_abbrev == "HBLL OUT N") %>% pull(species_common_name) %>% unique() %>% length() %>%
#   write_tex("hbllNSpp")
#
# filter(index, survey_abbrev != "HBLL OUT N") %>% pull(species_common_name) %>% unique() %>% length() %>%
#   write_tex("synNSpp")

# paste0("\n% average metric effects") |>
#   readr::write_lines("analysis/values.tex", append = TRUE)

# metrics_wide <- readRDS("data-generated/metrics-wide.rds")
metrics <- readRDS("data-generated/metrics-long2.rds")

# m <- metrics |>
#   mutate(survey_abbrev = as.character(survey_abbrev)) |>
#   mutate(est = ifelse(grepl("trend", measure), abs(est), est)) |>
#   group_by(restr_clean, survey_abbrev, measure) |>
#   summarise(
#     mean_est = mean(est),
#     .groups = "drop"
#   ) |>
#   mutate(digi = ifelse(grepl("precision", measure), 0, 2)) |>
#   mutate(mean_est = mround(mean_est, digi)) |>
#   mutate(restr = ifelse(grepl("^Shrunk", restr_clean), "shrunk", "extrap")) |>
#   select(-restr_clean) |>
#   mutate(measure = ifelse(grepl("precision", measure), "precision", measure)) |>
#   mutate(measure = ifelse(grepl("accuracy", measure), "accuracy", measure)) |>
#   mutate(measure = ifelse(grepl("bias", measure), "bias", measure)) |>
#   mutate(survey_abbrev = ifelse(grepl("QCS", survey_abbrev), "QCSHS", survey_abbrev)) |>
#   mutate(survey_abbrev = ifelse(grepl("WCHG", survey_abbrev), "WCHG", survey_abbrev)) |>
#   mutate(survey_abbrev = ifelse(grepl("HBLL", survey_abbrev), "HBLL", survey_abbrev)) |>
#   mutate(token = paste0(measure, survey_abbrev, restr))
#
# for (i in seq_len(nrow(m))) {
#   write_tex(m$mean_est[i], m$token[i])
# }
#
# paste0("\n% example metric effects") |>
#   readr::write_lines("analysis/values.tex", append = TRUE)
#
# m2 <- metrics |> filter(grepl("^Shrunk", restr_clean))


make_ex_metric <- function(sp, surv, .metric) {
  digi <- 2
  # digi <- if ("cv_perc" %in% .metric) 0 else 2
  m2 |> filter(species_common_name == sp) |>
    filter(grepl(surv, survey_abbrev), grepl(.metric, measure)) |>
    pull(est) |> mround(digi) |>
    write_tex(paste0(tolower(gsub("\\/", "", gsub(" ", "", sp))), surv, measure))
}


# make_ex_metric("Shortspine Thornyhead", "HBLL", "cv_perc")
# make_ex_metric("Rougheye/Blackspotted Rockfish", "WCHG", "precision")
# make_ex_metric("English Sole", "WCHG", "bias")
# make_ex_metric("China Rockfish", "HBLL", "bias")
# make_ex_metric("China Rockfish", "HBLL", "accuracy")
# make_ex_metric("Pacific Ocean Perch", "WCHG", "precision")
# make_ex_metric("Sharpchin Rockfish", "QCS", "accuracy")

paste0("\n % covariate slopes") |>
  readr::write_lines("analysis/values.tex", append = TRUE)

slopes <- readRDS("data-generated/metrics-slopes-table.rds")
slopes <- slopes |> mutate(measure = gsub("cv_perc", "precision", measure)) |>
  mutate(measure = gsub("slope_re", "trend", measure)) |>
  mutate(base_level = gsub("Design-based", "Design", base_level)) |>
  mutate(base_level = gsub("Geostatistical", "Geo", base_level)) |>
  mutate(measure = stringr::str_to_title(measure)) |>
  mutate(label = paste0("cov", measure, base_level))

apply(slopes, 1, function(.x) {
    write_tex(.x[["text"]], .x[["label"]])
  })

#
#   # filter(restr_clean == "Shrunk survey domain") |>
#   # select(-restr_clean)
#
# write_tex(m$mean_mare[m$survey_abbrev == "HBLL OUT N"], "mareHBLL")
# write_tex(m$mean_mare[m$survey_abbrev == "SYN WCHG"], "mareWCHG")
# write_tex(m$mean_mare[m$survey_abbrev == "SYN QCS, SYN HS"], "mareQCSHS")
#
# write_tex(m$mean_cv_ratio[m$survey_abbrev == "HBLL OUT N"], "cvHBLL")
# write_tex(m$mean_cv_ratio[m$survey_abbrev == "SYN WCHG"], "cvWCHG")
# write_tex(m$mean_cv_ratio[m$survey_abbrev == "SYN QCS, SYN HS"], "cvQCSHS")
#
# write_tex(m$mean_slope_re[m$survey_abbrev == "HBLL OUT N"], "slopeAbsReHBLL")
# write_tex(m$mean_slope_re[m$survey_abbrev == "SYN WCHG"], "slopeAbsReWCHG")
# write_tex(m$mean_slope_re[m$survey_abbrev == "SYN QCS, SYN HS"], "slopeAbsReQCSHS")
#
# filter(metrics_wide, species_common_name == "Shortspine Thornyhead",
#   restr_clean == "Shrunk survey domain", survey_abbrev == "HBLL OUT N") |>
#   select(cv_ratio, cv_lwr, cv_upr)
#   # mround(2)
#
# # metrics_wide |>
# #   filter(prop_mpa > 0.1) |>
# #   group_by(restr_clean, survey_abbrev) |>
# #   summarise(
# #     mean_mare = mround(mean(mare), 2),
# #     mean_cv_ratio = mround(mean(cv_ratio), 2),
# #     mean_slope_re = mround(mean(slope_re), 1), .groups = "drop"
# #   ) |>
# #   filter(restr_clean == "Shrunk survey domain") |>
# #   select(-restr_clean)

system("cp analysis/values.tex ~/src/gf-mpa-index-ms/values.tex")
# system("cp figs/*.pdf ~/src/overleaf/gf-mpa-index/figs/new/")

# table of spp ---------

survey_data <- readRDS("data-generated/dat_to_fit.rds")
hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
survey_data <- bind_rows(survey_data, hbll)
survey_data |> select(species_common_name, species_science_name) |>
  mutate(species_science_name = stringr::str_to_sentence(species_science_name)) |>
  mutate(species_common_name = stringr::str_to_title(species_common_name)) |>
  mutate(species_science_name = paste0("\\emph{", species_science_name, "}")) |>
  distinct() |>
  arrange(species_common_name) |>
  knitr::kable(col.names = c("Common name", "Scientific name"), format = "latex", escape = FALSE, booktabs = TRUE, longtable = TRUE, linesep = '', caption = "Groundfish species common and scientific names used throughout this analysis.", label = "spp-sci") |>
  readr::write_lines("../gf-mpa-index-ms/figs/spp-table.tex")
