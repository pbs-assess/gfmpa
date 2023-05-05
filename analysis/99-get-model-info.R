library(dplyr)

f <- list.files("data-generated/model-cache/", pattern = ".rds", )
f <- f[!grepl("model", f)]

get_info <- function(x) {
  spp <- x
  spp <- gsub("\\.rds", "", spp)
  spp <- gsub("-HBLL-OUT-N", "", spp)
  spp <- gsub("-SYN-HS", "", spp)
  spp <- gsub("-SYN-WCHG", "", spp)
  spp <- gsub("-SYN-QCS", "", spp)
  spp <- gsub("-", " ", spp)
  d <- readRDS(paste0("data-generated/model-cache/", x))

  survey <- if (grepl("SYN-HS", x)) {
    "SYN HS"
  } else if  (grepl("SYN-QCS", x)) {
    "SYN QCS"
  } else if  (grepl("SYN-WCHG", x)) {
    "SYN QCS"
  } else if  (grepl("HBLL-OUT-N", x)) {
    "HBLL OUT N"
  } else {
    stop("Survey not found")
  }

  tibble(
    family = if (isTRUE(d$family$delta)) "delta-Gamma" else "Tweedie",
    survey = survey,
    species_science_name = spp,
    share_range = paste(unlist(d$share_range), collapse = ", "),
    spatiotemporal_fields = paste(unlist(d$spatiotemporal), collapse = ", ")
  )
}

info <- purrr::map_dfr(f, get_info)
info$spatiotemporal_fields <- gsub("iid", "IID", info$spatiotemporal_fields)
info$share_range <- gsub("FALSE", "False", info$share_range)
info$share_range <- gsub("TRUE", "True", info$share_range)
info$spatiotemporal_fields <- gsub("off", "None", info$spatiotemporal_fields)
# View(info)

dat_hbll <- readRDS("data-generated/dat_to_fit_hbll.rds")
dat_syn <- readRDS("data-generated/dat_to_fit.rds")

d1 <- select(dat_syn, species_common_name, species_science_name) %>% distinct()
d2 <- select(dat_hbll, species_common_name, species_science_name) %>% distinct()
d <- bind_rows(d1, d2) %>% distinct()

d <- left_join(info, d)

d$species_common_name <- stringr::str_to_title(d$species_common_name)
d$species_science_name <- stringr::str_to_sentence(d$species_science_name)
d$species_science_name <- gsub("/", " ", d$species_science_name)

d <- arrange(d, species_common_name, survey)

View(d)
