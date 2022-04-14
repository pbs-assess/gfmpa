# get stats for manuscript
library(dplyr)

write_tex <- function(x, macro, ...) {
  paste0("\\newcommand{\\", macro, "}{", x, "}") %>%
    readr::write_lines("analysis/values.tex", append = TRUE)
}


suppressWarnings(file.remove("analysis/values.tex"))
paste0("% lost on each survey") %>% readr::write_lines("analysis/values.tex", append = TRUE)

.d <- readRDS("data-generated/hbll-n-grid-w-restr.rds")
write_tex(round((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHBLL")

grid <- readRDS("data-generated/syn-grid-w-restr.rds")
.d <- filter(grid, survey_abbrev == "SYN HS")
write_tex(round((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostHS")
.d <- filter(grid, survey_abbrev == "SYN QCS")
write_tex(round((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostQCS")
.d <- filter(grid, survey_abbrev == "SYN WCHG")
write_tex(round((
  sum(.d$restricted)/length(.d$restricted)
) , 2), "lostWCHG")

