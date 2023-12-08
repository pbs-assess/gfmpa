old_folder <- "~/src/gfmpa/gf-mpa-index-first-submission"
new_folder <- "~/src/gf-mpa-index-ms"

wd <- getwd()
setwd(new_folder)

f2 <- readLines("main.tex")
for (i in seq_along(f2)) {
  x <- f2[i]
  x <- gsub("\\\\R\\{[a-zA-Z0-9]+\\}", "", x)
  f2[i] <- x
}

# remove response letter at end:
resp <- grep("Response to reviewer comments", f2)
f2 <- c(f2[seq(1, resp-1)], "\\end{document}")

writeLines(f2, paste0(new_folder, "/main-clean.tex"))

system(paste0("latexdiff ", old_folder, "/main.tex main-clean.tex > diff.tex"))

system("pdflatex diff.tex")
system("bibtex diff")
system("pdflatex diff.tex")
system("pdflatex diff.tex")
system("pdflatex diff.tex")
system("pdflatex diff.tex")

system("pdflatex main-clean.tex")
system("bibtex main-clean")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")
system("pdflatex main-clean.tex")

system("open .")

setwd(wd)
