# restricted_cols <- RColorBrewer::brewer.pal(4, "Dark2")[-3][c(3, 2, 1)]
# plot(1:3, cex = 8, pch = 19, col = restricted_cols)

# restricted_cols <- RColorBrewer::brewer.pal(3, "Set2")
# names(restricted_cols) <- c("HBLL OUT N", "SYN QCS, SYN HS", "SYN WCHG")

pal <- RColorBrewer::brewer.pal(4, "Set2")
restricted_cols <- c(pal[c(1, 2, 3, 4)], pal[3])

# plot(1:4, col = pal, cex = 3, pch = 19)

names(restricted_cols) <- c(
  "SYN WCHG",
  "HBLL OUT N",
  "SYN HS",
  "SYN QCS",
  "SYN QCS, SYN HS"
)

# plot(1:5, cex = 8, pch = 19, col = restricted_cols)
# text(1:5, 1:5, labels = names(restricted_cols))

ggplot2::theme_set(ggsidekick::theme_sleek() +
  theme(panel.grid = element_line(colour = "grey95"),
    panel.grid.major = element_line(linewidth = rel(0.5)),
    panel.grid.minor = element_line(linewidth = rel(0.25)),
    tagger.panel.tag.text = element_text(colour = "grey30")
    )
  )
