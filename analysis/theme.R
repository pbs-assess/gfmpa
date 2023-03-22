# restricted_cols <- RColorBrewer::brewer.pal(4, "Dark2")[-3][c(3, 2, 1)]
# plot(1:3, cex = 8, pch = 19, col = restricted_cols)

restricted_cols <- RColorBrewer::brewer.pal(4, "Set2")[c(2, 3, 1)]
# plot(1:3, cex = 8, pch = 19, col = restricted_cols)


ggplot2::theme_set(ggsidekick::theme_sleek()) #+
    # theme(panel.grid = element_line(colour = "grey95"),
      # panel.grid.major = element_line(linewidth = rel(0.5)),
      # panel.grid.minor = element_line(linewidth = rel(0.25))))
