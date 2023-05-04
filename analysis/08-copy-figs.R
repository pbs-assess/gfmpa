f <-
  c(
    "slopes-wchg2.pdf",
    "geo-design-mare-trend.pdf",
    "geostat-vs-design2.pdf",
    "metrics-dotplot-design-geo.pdf",
    "metrics-dotplot-extrapolate.pdf",
    "metrics-dotplot-main.pdf",
    "grid-strata-restricted.pdf",
    "prop-mpa-vs-metrics-wide.pdf",
    "prop-mpa-vs-metrics.pdf",
    "cv-status-quo-mare.pdf",
    "sampled-dotplot-comparison.pdf",
    "index-geo-restricted-highlights.pdf",
    "upsample-example.pdf",
    "fig1.png",
    "fig1.pdf",
    "ts-qcs.pdf",
    "ts-wchg.pdf",
    "ts-hbll.pdf",
    "ts-hs.pdf",
    "ts-qcs-hs.pdf",
    "prop-mpa-vs-metrics-design.pdf",
    "metrics-dotplot-main-design.pdf",
    "sampled-dotplot-comparison2.pdf",
    "metrics-dotplot-by-model.pdf",
    "metrics-cross-plot1.pdf",
    "prop-mpa-vs-metrics-design2.pdf"
  )

# purrr::walk(f[length(f)], function(x) {
purrr::walk(f, function(x) {
  file.copy(file.path("figs", x), "../gf-mpa-index-ms/figs/", overwrite = TRUE)
})
