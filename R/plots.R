# Plotting helpers.

library(ggplot2)

theme_pub <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey60"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

# Black or white label depending on the luminance of the fill colour.
contrast_text <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum > 0.45, "black", "white")
}

fill_colours <- function(x, limits, option = "viridis") {
  pal <- viridisLite::viridis(256, option = option)
  r <- (pmin(pmax(x, limits[1]), limits[2]) - limits[1]) / diff(limits)
  pal[pmax(1, ceiling(r * 255) + 1)]
}

save_figure <- function(p, file, width, height, dpi = 600) {
  ggsave(paste0(file, ".tiff"), p, width = width, height = height, units = "cm",
         dpi = dpi, compression = "lzw", bg = "white")
  ggsave(paste0(file, ".pdf"), p, width = width, height = height, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(paste0(file, ".png"), p, width = width, height = height, units = "cm",
         dpi = 300, bg = "white")
}
