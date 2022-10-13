col_palette <- c(
  "day0" = rgb(177, 177, 177, maxColorValue = 255),
  "day10_CIS" = rgb(235, 194, 47, maxColorValue = 255),
  "day10_COCL2" = rgb(119, 228, 210, maxColorValue = 255),
  "day10_DABTRAM" = rgb(204, 124, 176, maxColorValue = 255),
  "week5_CIS" = rgb(235, 134, 47, maxColorValue = 255),
  "week5_COCL2" = rgb(70, 177, 70, maxColorValue = 255),
  "week5_DABTRAM" = rgb(122, 49, 126, maxColorValue = 255)
)
col_palette <- col_palette[order(names(col_palette))]

# plot(1:length(col_palette), pch = 16, col = col_palette, cex = 5)

