# modifies the code
library(d3heatmap)
myd3 <- function (x, cluster = !any(is.na(x)), theme = NULL, colors = "RdYlBu", 
                  width = NULL, height = NULL, xaxis_height = 120, yaxis_width = 120, 
                  xaxis_font_size = NULL, yaxis_font_size = NULL, brush_color = "#0000FF", 
                  show_grid = TRUE, anim_duration = 500, heatmap_options = list()) 
{
  matrix <- as.matrix(x)
  rng <- range(matrix, na.rm = TRUE)
  rowDend <- NULL
  colDend <- NULL
  options <- NULL
  tmp <- tempfile()
  png(tmp)
  heatmap_options$x = quote(matrix)
  heatmap_options$keep.dendro = TRUE
  hm <- do.call(stats::heatmap, heatmap_options)
  dev.off()
  unlink(tmp)
  matrix <- matrix[hm$rowInd, hm$colInd]
  options <- c(options, list(xaxis_height = xaxis_height, yaxis_width = yaxis_width, 
                             xaxis_font_size = xaxis_font_size, yaxis_font_size = yaxis_font_size, 
                             brush_color = brush_color, show_grid = show_grid, anim_duration = anim_duration))
  domain <- seq.int(rng[1], rng[2], length.out = 100)
  colors <- (leaflet::colorNumeric(colors, 1:100))(1:100)
  matrix <- list(data = as.numeric(t(matrix)), dim = dim(matrix), 
                 rows = row.names(matrix) %||% paste(1:nrow(matrix)), 
                 cols = colnames(matrix) %||% paste(1:ncol(matrix)), colors = colors, 
                 domain = domain)
  x <- list(rows = rowDend, cols = colDend, matrix = matrix, 
            theme = theme, options = options, hm = hm)
  htmlwidgets::createWidget(name = "d3heatmap", x, width = width, 
                            height = height, package = "d3heatmap", sizingPolicy = htmlwidgets::sizingPolicy(browser.fill = TRUE))
}


tmpfun <- get("d3heatmap", envir = asNamespace("d3heatmap"))
environment(myd3) <- environment(tmpfun)
attributes(myd3) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("d3heatmap", myd3, ns="d3heatmap")