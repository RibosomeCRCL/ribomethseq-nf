plot_heatmap_corr <- function(count_matrix = NULL, use_triangle = FALSE) {
  count_matrix["position"] <- NULL
  corr_matrix <- 1 - cor(count_matrix)
  dist_cor <- as.dist(corr_matrix)

  white_red <- colorRampPalette(c("white", "red"), interpolate = "linear")(100)
  if (use_triangle) {
    corr_matrix[lower.tri(corr_matrix)] <- NA

    htmap <- pheatmap::pheatmap(corr_matrix,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      color = white_red,
      breaks = seq(0, 1, by = 0.01),
      main = "correlation-based distance heatmap"
    )
  } else {
    htmap <- pheatmap::pheatmap(corr_matrix,
      clustering_method = "ward.D2",
      cluster_rows = FALSE,
      clustering_distance_cols = dist_cor,
      cutree_cols = 4,
      color = white_red,
      breaks = seq(0, 1, by = 0.01),
      main = "correlation-based distance heatmap"
    )
  }
  return(htmap)
}
