plot_boxplot_samples <- function(matrix,show_outlier) {
  matrix["position"] <- NULL


  matrix_melted <- tidyr::gather(matrix, key = sample, value = counts)
  shape_outlier <- NA
  if(show_outlier) shape_outlier <- 19

  p <- ggplot2::ggplot(matrix_melted, aes_string(x = "sample", y = "log10(counts)")) +
      geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  p <-  p +  geom_hline(yintercept = 2, colour = "blue") 

  return(p)
}