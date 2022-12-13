plot_boxplot_samples <- function(matrix,values_col_name,values_to_plot,show_outlier) {
  id_vars <- "sample"
  matrix["position"] <- NULL

  matrix_inv <- as.data.frame(t(matrix))
  matrix_inv <- tibble::rownames_to_column(matrix_inv,"sample")

  matrix_melted <- reshape2::melt(matrix_inv, id.vars = id_vars, value.name = values_col_name)  
  shape_outlier <- NA
  if(show_outlier) shape_outlier <- 19
  

  p <- ggplot2::ggplot(matrix_melted, aes_string(x = "sample", y = values_to_plot)) +
      geom_boxplot(outlier.shape = shape_outlier) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  p <-  p +  geom_hline(yintercept = 2, colour = "blue") 

  return(p)
}
