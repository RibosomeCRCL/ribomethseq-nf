plot_coa <- function(count_matrix,col_for_color = NULL, axis = c(1,2)) {
  coa_calculated <- .compute_coa(count_matrix)
  return(.plot_coa(coa_calculated))
}

.compute_coa <- function(counts = NULL){
  counts["position"] <- NULL
  res_coa <- ade4::dudi.coa(counts[complete.cases(counts), ],
                            scannf = FALSE,
                            nf = 5)
  return(res_coa)
}

.plot_coa <- function(dudi.coa = NULL, axis = c(1,2)) {
  plot.coa <- fviz_ca_col(dudi.coa, 
                          repel = T,
                          col.col = "black",
                          pointsize = 2,  #TODO : reduce size
                          labelsize = 4, 
                          axes = axis,
                          title = paste("Correspondence analysis of the raw counts on all genomic positions (", nrow(dudi.coa$tab),")")) 
  plot.coa <- plot.coa + theme(text = element_text(size = 12)) 
  plot.coa <- plot.coa + labs(color = "black", subtitle = paste(ncol(dudi.coa$tab), "samples"))  # col.by.col returns "e"
  return(plot.coa)  
} 
