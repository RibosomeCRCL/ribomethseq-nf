
plotReadTotalsPercentage <- function(table, col_sample, col_uniqueReads, col_duplicatedReads) {
  #' Show the percentage of unique and duplicated reads in a ggplo2 barplot
  #' @param table a table containing data for all samples
  #' @param col_sample name of the column containing samples' name
  #' @param col_uniqueReads name of the column containing the percentage of unique reads by sample
  #' @param col_duplicatedReads name of the column containing the percentage of duplicated reads by sample

  ggplot_table = table[c(col_sample,col_uniqueReads,col_duplicatedReads)]
  Type = c(col_uniqueReads, col_duplicatedReads)
  ggplot_table = tidyr::gather(ggplot_table, key = Type, value ="Total", one_of(Type))
  
  p = ggplot(data=ggplot_table, aes_string(x=col_sample, y="Total", fill="Type", label = "Total"))+
    geom_bar( stat="identity") +   geom_text(position = position_stack(vjust = 0.5))
  p = p+coord_flip()

  return(p)
}

plotBowtieReport <- function(bowtie_table, percentage = F) {
  
  Type = c("Unique_Unpaired","Multiple_Unpaired","Not_Aligned")
  bowtie_table = bowtie_table[c("Filename",Type)]

  bowtie_table = tidyr::gather(bowtie_table, key = Type, value ="Total", one_of(Type))
  
  if(!percentage) {
  p = ggplot(data=bowtie_table, aes_string(x="Filename", y="Total", fill="Type", label = "Total"))+
    geom_bar( stat="identity") 
  }
  else {
    p = ggplot(data=bowtie_table, aes_string(x="Filename", y="Total", fill="Type", label = "Total"))+
      geom_bar( position = "fill", stat="identity") 
    
  }
  p = p+coord_flip()
  
  return(p)
  
}
