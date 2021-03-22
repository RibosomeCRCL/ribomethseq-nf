
plotReadTotalsPercentage <- function(table, col_sample, col_uniqueReads, col_duplicatedReads, sample.order) {
  #' Show the percentage of unique and duplicated reads in a ggplo2 barplot
  #' @param table a table containing data for all samples
  #' @param col_sample name of the column containing samples' name
  #' @param col_uniqueReads name of the column containing the percentage of unique reads by sample
  #' @param col_duplicatedReads name of the column containing the percentage of duplicated reads by sample

  ggplot_table = table[c(col_sample,col_uniqueReads,col_duplicatedReads)]
  colnames(ggplot_table) = c("Filename","Unique Reads","Duplicated Reads")
  Type = c("Unique Reads", "Duplicated Reads")
  ggplot_table = tidyr::gather(ggplot_table, key = Type, value ="Total", one_of(Type))
  
  p = ggplot(data=ggplot_table, aes_string(x=sample.order, y="Total", fill="Type", label = "Total"))+
    geom_bar( stat="identity") + ylab("percentage of reads") + xlab("Sample") 
  p = p+coord_flip()

  return(p)
}

plotBowtieReport <- function(bowtie_table, percentage = F, sample.order) {
  
  Type = c("Unique Unpaired","Multiple Unpaired","Not Aligned")
  
  #renaming some columns
  colnames(bowtie_table) <- gsub("_", " ", colnames(bowtie_table))

  bowtie_table = bowtie_table[c("Filename",Type)]
  

  #reformat the table to make it readable for ggplot2
  bowtie_table = tidyr::gather(bowtie_table, key = Type, value ="Total", one_of(Type))
  
  if(!percentage) {
  p = ggplot(data=bowtie_table, aes_string(x="Filename", y="Total", fill="Type", label = "Total"))+
    geom_bar( stat="identity") + ylab("Number of reads") + xlab("Sample") + xlim(levels(sample.order))
  }
  else {
    p = ggplot(data=bowtie_table, aes_string(x="Filename", y="Total", fill="Type", label = "Total"))+
      geom_bar( position = "fill", stat="identity")+ ylab("Percentage of reads") + xlab("Sample") + xlim(levels(sample.order))
    
  }
  p = p+coord_flip()
  
  return(p)
  
}

plotTrimmomaticDropped <- function(table, percentage = F, sample.order) {
  Type = c("Surviving","Dropped")
  
  trim_table = table[c("Filename",Type)]
  
  trim_table = tidyr::gather(trim_table, key = Type, value ="Total", one_of(Type))
  
  if(!percentage) {
    p = ggplot(data=trim_table, aes_string(x=sample.order, y="Total", fill="Type", label = "Total"))+
      geom_bar( stat="identity") + ylab("Number of reads") + xlab("Sample") 
  }
  else {
    p = ggplot(data=trim_table, aes_string(x=sample.order, y="Total", fill="Type", label = "Total"))+
      geom_bar( position = "fill", stat="identity")+ ylab("Percentage of reads") + xlab("Sample") 
    
  }
  p = p+coord_flip()
  
  return(p)
  
  
}

plotTrimBowtie <- function(d,sample.order) {
  message("trimbowtie")

  
 p = ggplot(d, aes(sample.order,value, col=variable)) + 
    geom_point() + 
    stat_smooth() + coord_flip() + xlab("Sample") + ylab("% of reads") + ylim(0,100)
  
  return(p)
}

plotSummaryNewFlags <- function(table,fdl,sample.order) {

tbl_outlier_tot = FlagDroppedNotAligned(table)

  #Get FastQC's summary and shorten the filenames
  fastqc_summary = getSummary(fdl)
  fastqc_summary = fastqc_summary[c("Filename","Category","Status")] 
  fastqc_summary["Filename"]= lapply(fastqc_summary["Filename"],BeautifyFilename)
  
    tbl_outlier_tot = rbind(tbl_outlier_tot, fastqc_summary)
  
    # tbl_outlier_tot = tbl_outlier_tot[order(match(tbl_outlier_tot["Filename"],sample.order))]
    # print(head(tbl_outlier_tot))
  #plot our new summary
  sumPlot <- ggplot(tbl_outlier_tot, aes_string("Category", "Filename", fill = "Status")) +
    geom_tile(colour = "black", size = 0.5) +
    scale_fill_manual(breaks = c("FAIL","WARN","PASS"),values = c("coral","yellow","palegreen")) +
    labs(x = "QC Category", y = "Sample") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    xlim(c("Adapter Content","Per base sequence content","Per base sequence quality","Per sequence quality scores", "Per base N content","Sequence Length Distribution","Dropped reads","Non-Aligned reads")) +
    ylim(levels(sample.order)) +
    theme_bw() +
    theme(axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5
    ))

  return(sumPlot)  
  
}

BeautifyFilename <- function(name) {
  
  return(sub("_.*", "\\1", name, perl = T))

}
