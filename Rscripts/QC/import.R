
GenerateMergedTable <- function(fastQC_name, bowtie_name, trimmomatic_name) {
  
  #Import FastQC files
  
  #Import trimmomatic
  trim_logs_list <- list.files(path = input_datadir, pattern = "*.trimmomatic.stats.log", full.names = T, recursive = T)
  data_trimmomatic <- suppressWarnings(lapply(trim_logs_list, readLines)) # load all lines for all trimmomatic logs
  names(data_trimmomatic) <- basename(trim_logs_list)
  
  trim_logs <- tryCatch({
    trim_logs <- parseTrimmomaticLogs(data_trimmomatic)
    # Inside trim_logs, we have user-specified parameters used by trimmomatic. We will keep them in trim_parameters
    trim_parameters <- t(trim_logs[1, -c(1:9)])
    trim_parameters[is.na(trim_parameters)] <- "Not specified"
    
    print("Trimmomatic logs have been successfully imported")
    
    trim_logs <- trim_logs[, colSums(is.na(trim_logs)) != nrow(trim_logs)]
    trim_logs["%Dropped"] <- trim_logs["Dropped"] / trim_logs["Input_Reads"] * 100
    trim_logs["%Dropped"] <- round(trim_logs["%Dropped"], n_decimals)
    trim_logs["%Surviving"] <- 100 - trim_logs["%Dropped"]
    trim_logs["%Surviving"] <- round(trim_logs["%Surviving"], n_decimals)
    
    trim_logs <- trim_logs[trimmomatic_selected_columns]
  },
  error = function(cond) {
    message(paste("Trimmomatic parser failed to read logs !"))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning = function(cond) {
    message(paste("Trimmomatic parser issued a warning"))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NA)
  },
  finally = {
  }
  )
  
  #Import bowtie 2 files
  bowtie_logs_list <- list.files(path = input_datadir, pattern = "*.bowtie2.stats.log", full.names = T, recursive = T)
  bowtie_logs <- importNgsLogs(bowtie_logs_list, type = "bowtie2")
  
  
}