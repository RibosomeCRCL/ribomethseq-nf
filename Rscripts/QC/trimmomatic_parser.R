# modified from https://github.com/steveped/ngsReports/blob/master/R/importNgsLogs.R

parseTrimmomaticLogs <- function(data, ...) {
  parseTrimmoSingle <- function(x) {

    ## Initialise values which may or may not be present
    Input_Reads <- Input_Read_Pairs <- Surviving <- Both_Surviving <- NA
    Forward_Only_Surviving <- Reverse_Only_Surviving <- NA

    readType <- gsub(".+(PE|SE).+", "\\1", x[[1]])
    Illumina_Clip <- ifelse(
      grepl("ILLUMINACLIP", x[[2]]),
      gsub(".+ILLUMINACLIP:([^ ]+).+", "\\1", x[[2]]),
      NA
    )
    Leading <- ifelse(
      grepl("LEADING", x[[2]]),
      as.integer(gsub(".+LEADING:([0-9:]+).+", "\\1", x[[2]])),
      NA_integer_
    )
    Trailing <- ifelse(
      grepl("TRAILING", x[[2]]),
      as.integer(gsub(".+TRAILING:([0-9:]+).+", "\\1", x[[2]])),
      NA_integer_
    )
    Crop <- ifelse(
      grepl(" CROP", x[[2]]),
      as.integer(gsub(".+ CROP:([0-9:]+).+", "\\1", x[[2]])),
      NA_integer_
    )
    Head_Crop <- ifelse(
      grepl("HEADCROP", x[[2]]),
      as.integer(gsub(".+HEADCROP:([0-9:]+).+", "\\1", x[[2]])),
      NA_integer_
    )
    Sliding_Window <- ifelse(
      grepl("SLIDINGWINDOW", x[[2]]),
      gsub(".+SLIDINGWINDOW:([0-9:]+).+", "\\1", x[[2]]),
      NA
    )
    Min_Len <- ifelse(
      grepl("MINLEN", x[[2]]),
      gsub(".*MINLEN:([0-9:]+).+", "\\1", x[[2]]),
      NA
    )
    Max_Info <- ifelse(
      grepl("MAXINFO", x[[2]]),
      gsub(".+MAXINFO:([^ ]+).+", "\\1", x[[2]]),
      NA
    )
    Avg_Qual <- ifelse(
      grepl("AVGQUAL", x[[2]]),
      as.integer(gsub(".*AVGQUAL:([0-9]+).*", "\\1", x[[2]])),
      NA_integer_
    )
    # Quality_Encoding <- grep("Quality encoding", x, value = TRUE)
    # Quality_Encoding <-
    # gsub("Quality encoding detected as ", "", Quality_Encoding)

    ## Get the line with the summary values
    valLine <- x[[length(x) - 1]]

    if (readType == "SE") {
      Input_Reads <- gsub("Input Reads: ([0-9]+).+", "\\1", valLine)
      Input_Reads <- as.integer(Input_Reads)

      Surviving <- gsub(".+Surviving: ([0-9]+).+", "\\1", valLine)
      Surviving <- as.integer(Surviving)
    }
    if (readType == "PE") {
      Input_Read_Pairs <-
        gsub("Input Read Pairs: ([0-9]+).+", "\\1", valLine)
      Input_Read_Pairs <- as.integer(Input_Read_Pairs)

      Both_Surviving <-
        gsub(".+Both Surviving: ([0-9]+).+", "\\1", valLine)
      Both_Surviving <- as.integer(Both_Surviving)

      Forward_Only_Surviving <-
        gsub(".+Forward Only Surviving: ([0-9]+).+", "\\1", valLine)
      Forward_Only_Surviving <- as.integer(Forward_Only_Surviving)

      Reverse_Only_Surviving <-
        gsub(".+Reverse Only Surviving: ([0-9]+).+", "\\1", valLine)
      Reverse_Only_Surviving <- as.integer(Reverse_Only_Surviving)
    }
    Dropped <- gsub(".+Dropped: ([0-9]+).+", "\\1", valLine)
    Dropped <- as.integer(Dropped)

    return(tibble(
      Type = readType,
      Input_Reads,
      Input_Read_Pairs,
      Surviving,
      Both_Surviving,
      Forward_Only_Surviving,
      Reverse_Only_Surviving,
      Dropped,
      Illumina_Clip,
      Sliding_Window,
      Max_Info,
      Leading,
      Trailing,
      Crop,
      Head_Crop,
      Min_Len,
      Avg_Qual,
      #    Quality_Encoding # we don't have the quality encoding line in our logs because we specified it in our command
    ))
  }

  out <- lapply(data, parseTrimmoSingle)
  out <- dplyr::bind_rows(out)
  out$Filename <- names(data)

  ## Many of the above values may be missing.
  ## Remove them if so using a quick tidy
  # value <- c() # Avoiding an R CMD check NOTE
  # out <- tidyr::gather(out, "key", "value", -1)
  # out <- dplyr::filter(out, !is.na(value))
  # out <- tidyr::spread(out, "key", "value")

  ## Return the final output
  dplyr::select(
    out,
    "Filename",
    "Type",
    starts_with("Input"),
    contains("Surviving"),
    everything()
  )
}
