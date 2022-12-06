#!/usr/bin/env Rscript

get_option <- function(args, opt, dft=NA) {
  ok <- grepl(paste0("^", opt, "="), args)
  if (any(ok)) strsplit(head(args[ok], 1), "=")[[1]][2]
  else dft
}
script_path <- normalizePath(get_option(commandArgs(trailingOnly=F), "--file"))
script_dir  <- dirname(script_path)

args <- commandArgs(trailing = TRUE)
datadir <- as.character(args[1])
 outdir <- as.character(args[2])

if (is.na(outdir)) outdir <- datadir

message(paste("installation directory :", script_dir))
message(paste("data directory :", normalizePath(datadir)))
message(paste("output directory :", normalizePath(outdir)))
message(paste("working directory :", getwd()))

source(file.path(script_dir, "trimmomatic_parser.R"))
source(file.path(script_dir, "report_helper.R"))
source(file.path(script_dir, "r_flagging.R"))

library(ngsReports)
library(reshape2)
library(rmdformats)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)

n_decimals <- 2 # number of decimals

###################### SELECTED COLUMNS FOR THE FINAL TABLE ######################

bowtie_selected_columns <- c("Filename", "Unique_Unpaired", "%Unique_Unpaired", "Multiple_Unpaired", "%Multiple_Unpaired", "Alignment_Rate", "Not_Aligned", "%Not_Aligned")
trimmomatic_selected_columns <- c("Filename", "Dropped", "%Dropped", "Surviving", "%Surviving")
##################################################################################

# First, we import FastQC logs with NgsReport and extract the Basic stats.
fdl <- list.files(path = datadir, pattern = "fastqc.zip", full.names = TRUE, recursive = TRUE)

# import FastQC zip files
base_table <- getModule(fdl, "Basic_Statistics")
deduplicated_percentage <- getModule(fdl, "Total_Deduplicated_Percentage")
deduplicated_percentage["Total_Deduplicated_Percentage"] <- round(deduplicated_percentage["Total_Deduplicated_Percentage"], n_decimals)

base_table <- merge(base_table, deduplicated_percentage, by = "Filename")

base_table["Duplicated_Reads"] <- as.integer(unlist(((100 - base_table["Total_Deduplicated_Percentage"]) / 100) * base_table["Total_Sequences"]))
base_table["Unique_Reads"] <- base_table["Total_Sequences"] - base_table["Duplicated_Reads"]
base_table["%Duplicated_Reads"] <- round(base_table["Duplicated_Reads"] / base_table["Total_Sequences"] * 100, n_decimals)
base_table["%Unique_Reads"] <- round(base_table["Unique_Reads"] / base_table["Total_Sequences"] * 100, n_decimals)

# Adapter content
base_table["%Adapter_content"] <- NA
adapter_content <- getModule(fdl, "Adapter_Content")
write.csv(adapter_content, file = paste(datadir, "/adapter.csv", sep = ""))
adapter_content["Filename"] <- as.factor(adapter_content$Filename)
for (sample_adapter in levels(adapter_content$Filename)) {
  adapter_by_sample <- adapter_content[which(adapter_content["Filename"] == sample_adapter), ]
  adapter_by_sample <- tail(adapter_by_sample, 1)
  base_table["%Adapter_content"][base_table["Filename"] == sample_adapter] <- round(sum(adapter_by_sample[, 3:7]), n_decimals)
}



sequence_quality_scores <- getModule(fdl, "Per_sequence_quality_scores")
base_table["Per_sequence_quality_scores"] <- NA
for (sample_sq in levels(as.factor(sequence_quality_scores$Filename))) {
  sq_bySample <- sequence_quality_scores[which(sequence_quality_scores["Filename"] == sample_sq), ]

  quality_maxCount <- sq_bySample$Quality[which(sq_bySample$Count == max(sq_bySample$Count))]

  base_table["Per_sequence_quality_scores"][base_table["Filename"] == sample_sq] <- quality_maxCount
}

############ TRIMMOMATIC ############

# Then, we import Trimmomatic logs with our modified parser in trimmomatic_parser.R
trim_logs_list <- list.files(path = datadir, pattern = "*.trimmomatic.stats.log", full.names = T, recursive = T)

data_trimmomatic <- suppressWarnings(lapply(trim_logs_list, readLines)) # load all lines for all trimmomatic logs
names(data_trimmomatic) <- basename(trim_logs_list)

# The parser can easily fail if the logs have not a correct structure
trim_logs <- tryCatch(
  {
    trim_logs <- parseTrimmomaticLogs(data_trimmomatic)
    # Inside trim_logs, we have user-specified parameters used by trimmomatic. We will keep them in trim_parameters
    trim_parameters <- t(trim_logs[1, -c(1:9)])
    trim_parameters[is.na(trim_parameters)] <- "Not specified"

    print(paste("Trimmomatic logs have been successfully imported, number of logs :",length(trim_logs_list)))

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



############ BOWTIE 2 ############
bowtie_logs_list <- list.files(path = datadir, pattern = "*.bowtie2.stats.log", full.names = T, recursive = T)

bowtie_logs <- importNgsLogs(bowtie_logs_list, type = "bowtie2")
bowtie_logs["Alignment_Rate"] <- round(bowtie_logs["Alignment_Rate"] * 100, n_decimals)
bowtie_logs["%Unique_Unpaired"] <- round(bowtie_logs["Unique_Unpaired"] / bowtie_logs["Total_Reads"] * 100, n_decimals)
bowtie_logs["%Multiple_Unpaired"] <- round(bowtie_logs["Multiple_Unpaired"] / bowtie_logs["Total_Reads"] * 100, n_decimals)
bowtie_logs["%Not_Aligned"] <- round(bowtie_logs["Not_Aligned"] / bowtie_logs["Total_Reads"] * 100, n_decimals) # P of reads not aligned
bowtie_logs <- bowtie_logs[bowtie_selected_columns]
print(paste("Bowtie logs have been successfully imported, number of logs :",length(bowtie_logs_list)))

############ combine data ############

# first, we rename the filename and then with merge the tables
base_table["Filename"] <- lapply(base_table["Filename"], sub, pattern = "[.].*", replacement = "")
bowtie_logs["Filename"] <- lapply(bowtie_logs["Filename"], sub, pattern = "[.].*", replacement = "")

if (length(trim_logs) > 1) {
  trim_logs["Filename"] <- lapply(trim_logs["Filename"], sub, pattern = "[.].*", replacement = "")
  total_table <- merge(base_table, trim_logs, by = "Filename")
}

total_table <- merge(total_table, bowtie_logs, by = "Filename")

total_table["%Used_reads_for_counting"] <- round(total_table["Unique_Unpaired"] / total_table["Total_Sequences"] * 100, n_decimals)

ggplot_total_table <- total_table
colnames(ggplot_total_table) <- gsub("%", "P", colnames(ggplot_total_table))

############ Generate flags ############
total_table_corrected <- total_table
colnames(total_table_corrected) <- make.names(colnames(total_table_corrected))
flags_summary_list <- GenerateFlags(total_table_corrected, file.path(script_dir, "RULES.CSV"))
flags_summary <- as.data.frame(bind_rows(flags_summary_list))
flags_summary_spread <- spread(flags_summary, "Category", "Status")

tbl_test <- GetOutliers(tbl = total_table, sample_col = "Filename", val_col = "Surviving", up_low = "Lower", threshold = 2)
write.csv(tbl_test, "outlier.csv")
############ export data ############
total_table <- merge(total_table, flags_summary_spread, by = "Filename")
write.csv(total_table, file = paste(datadir, "/QCtable.csv", sep = ""))
fdl <- FastqcDataList(fdl)
fqName(fdl) <- BeautifyFilename(fqName(fdl))
overRep2Fasta(fdl, n = 20, path = "overrep.fasta")

rmarkdown::render(file.path(script_dir, "template_report.rmd"), output_file = "pipeline_report.html", output_dir = outdir, intermediates_dir = outdir)
