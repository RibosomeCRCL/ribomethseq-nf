get_option <- function(args, opt, dft=NA) {
  ok <- grepl(paste0("^", opt, "="), args)
  if (any(ok)) strsplit(head(args[ok], 1), "=")[[1]][2]
  else dft
}

script_dir <- dirname(normalizePath(get_option(commandArgs(trailingOnly=F), "--file")))

###
### init
###
library(ggplot2)
library(dplyr)
library(factoextra)
library(pheatmap)

args <- commandArgs(trailing = TRUE)
datadir <- as.character(args[1])

message(paste("installation directory :", script_dir))
message(paste("data directory :", normalizePath(datadir)))
message(paste("working directory :", getwd()))

source(file.path(script_dir, "plot_boxplot.R"))
source(file.path(script_dir, "plot_RLE.R"))
source(file.path(script_dir, "plot_heatmap.R"))
source(file.path(script_dir, "plot_COA.R"))

###
### load count files
###
fcount_list  <- list.files(datadir, pattern = ".csv", full.names = TRUE)
fcount_data  <- lapply(fcount_list, read.csv, sep = "\t", header = FALSE)
names(fcount_data) <- sub('.5_counts.csv', '', basename(fcount_list))
#fcount_matrix <- merge_counts(fcount_data, rna_col, position_col, count_col)
fcount_matrix <-
  do.call('cbind',
           lapply(names(fcount_data), function(name) {
             tab <- fcount_data[[name]]
             rownames(tab) <- paste0(tab[,1],tab[,2])
             tab <- tab[,3,drop=FALSE]
             colnames(tab) <- name
             tab
           })
        )

###
### render report
###
rmarkdown::render(file.path(script_dir, "template.rmd"), output_file = "rms_report.html", output_dir = getwd(), intermediates_dir = getwd())
