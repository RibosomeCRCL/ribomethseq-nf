# RiboMethSeq Quality Control

This directory contains a small Rmarkdown script. It is run during the pipeline execution.

> **Warning**
> This script has been written for GenomeCov-style output. It expects three columns : RNA, position on the RNA and count for this position.

## Requirements

For a conda environment:
  - xorg-libxrender
  - pandoc
  - r-ggplot2
  - r-pheatmap
  - r-factoextra
  - r-reshape2
  - r-dplyr
  - r-tibble
  - r-rmarkdown
 

## Manual run

To run the script manually, the following command can be used :

`Rscript generate_report.R "./path/to/countfiles/" 1 2 3`

Where
* `"./path/to/countfiles/"` corresponds to the path to the directory where the count CSV files are saved.
* `1` corresponds to the column index for the **RNA names** in the count files.
* `2` corresponds to the column index for the **base position** in the count files.
* `3` corresponds to the column index for the **read end count** in the count files.