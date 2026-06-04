#!/usr/bin/env Rscript

# Merge 5' and 3' read-end counts per reference RNA and annotate every genomic
# position with its nucleotide. Generalised over an arbitrary number of
# references (rRNA, snRNA, ...) found in the count files -- no per-RNA code.
#
# Usage:   r_refine.R <counts_5p> <counts_3p> <reference_fasta> <out_prefix>
# Outputs: one CSV per reference -> <out_prefix>.<reference>.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: r_refine.R <counts_5p> <counts_3p> <reference_fasta> <out_prefix>")
}

counts_5_end    <- args[1]
counts_3_end    <- args[2]
reference_fasta <- args[3]
out_prefix      <- args[4]

# Read a (multi-)FASTA into a named list of per-position data.frames, keyed by
# the sequence id (first token of the header, i.e. the name bedtools/bowtie use).
read_multifasta <- function(file) {
  lines   <- readLines(file, warn = FALSE)
  headers <- grep("^>", lines)
  if (length(headers) == 0L) stop("No FASTA header found in ", file)

  ids    <- sub("\\s.*$", "", sub("^>\\s*", "", lines[headers]))
  starts <- headers + 1L
  ends   <- c(headers[-1] - 1L, length(lines))

  seqs <- Map(function(s, e) {
    if (s > e) "" else paste(lines[s:e], collapse = "")
  }, starts, ends)

  setNames(
    lapply(seqs, function(seq) {
      chars <- strsplit(seq, "")[[1]]
      data.frame(position = seq_along(chars), sequence = chars,
                 stringsAsFactors = FALSE)
    }),
    ids
  )
}

fasta_by_ref <- read_multifasta(reference_fasta)

counts5 <- read.csv(counts_5_end, sep = "\t", header = FALSE)
colnames(counts5) <- c("sequence_info", "position", "counts_5end")

counts3 <- read.csv(counts_3_end, sep = "\t", header = FALSE)
colnames(counts3) <- c("sequence_info", "position", "counts_3end")
# 3'end counts are shifted by +1 so their gap aligns with the 5'end gap.
counts3$position <- counts3$position + 1L

# One output file per reference present in the 5'end counts.
for (ref in unique(counts5$sequence_info)) {
  c5 <- counts5[counts5$sequence_info == ref, ]
  c3 <- counts3[counts3$sequence_info == ref, c("position", "counts_3end")]

  merged <- merge(c5, c3, by = "position", all = TRUE)
  merged$counts <- rowSums(merged[, c("counts_5end", "counts_3end")], na.rm = TRUE)

  fa <- fasta_by_ref[[ref]]
  if (is.null(fa)) {
    warning(sprintf("Reference '%s' absent from %s; 'sequence' column left empty.",
                    ref, reference_fasta))
  } else {
    merged <- merge(merged, fa, by = "position", all = TRUE)
  }

  write.table(merged, paste0(out_prefix, ".", ref, ".csv"),
              col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
}
