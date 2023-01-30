#!/usr/bin/env Rscript

args <- commandArgs(trailing = TRUE)

counts_5_end <- as.character(args[1])
counts_3_end <- as.character(args[2])
   fasta_58S <- as.character(args[3])
   fasta_18S <- as.character(args[4])
   fasta_28S <- as.character(args[5])
    fasta_5S <- as.character(args[6])
  out_prefix <- as.character(args[7])

# Read fasta sequence file and
# convert it to a data.frame
fasta2df <- function(file = NULL) {
  if(is.null(file))
    stop("No sequence file provided")

  sequence <- readLines(file, warn = FALSE)
  sequence <- as.character(paste(sequence[-1], collapse = ""))
  sequence <- unlist(strsplit(sequence, split = ""))
  sequence_df <- data.frame(position = seq(1,length(sequence)), sequence = sequence)
  return(sequence_df)
}

fasta_df_58S <- fasta2df(file = fasta_58S)
fasta_df_18S <- fasta2df(file = fasta_18S)
fasta_df_28S <- fasta2df(file = fasta_28S)
fasta_df_5S  <- fasta2df(file = fasta_5S)

# Read counts files
counts5end <- read.csv(counts_5_end, sep = "\t", header = FALSE)
colnames(counts5end) <- c("sequence_info", "position", "counts_5end")

counts3end <- read.csv(counts_3_end, sep = "\t", header = FALSE)
colnames(counts3end) <- c("sequence_info", "position", "counts_3end")
# 3'end counts are shifted to +1 position in order that the gap of counts == to gap of 5'end
counts3end$position <- counts3end$position + 1

# 5.8S rRNA counts
cnt53_58S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[1]),],
                   counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[1]),c(2,3)],
                   by = "position", all = TRUE)
cnt53_58S$counts <- rowSums(cnt53_58S[,c(3,4)], na.rm = TRUE)
cnt53_58S <- merge(cnt53_58S, fasta_df_58S, by = "position", all = TRUE)

# 18S rRNA counts
cnt53_18S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[2]),],
                   counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[2]),c(2,3)],
                   by = "position", all = TRUE)
cnt53_18S$counts <- rowSums(cnt53_18S[,c(3,4)], na.rm = TRUE)
cnt53_18S <- merge(cnt53_18S, fasta_df_18S, by = "position", all = TRUE)

# 28S rRNA counts
cnt53_28S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[3]),],
                   counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[3]),c(2,3)],
                   by = "position", all = TRUE)
cnt53_28S$counts <- rowSums(cnt53_28S[,c(3,4)], na.rm = TRUE)
cnt53_28S <- merge(cnt53_28S, fasta_df_28S, by = "position", all = TRUE)

# 5S rRNA counts
cnt53_5S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[4]),],
                  counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[4]),c(2,3)],
                  by = "position", all = TRUE)
cnt53_5S$counts <- rowSums(cnt53_5S[,c(3,4)], na.rm = TRUE)
cnt53_5S <- merge(cnt53_5S, fasta_df_5S, by = "position", all = TRUE)

#
write.table(cnt53_58S, paste0(out_prefix, ".58S.csv"), col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_18S, paste0(out_prefix, ".18S.csv"), col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_28S, paste0(out_prefix, ".28S.csv"), col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_5S,  paste0(out_prefix, ".5S.csv"),  col.names=T, row.names=F, sep=",", quote=F)
