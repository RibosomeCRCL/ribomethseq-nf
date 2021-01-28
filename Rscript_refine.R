#!Rscript

args <- commandArgs(trailing=TRUE);

input_file_counts5 <- as.character(args[1])
input_file_counts3 <- as.character(args[2])
fasta_path_58S <- as.character(args[3])
fasta_path_18S <- as.character(args[4])
fasta_path_28S <- as.character(args[5])
fasta_path_5S <- as.character(args[6])
output_file <- as.character(args[7])

# read sequence file

sequence_df <- function(path_to_seq = NULL){
  if(is.null(path_to_seq)){stop("No sequence path")}
  
  read_seq <- readLines(path_to_seq, warn = FALSE)
  sequence <- as.character(paste(read_seq[-1], collapse = ""))
  sequence <- unlist(strsplit(sequence, split = ""))
  sequence_df <- data.frame(position = seq(1,length(sequence)),sequence = sequence)
  return(sequence_df)
}

sequence_df_58S <- sequence_df(path_to_seq = fasta_path_58S)
sequence_df_18S <- sequence_df(path_to_seq = fasta_path_18S)
sequence_df_28S <- sequence_df(path_to_seq = fasta_path_28S)
sequence_df_5S <- sequence_df(path_to_seq = fasta_path_5S)

# read count file

counts5end <- read.csv(input_file_counts5, sep = "\t", header=F)
colnames(counts5end) <- c("sequence_info", "position", "counts_5end")

counts3end <- read.csv(input_file_counts3, sep = "\t", header=F)
colnames(counts3end) <- c("sequence_info", "position", "counts_3end")

# counts of 3'end are shifted to +1 position in order that the gap of counts == to gap of 5'end
counts3end$position <- counts3end$position + 1

# 5.8S rRNA counts

cnt53_58S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[1]),], counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[1]),c(2,3)], 
                   by = "position", all = T)

cnt53_58S$counts <- rowSums(cnt53_58S[,c(3,4)], na.rm = T)
cnt53_58S <- merge(cnt53_58S, sequence_df_58S, by = "position", all = T)


# 18S rRNA counts

cnt53_18S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[2]),], counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[2]),c(2,3)], 
                   by = "position", all = T)

cnt53_18S$counts <- rowSums(cnt53_18S[,c(3,4)], na.rm = T)
cnt53_18S <- merge(cnt53_18S, sequence_df_18S, by = "position", all = T)


# 28S rRNA counts

cnt53_28S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[3]),], counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[3]),c(2,3)], 
                   by = "position", all = T)

cnt53_28S$counts <- rowSums(cnt53_28S[,c(3,4)], na.rm = T)
cnt53_28S <- merge(cnt53_28S, sequence_df_28S, by = "position", all = T)

# 5S rRNA counts

cnt53_5S <- merge(counts5end[which(counts5end$sequence_info == unique(counts5end$sequence_info)[4]),], counts3end[which(counts3end$sequence_info == unique(counts3end$sequence_info)[4]),c(2,3)], 
                   by = "position", all = T)

cnt53_5S$counts <- rowSums(cnt53_5S[,c(3,4)], na.rm = T)
cnt53_5S <- merge(cnt53_5S, sequence_df_5S, by = "position", all = T)

write.table(cnt53_58S, paste("Treatment_5.8S_Sample_",output_file, ".csv", sep = "") , col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_18S, paste("Treatment_18S_Sample_",output_file, ".csv", sep = "") , col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_28S, paste("Treatment_28S_Sample_",output_file, ".csv", sep = "") , col.names=T, row.names=F, sep=",", quote=F)
write.table(cnt53_5S, paste("Treatment_5S_Sample_",output_file, ".csv", sep = "") , col.names=T, row.names=F, sep=",", quote=F)
