merge_counts <- function(count_list, rna_col, position_col,count_col) {
    
  sample_list <- count_list
  sample_list_nm <- names(sample_list)
  
  # Get the positions from on sample and use them to prepare our merged dataframe
  # (the positions are identical to all samples, as they are aligned on same reference)
  position_list <- sample_list[[1]][,c(rna_col,position_col)]
  position_list <- paste0(position_list[,1],"_",position_list[,2])
  matrix_all <- data.frame(position = position_list)

    for(sample_nm in sample_list_nm) {
    sample_df <- sample_list[[sample_nm]]
    position_list <- sample_df[,c(rna_col,position_col)]
    position_list <- paste0(position_list[,1],"_",position_list[,2])
    sample_df["position"] <- position_list
    count_col <- 3
    sample_df <- sample_df[,c(4,count_col)]
    matrix_all <- dplyr::full_join(matrix_all,sample_df,by="position")
    names(matrix_all)[length(names(matrix_all))] <- sample_nm 
    matrix_all <- matrix_all[match(position_list,matrix_all[,"position"]),]
    }

    return(matrix_all)
  
  }
