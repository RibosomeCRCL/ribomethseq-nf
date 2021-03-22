

################# GENERIC FUNCTIONS ################# 

GetOutliers <- function(tbl,sample_col,val_col, up_low) {
  
  #the list cannot be coerced to type "double". So we unlist first
  tbl_val_unlisted = as.numeric(unlist(tbl[val_col]))
  
  #calcultate IQR
  tbl_IQR = IQR(tbl_val_unlisted)
  tbl_Q1 = quantile(tbl_val_unlisted,1/4)
  tbl_Q3 = quantile(tbl_val_unlisted,3/4)
  
  #calcultate upper and lower bounds
  tbl_LB = (tbl_Q1 - 1.5 * tbl_IQR)
  tbl_UB = (tbl_Q3 + 1.5 * tbl_IQR)
  
  #We flag with "FAIL" outliers...
  #...above upper bound
  tbl_outliers_UB_name = paste(val_col,"_outliers_UB",sep = "")
  tbl[tbl_outliers_UB_name] = tbl[val_col] > tbl_UB
  #...Below lower bound 
  tbl_outliers_LB_name = paste(val_col,"_outliers_LB",sep = "")
  tbl[tbl_outliers_LB_name] = tbl[val_col] < tbl_LB
  
  #according to up_low parameter, three different dataframes can be returned :
  # 1) sample + outliers below lower bound ("Lower")
  # 2) sample + outliers above upper bound ("Upper")
  # 3) sample + outliers above upper bound + outliers below lower bound
  col_to_keep = switch (up_low,
                        "Upper" = tbl_outliers_UB_name,
                        "Lower" = tbl_outliers_LB_name,
                        "Both"  = c(tbl_outliers_LB_name,tbl_outliers_UB_name)
  )
 
  tbl = tbl[c("Filename",col_to_keep)]
  
  return(tbl)
}

FlagOutliers <- function(tbl,sample_col,val_col, up_low) {
# The sole purpose of this function is to replace "TRUE" and "FALSE" from GetOutliers to "FAIL" and "PASS".
  
tbl = GetOutliers(tbl,sample_col, val_col,up_low)

tbl[tbl=="TRUE"] = "WARN"
tbl[tbl=="FALSE"] = "PASS"
 
  return(tbl)

}


################# FLAGS ################# 

FlagDroppedNotAligned <-  function(tbl, fail.trimmomatic = 25, fail.bowtie = 25) {
  # Rule :
  # WARN => the sample has an outlying value for dropped/non-aligned reads but it does not exceed fail.trimmomatic or fail.bowtie
  # FAIL => the percentage of dropped and/or non-aligned exceed fail.trimmomatic or fail.bowtie
  
  #first, let's get outliers
  tbl_outlier_dropped = FlagOutliers(tbl,"Filename","%Dropped", "Upper")

  tbl_outlier_Naligned = FlagOutliers(tbl,"Filename","%Not_Aligned", "Upper") 
  
  #let's get the samples that have reached the "FAIL" threshold
  tbl["fail_dropped"] = tbl["%Dropped"] > fail.trimmomatic
  tbl_fail_dropped = tbl[c("Filename","fail_dropped")]
  
  
  tbl["fail_Naligned"]= tbl["%Not_Aligned"] > fail.bowtie 
  tbl_fail_Naligned = tbl[c("Filename","fail_Naligned")]
  
  #let's merge our outliers and our fails together
  
  tbl_outlier = merge(tbl_outlier_dropped,tbl_outlier_Naligned,by = "Filename")
  tbl_outlier = merge(tbl_outlier,tbl_fail_dropped, by = "Filename")
  tbl_outlier = merge(tbl_outlier,tbl_fail_Naligned, by = "Filename")
  
  #For each sample, we first assign the value from the outlier check ("PASS" or "WARN")
  tbl_outlier["Dropped reads"] = tbl_outlier["%Dropped_outliers_UB"]
  tbl_outlier["Non-Aligned reads"] = tbl_outlier["%Not_Aligned_outliers_UB"]
  #then we check if the sample has reached the "fail" threshold
  tbl_outlier["Dropped reads"][tbl_outlier["fail_dropped"] == TRUE] = "FAIL"
  tbl_outlier["Non-Aligned reads"][tbl_outlier["fail_Naligned"] == TRUE] = "FAIL"
  
  #we only need the "dropped reads" and "Non-Aligned reads" cols, which contain our final results
  tbl_outlier = tbl_outlier[c("Filename","Dropped reads","Non-Aligned reads")]
  
  
  
  #The filename is not a category
  Category = colnames(tbl_outlier)[-1]
  tbl_outlier_tot = tidyr::gather(tbl_outlier, key = Category, value ="Status", one_of(Category))
  
  return(tbl_outlier_tot)
  #A dataframe is returned. 
  # Example output : 
  # Filename    Category          Status
  # 1         Dropped reads       FAIL
  # 10        Dropped reads       PASS
  # 12        Dropped reads       WARN
  # 2         Dropped reads       FAIL
  # 3         Dropped reads       PASS
  # 7         Dropped reads       PASS
  # ...           ...             ...
  # RNA-4   Non-Aligned reads     PASS
  
}



