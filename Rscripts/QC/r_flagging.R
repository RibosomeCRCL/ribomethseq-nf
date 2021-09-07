

################# GENERIC FUNCTIONS ################# 

GetOutliers <- function(tbl,sample_col,val_col, up_low, threshold) {
  
  #the list cannot be coerced to type "double". So we unlist first
  tbl_val_unlisted = as.numeric(unlist(tbl[val_col]))
  
  tbl_median = median(tbl_val_unlisted)
  tbl_mad <- threshold*mad(tbl_val_unlisted)
  tbl_UB = tbl_median + tbl_mad
  tbl_LB = tbl_median - tbl_mad
  

  
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
  tbl["is_outlier"] = tbl[tbl_outliers_UB_name] + tbl[tbl_outliers_LB_name] >= 1
  tbl = tbl[c("Filename",col_to_keep,"is_outlier")]
  
  return(tbl)
}

FlagOutliers <- function(tbl,sample_col,val_col, up_low) {
# The sole purpose of this function is to replace "TRUE" and "FALSE" from GetOutliers to "FAIL" and "PASS".
  
tbl = GetOutliers(tbl,sample_col, val_col,up_low)

tbl[tbl=="TRUE"] = "WARN"
tbl[tbl=="FALSE"] = "PASS"
 
  return(tbl)

}


FlagByValue <- function(tbl,flag_name,sample_col,val_col,warning_condition,fail_condition) {
  
  tbl_vals = as.vector(unlist(tbl[val_col]))
  warning_condition = eval(parse(text = paste("function(x) x", warning_condition)))
  fail_condition = eval(parse(text = paste("function(x) x", fail_condition)))
  
  tbl["warning"] = warning_condition(tbl_vals)
  tbl["fail"] = fail_condition(tbl_vals)
  tbl[flag_name] = "PASS"
  tbl[flag_name][tbl["warning"] == T] = "WARN"
  tbl[flag_name][tbl["fail"] == T] = "FAIL"
  
  tbl = tbl[,c("Filename",flag_name)]
  
  Category = flag_name
  tbl_output = tidyr::gather(tbl, key = Category, value ="Status", one_of(Category))
  return(tbl_output)
  
  # Example call : Check if the number of unpaired reads is either under 30000000 (=Warning) or 15000000 (=Fail)
  # FlagByValue(total_table,"Per sequence test","Filename","Unique_Unpaired",function(x) x<30000000, function(x) x<15000000)
  
  # Example output
#  Filename     Category            Status
#        1     Per sequence test     FAIL
#       10     Per sequence test   WARNING
#       12     Per sequence test   WARNING
#        2     Per sequence test     FAIL
}

FlagByMAD <- function(tbl,flag_name,sample_col,val_col,rule_tbl) {
  write.csv(rule_tbl,"rules.csv")
  message(paste("calling flagbymad.", val_col))
  tmp_tbl <- data.frame(filename <-  tbl[sample_col],val_col= tbl[val_col] )

  warning.lower = rule_tbl[which(rule_tbl["QC.name"] == val_col),"warning.lower"]
  warning.upper = rule_tbl[which(rule_tbl["QC.name"] == val_col),"warning.upper"]
  
  fail.lower = rule_tbl[which(rule_tbl["QC.name"] == val_col),"fail.lower"]
  fail.upper = rule_tbl[which(rule_tbl["QC.name"] == val_col),"fail.upper"]
  
  
  
  message("rules loaded")
  tmp_tbl["is.warning"] <- lapply(tmp_tbl[val_col], function(x) x <= as.numeric(warning.lower) | x  >= as.numeric(warning.upper))
  tmp_tbl["is.fail"] = lapply(tmp_tbl[val_col], function(x) x <= as.numeric(fail.lower) | x  >= as.numeric(fail.upper))
  message("tests done")
  tmp_tbl[flag_name] = "PASS"
  tmp_tbl[flag_name][tmp_tbl["is.warning"] == T] = "WARN"
  tmp_tbl[flag_name][tmp_tbl["is.fail"] == T] = "FAIL"

  tmp_tbl <- tmp_tbl[,c("Filename",flag_name)]

  
  Category <- flag_name
  tbl_output <- tidyr::gather(tmp_tbl, key = Category, value ="Status", one_of(Category))

  return(tbl_output)
  
  # Example call : Check if the number of unpaired reads is either under 30000000 (=Warning) or 15000000 (=Fail)
  # FlagByValue(total_table,"Per sequence test","Filename","Unique_Unpaired",function(x) x<30000000, function(x) x<15000000)
  
  # Example output
  #  Filename     Category            Status
  #        1     Per sequence test     FAIL
  #       10     Per sequence test   WARNING
  #       12     Per sequence test   WARNING
  #        2     Per sequence test     FAIL
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

GenerateFlags <- function(tbl, config_table_path) {
  flags_file = read.csv(config_table_path,sep = ",")
  flag_list = list()
  for(row in 1:nrow(flags_file))  {
   flag_list[[row]] = FlagByMAD(tbl,sample_col = "Filename",val_col = flags_file[row,"QC.name"],rule_tbl = flags_file,flag_name = flags_file[row,"QC.title"])
  }
  return(flag_list)
}




