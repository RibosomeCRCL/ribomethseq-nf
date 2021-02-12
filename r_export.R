#!Rscript
source("trimmomatic_parser.R")
source("report_helper.R")
library(ngsReports)



args <- commandArgs(trailing=TRUE);

input_resultFolder <- as.character(args[1])
output_file <- as.character(args[2])

n_decimals = 2 #number of places for decimals


###################### SELECTED COLUMNS FOR THE FINAL TABLE ###################### 

Bowtie_selected_columns = c("Filename","Unique_Unpaired", "Multiple_Unpaired", "Alignment_Rate","Not_Aligned","PNot_Aligned","PUnique_Unpaired","PMultiple_Unpaired")
trimmomatic_selected_columns = c("Filename","Dropped","PDropped","Surviving","PSurviving")
################################################################################## 

#First, we import FastQC logs with NgsReport and extract the Basic stats.

fdl <- list.files(paste(input_resultFolder,"/fastqc_before",sep=""), pattern = "fastqc.zip$", full.names = TRUE)

#import FastQC zip files
base_table = getModule(fdl, "Basic_Statistics")
deduplicated_percentage = getModule(fdl, "Total_Deduplicated_Percentage" )
deduplicated_percentage$Total_Deduplicated_Percentage = round(deduplicated_percentage$Total_Deduplicated_Percentage, n_decimals)

base_table = merge(base_table,deduplicated_percentage,by = "Filename")

base_table["Duplicated_Reads"] = as.integer(((100 - base_table$Total_Deduplicated_Percentage)/100) * base_table$Total_Sequences)
base_table["Unique_Reads"] = base_table["Total_Sequences"] - base_table["Duplicated_Reads"]
base_table["PDuplicated_Reads"] = base_table["Duplicated_Reads"] / base_table["Total_Sequences"] *100
base_table["PUnique_Reads"] = base_table["Unique_Reads"] / base_table["Total_Sequences"] * 100

############ TRIMMOMATIC ############ 

#Then, we import Trimmomatic logs with our modified parser in trimmomatic_parser.R
trim_logs_list = list.files(paste(input_resultFolder,"/trim_logs",sep = ""), full.names = T)

data_trimmomatic <- suppressWarnings(lapply(trim_logs_list, readLines)) #load all lines for all trimmomatic logs
names(data_trimmomatic) <- basename(trim_logs_list)

#The parser can easily fail if the logs have not a correct structure
trim_logs <- tryCatch( 
  {
    trim_logs <- parseTrimmomaticLogs(data_trimmomatic)
    print("Trimmomatic logs have been successfully imported")
    
    trim_logs <- trim_logs[, colSums(is.na(trim_logs)) != nrow(trim_logs)]
    trim_logs["PDropped"] = trim_logs$Dropped/trim_logs$Input_Reads *100
    trim_logs["PDropped"] = round(trim_logs["PDropped"],n_decimals)
    trim_logs["PSurviving"] = 100 - trim_logs["PDropped"]
    trim_logs["PSurviving"] = round(trim_logs["PSurviving"],n_decimals)
    
    trim_logs = trim_logs[trimmomatic_selected_columns]
    
    

  },
  error=function(cond) {
    message(paste("Trimmomatic parser failed to read logs !", url))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    message(paste("Trimmomatic parser issued a warning", url))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NA)
  },
  finally={
  }
)



############ BOWTIE 2 ############ 
bowtie_logs_list = list.files(paste(input_resultFolder,"/bowtie2_logs", sep = ""),full.names = T)

bowtie_logs = importNgsLogs(bowtie_logs_list, type = "bowtie2")
print("bowtie logs have been successfully imported!")

bowtie_logs["PUnique_Unpaired"] = round(bowtie_logs["Unique_Unpaired"]/ bowtie_logs["Total_Reads"] *100,n_decimals)
bowtie_logs["PMultiple_Unpaired"] = round(bowtie_logs["Multiple_Unpaired"]/ bowtie_logs["Total_Reads"] * 100, n_decimals)
bowtie_logs["PNot_Aligned"] = round(bowtie_logs["Not_Aligned"] / bowtie_logs["Total_Reads"] * 100, n_decimals) # P of reads not aligned

bowtie_logs = bowtie_logs[Bowtie_selected_columns]

print("Bowtie logs have been successfully imported")

############ combine data ############ 

#first, we rename the filename and then with merge the tables
base_table$Filename = sub("_.*", "",  base_table$Filename)
bowtie_logs$Filename = sub("[.].*", "",  bowtie_logs$Filename)

if(!is.na(trim_logs)) {
  trim_logs$Filename = sub("[.].*", "",  trim_logs$Filename)
  total_table = merge(base_table,trim_logs,by = "Filename")
}

total_table = merge(total_table,bowtie_logs,by = "Filename")

############ export data ############ 
write.csv(total_table, file = paste(input_resultFolder,output_file,sep = ""))
rmarkdown::render("template_report.Rmd",output_file = paste(input_resultFolder,"/pipeline_report.html",sep=""))