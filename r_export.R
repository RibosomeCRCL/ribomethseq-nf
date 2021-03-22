#!Rscript
library(ngsReports)
library(reshape2)
library(rmdformats)

args <- commandArgs(trailing=TRUE);

input_installDir <- as.character(args[1]) # get directory where this script and its dependencies are installed
output_file <- as.character(args[2]) 

source(paste(input_installDir,"/trimmomatic_parser.R",sep=""))
source(paste(input_installDir,"/report_helper.R",sep=""))
source(paste(input_installDir,"/r_flagging.R",sep=""))

n_decimals = 2 #number of decimals

current_wd = getwd() # get the temp directory assigned by nextflow
print(current_wd)

###################### SELECTED COLUMNS FOR THE FINAL TABLE ###################### 

Bowtie_selected_columns = c("Filename","Unique_Unpaired","%Unique_Unpaired", "Multiple_Unpaired","%Multiple_Unpaired", "Alignment_Rate","Not_Aligned","%Not_Aligned")
trimmomatic_selected_columns = c("Filename","Dropped","%Dropped","Surviving","%Surviving")
################################################################################## 

#First, we import FastQC logs with NgsReport and extract the Basic stats.

fdl <- list.files(pattern = "fastqc.zip$", full.names = TRUE)

#import FastQC zip files
base_table = getModule(fdl, "Basic_Statistics")
deduplicated_percentage = getModule(fdl, "Total_Deduplicated_Percentage" )
deduplicated_percentage$Total_Deduplicated_Percentage = round(deduplicated_percentage$Total_Deduplicated_Percentage, n_decimals)

base_table = merge(base_table,deduplicated_percentage,by = "Filename")

base_table["Duplicated_Reads"] = as.integer(((100 - base_table$Total_Deduplicated_Percentage)/100) * base_table$Total_Sequences)
base_table["Unique_Reads"] = base_table["Total_Sequences"] - base_table["Duplicated_Reads"]
base_table["%Duplicated_Reads"] = round(base_table["Duplicated_Reads"] / base_table["Total_Sequences"] *100, n_decimals)
base_table["%Unique_Reads"] = round(base_table["Unique_Reads"] / base_table["Total_Sequences"] * 100, n_decimals)

############ TRIMMOMATIC ############ 

#Then, we import Trimmomatic logs with our modified parser in trimmomatic_parser.R
trim_logs_list = list.files(pattern = "*.trimmomatic.stats.log", full.names = T)

data_trimmomatic <- suppressWarnings(lapply(trim_logs_list, readLines)) #load all lines for all trimmomatic logs
names(data_trimmomatic) <- basename(trim_logs_list)

#The parser can easily fail if the logs have not a correct structure
trim_logs <- tryCatch( 
  {
    trim_logs <- parseTrimmomaticLogs(data_trimmomatic)
    #Inside trim_logs, we have user-specified parameters used by trimmomatic. We will keep them in trim_parameters
    trim_parameters <- t(trim_logs[1,-c(1:9)])
    trim_parameters[is.na(trim_parameters)] <- "Not specified"

    print("Trimmomatic logs have been successfully imported")
    
    trim_logs <- trim_logs[, colSums(is.na(trim_logs)) != nrow(trim_logs)]
    trim_logs["%Dropped"] = trim_logs$Dropped/trim_logs$Input_Reads *100
    trim_logs["%Dropped"] = round(trim_logs["%Dropped"],n_decimals)
    trim_logs["%Surviving"] = 100 - trim_logs["%Dropped"]
    trim_logs["%Surviving"] = round(trim_logs["%Surviving"],n_decimals)
    
    trim_logs = trim_logs[trimmomatic_selected_columns]
    
    

  },
  error=function(cond) {
    message(paste("Trimmomatic parser failed to read logs !"))
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    message(paste("Trimmomatic parser issued a warning"))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NA)
  },
  finally={
  }
)



############ BOWTIE 2 ############ 
bowtie_logs_list = list.files(pattern = "*.bowtie2.stats.log",full.names = T)

bowtie_logs = importNgsLogs(bowtie_logs_list, type = "bowtie2")
print("bowtie logs have been successfully imported!")
bowtie_logs["Alignment_Rate"] = round(bowtie_logs["Alignment_Rate"]*100, n_decimals)
bowtie_logs["%Unique_Unpaired"] = round(bowtie_logs["Unique_Unpaired"]/ bowtie_logs["Total_Reads"] *100,n_decimals)
bowtie_logs["%Multiple_Unpaired"] = round(bowtie_logs["Multiple_Unpaired"]/ bowtie_logs["Total_Reads"] * 100, n_decimals)
bowtie_logs["%Not_Aligned"] = round(bowtie_logs["Not_Aligned"] / bowtie_logs["Total_Reads"] * 100, n_decimals) # P of reads not aligned
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

total_table["%Used_reads_for_counting"] = round((total_table["Total_Sequences"] - total_table["Dropped"] - total_table["Not_Aligned"]) / total_table["Total_Sequences"] * 100, n_decimals)

ggplot_total_table = total_table
colnames(ggplot_total_table) <- gsub("%", "P", colnames(ggplot_total_table))
############ export data ############ 
write.csv(total_table, file = output_file)
fdl = FastqcDataList(fdl)
fqName(fdl) = BeautifyFilename(fqName(fdl))
overRep2Fasta(fdl,n = 20, path = "overrep.fasta")
rmarkdown::render(paste(input_installDir,"/template_report.rmd",sep=""),output_file = paste(current_wd,"/pipeline_report.html",sep=""))