library(tidyverse)
library(optparse)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)

option_list = list(
  make_option(c("--output_file"), action="store", help = "path to kraken output file"),
  make_option(c("--kraken_report"), action="store", help = "path to kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to standard kraken mpa report"),
  make_option(c("--extract_file"), action="store", help = "extract filename path"),
  make_option(c("--keep_original"), action="store", default=T, help ="delete original fastq file? T/F"),
  make_option(c("--ntaxid"), action="store", default = 8000, help = "number of taxids to extract at a time"),
  make_option(c("--cores"), action="store", default = 8, help = "number of cores at a time")
)

## ntaxid: There is a certain limit on the length of command line in Linux system. The maximum length of the command line is determined by the system kernel parameter `ARG_MAX`, 
## which is usually several thousand or tens of thousands of characters. If you enter more characters than this limit on the command line, the system will not be able 
## to handle the command properly. Therefore, when searching for a large number of tax IDs, it is important to avoid exceeding the limit of command line length for 
## successful execution of the command.

opt = parse_args(OptionParser(option_list = option_list))

kr = read.delim(opt$kraken_report, header = F)
# removing root and unclassified taxa
kr = kr[-c(1:2), ]
mpa = read.delim(opt$mpa_report, header = F)

n = str_which(mpa$V1, 'k__Bacteria|k__Fungi|k__Viruses')
taxid = kr$V7[n]

taxid.list = split(taxid, ceiling(seq_along(taxid)/opt$ntaxid))



if(file.exists(opt$extract_file)){
  system(paste0('rm ', opt$extract_file))
}

# Here, we construct the output file by capturing the system output, so that Snakemake can receive the output thrown by this script.
n_cores <- opt$cores
registerDoParallel(cores = n_cores)

results <- foreach(i = 1:length(taxid.list), .packages = c("stringr"), .combine = "c") %dopar% {
  cat("Extracting output data", i, "/", length(taxid.list), "\n")
  
  taxid = paste0("(taxid ", taxid.list[[i]], ")", collapse = "\\|")
  taxid = paste0("'", taxid, "'")
  str = paste0("grep -w ", taxid, " ", opt$output_file)
  output = system(str, intern = TRUE)
  return(output)
}


cat(results, file = opt$extract_file, sep = "\n")

# results <- foreach(i = 1:length(taxid.list), .packages = c("stringr")) %dopar% {
#   cat("Extracting output data", i, "/", length(taxid.list), "\n")
  
#   taxid = paste0("(taxid ", taxid.list[[i]], ")", collapse = "\\|")
#   taxid = paste0("'", taxid, "'")
#   str = paste0("grep -w ", taxid, " ", opt$output_file, " >> ", opt$extract_file)
#   result <- try(system(paste0(str, " 2>&1"), intern = TRUE, ignore.stderr = FALSE), silent = TRUE)

#   if (inherits(result, "try-error")) {
#     message("Command execution failed with error message: ", result)
#   } else {
#     # 命令执行成功，打印输出结果
#     cat("ok")
#   }
# }

# ## TODO:use foreach package to do parrael run
# for(i in 1:length(taxid.list)){
#   print(paste('Extracting output data', i, '/', length(taxid.list)))
  
#   taxid = paste0("(taxid ", taxid.list[[i]], ")", collapse = "\\|")
#   taxid = paste0("'", taxid, "'")
#   str = paste0("grep -w ", taxid, " ", opt$output_file, " >> ", opt$extract_file)
#   result <- try(system(paste0(str, " 2>&1"), intern = TRUE, ignore.stderr = FALSE), silent = TRUE)

#   if (inherits(result, "try-error")) {
#     message("Command execution failed with error message: ", result)
#   } else {
#     # 命令执行成功，打印输出结果
#     cat("ok")
#   }
# }

if(opt$keep_original == F){
  system(paste('rm', opt$output_file))
}

print('Done')