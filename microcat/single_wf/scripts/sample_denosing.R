library(optparse)
library(data.table)
library(tibble)
library(ggplot2)
library(dplyr)
# library(tidyverse)
library(scales)

read_kraken_reports = function(files, sample_names = NULL, study_name = NULL, min_reads = 2, min_uniq = 2,path){
  
  if(is.null(sample_names)){sample_names = files}
  if(is.null(study_name)){study_name = NA}
  if(length(study_name) == 1){study_name = rep(study_name, length(files))}
  
  df = list()
  n = 0
  for(i in 1:length(files)){
    if(round(i/length(files)*100, 2) > n){n = round(i/length(files)*100, 2); cat(paste0('\r',n,'% done   '))}
    x = read.delim(paste0(path,"/",files[i]), header = F)
    x$V8 = trimws(x$V8)
    total_reads = x$V2[1] + x$V2[2]
    n_microbiome_reads = sum(x$V2[x$V8 %in% c('Bacteria', 'Fungi', 'Viruses')])
    df[[i]] = data.frame(study = study_name[i], sample = sample_names[i],
                         rank = x$V6, taxid = x$V7, name = x$V8, 
                         reads = x$V2, min = x$V4, uniq = x$V5, 
                         rpm = x$V2/total_reads*10^6,
                         rpmm = x$V2/n_microbiome_reads*10^6)
  }
  df = rbindlist(df) %>% tibble()
  cat('\n')
  
  df = subset(df, reads >= min_reads & uniq >= min_uniq) 
  df
}

# 创建选项列表
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default = NULL,
              help = "One or more file path containing custom style Kraken reports"),
  make_option(c("-m", "--kmer_data"), type = "character", default = NULL,
              help = "single kmer report path"),
  make_option(c("-o", "--out_path"), type = "character", default = NULL,
              help = "Result output path"),
  make_option(c("-s", "--sample_names"), type = "character", default = NULL,
              help = "One or more sample names corresponding to the input files"),
  make_option(c("-n", "--study_name"), type = "character", default = NULL,
              help = "Name of the study"),
  make_option(c("-r", "--min_reads"), type = "integer", default = 2,
              help = "Minimum number of reads per taxon"),
  make_option(c("-u", "--min_uniq"), type = "integer", default = 2,
              help = "Minimum number of unique sequences per taxon")           
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


files <- list.files(opt$path,recursive = TRUE)
stopifnot(length(files) > 0) # make sure the input dir path is set correctly
min_reads <- opt$min_reads
min_uniq <- opt$min_uniq

kraken_report_sample <- read_kraken_reports(files = files,
                                            path = opt$path,
                                            min_reads = min_reads, 
                                            min_uniq = min_uniq)

kraken_report_sample_filtered = kraken_report_sample %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  select(-nn)

kraken_report_sample_filtered$sample <- gsub("/.*", "", kraken_report_sample_filtered$sample)

# 假设你的数据赋值为 min 和 uniq
# 先尝试进行 spearman 相关性检验
corr_test <- function(x, y) {
  result <- tryCatch(cor.test(x, y, method = 'spearman'), error = function(e) e)
  if (inherits(result, "error")) {
    return(1)  # 返回1，表示无法计算p值
  } else {
    return(result$p.value)
  }
}
## since some data only have one
# run correlations 
corr_data = kraken_report_sample_filtered %>%
  subset(rank %in% c('G', 'S')) %>% 
  group_by(name) %>%
  summarize(r1 = cor(min,uniq,method='spearman'),
            r2 = cor(min,reads,method='spearman'),
            r3 = cor(reads,uniq,method='spearman'),
            p1 = corr_test(min,uniq),
            p2 = corr_test(min,reads),
            p3 = corr_test(reads,uniq)
            )

kmer_barcode <- read.table(opt$kmer_data,header=T,sep="\t")

kmer_test = left_join(kmer_barcode, corr_data, by = 'name') %>% 
  left_join(select(kraken_report_sample, rank, name) %>% distinct())%>% 
  subset(r1>0 & r2>0 & r3>0 & p_adjust<0.05 & p1<0.05 & p2<0.05 & p3<0.05 & rank == 'S')

cell.lines = read.delim('/data/cell.lines.txt', header = T) %>% tibble()

df = cell.lines[,1:11] %>% mutate(study = 'cell lines'); df = df[, -2]


df = rbind(df, kraken_report_sample)

qtile = 0.99
q_df = cell.lines %>%
  group_by(name, rank) %>% 
  summarize(CLrpmm = 10^quantile(log10(rpmm), qtile, na.rm = T), 
            .groups = 'keep')

kr_rpmm <- subset(kraken_report_sample_filtered, sample == opt$sample_names)


sample_kmer_test <- left_join(kmer_test, q_df, by = c('name', 'rank')) %>% 
  left_join(kr_rpmm, by = c('name','taxid'))



write.table(sample_kmer_test, file = opt$out_path, quote = F, row.names = F,sep="\t")


