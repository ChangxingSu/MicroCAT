library(optparse)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)


option_list = list(
  make_option(c("--output_file"), action="store", help = "path to kraken output file"),
  make_option(c("--taxa"),type = "character", help = "tsv containing taxa to extract"),
  make_option(c("--out_krak2_denosing"), action="store", help = "output path to save files")
)

opt = parse_args(OptionParser(option_list = option_list))


krak2_output_file = read.delim(opt$output_file, header = F)

krak2_output_file_taxid = krak2_output_file %>% 
  select(-V1) %>% 
  separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
  mutate(taxid = str_remove(taxid, '\\)') %>% trimws(),
        name = trimws(name))


taxa = read.table(opt$taxa,sep="\t",header=T)


krak2_output_file_taxid_denosing <- krak2_output_file_taxid %>% subset(taxid %in% taxa$taxid)


krak2_output_file_to_save <- krak2_output_file %>% subset(V2 %in% krak2_output_file_taxid_denosing$V2)

print(opt$out_krak2_denosing)

write.table(krak2_output_file_to_save,file = opt$out_krak2_denosing,col.names=F,row.names=F,quote=F,sep="\t")