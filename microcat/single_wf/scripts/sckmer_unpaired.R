library(optparse)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ShortRead)

option_list = list(
  make_option(c("--sample_name"), action="store", help = "sample name"),
  make_option(c("--fa1"), action="store", help = "path to fasta 1 file"),
  make_option(c("--microbiome_output_file"), action="store", help = "path to microbiome output file"),
  make_option(c("--cb_len"), action="store", default=16, help = "nucleutide length of cell barcodes"),
  make_option(c("--umi_len"), action="store", default=10, help = "nucleutide length of umis"),
  make_option(c("--kraken_report"), action="store", help = "path to standard kraken report"),
  make_option(c("--mpa_report"), action="store", help = "path to standard kraken mpa report"),
  make_option(c("--ranks"), action="store", default = c('G', 'S'), help = "taxa ranks to analyze"),
  make_option(c("--host"), action="store", default = 9606, help = "host taxid to exclude"),
  make_option(c("--min_frac"), action="store", default = 0.5, help = "minimum fraction of kmers directly assigned to taxid to use read"),
  make_option(c("--kmer_len"), action="store", default = 35, help = "Kraken kmer length"),
  make_option(c("--nsample"), action="store", default = 1000, help = "Max number of barcodes to sample per taxa"),
  make_option(c("--kmer_test_file"), action="store", help = "kmer test output filename"),
  make_option(c("--kmer_file"), action="store", help = "output kmer filename")
)

opt = parse_args(OptionParser(option_list = option_list))
# opt = list(fa1 = '/data/scRNA/EarlyGastricCancer-EGC/results/03.classifier/rmhost_extracted_classified_output/SRR9713114/SRR9713114_kraken2_extracted_classified.fq',
#            microbiome_output_file = '/data/scRNA/EarlyGastricCancer-EGC/results/03.classifier/rmhost_extracted_classified_output/SRR9713114/SRR9713114_kraken2_extracted_classified_output.txt',
#            cb_len = 16,
#            umi_len = 10,
#            kraken_report = '/data/scRNA/EarlyGastricCancer-EGC/results/03.classifier/rmhost_kraken2_report/custom/SRR9713114/SRR9713114_kraken2_report.txt',
#            mpa_report = '/data/scRNA/EarlyGastricCancer-EGC/results/03.classifier/rmhost_kraken2_report/mpa/SRR9713114/SRR9713114_kraken2_mpa_report.txt',
#            ranks = c('G', 'S'),
#            host = 9606, 
#             min_frac = 0.5,
#            kmer_len = 35,
#            nsample = 1000)

reads1 = readFasta(opt$fa1)
sequences1 = sread(reads1) 
# reads2 = readFasta(opt$fa2)
# sequences2 = sread(reads2) 

headers = ShortRead::id(reads1)
# [1]    39 SRR9713114.29152797 kraken:taxid|311790
# [2]    38 SRR9713114.89189781 kraken:taxid|54180
# [3]    40 SRR9713114.17845088 kraken:taxid|1675787
# [4]    40 SRR9713114.71942326 kraken:taxid|1547445
# [5]    39 SRR9713114.14176755 kraken:taxid|198719
barcode = substr(sequences1, 1, opt$cb_len)
umi = substr(sequences1, opt$cb_len+1, opt$cb_len + opt$umi_len)
taxid = gsub('.*taxid\\|', '', headers)
#311790
id = gsub('\\s.*', '', headers)
#SRR9713114.29152797

kr = read.delim(opt$kraken_report, header = F)
#      V1      V2      V3       V4      V5 V6     V7                   V8
# 1 59.19 3396557 3396557        0       0  U      0         unclassified
# 2 40.81 2341832  428499 13572948 4666916  R      1                 root
# 3 33.31 1911265   51602  6933013 3638862 R1 131567   cellular organisms
# 4 30.90 1773412   18111  6285824 3372655  D   2759            Eukaryota
# 5 29.95 1718595    3993  6048722 3213287 D1  33154         Opisthokonta
# 6 29.73 1706178      34  5997474 3173170  K  33208              Metazoa
kr = kr[-c(1:2), ] %>% mutate(V8 = trimws(V8)) %>% mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
#      V1      V2    V3      V4      V5 V6     V7                 V8
# 3 33.31 1911265 51602 6933013 3638862 R1 131567 cellular_organisms
# 4 30.90 1773412 18111 6285824 3372655  D   2759          Eukaryota
# 5 29.95 1718595  3993 6048722 3213287 D1  33154       Opisthokonta
# 6 29.73 1706178    34 5997474 3173170  K  33208            Metazoa
# 7 29.73 1706066   717 5996906 3173170 K1   6072          Eumetazoa
# 8 29.69 1703787 30772 5988194 3170988 K2  33213          Bilateria
mpa = read.delim(opt$mpa_report, header = F)
#                                                                         V1
# 1                                                             k__Eukaryota
# 2                                                  k__Eukaryota|k__Metazoa
# 3                                      k__Eukaryota|k__Metazoa|p__Chordata
# 4                          k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia
# 5              k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates
# 6 k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates|f__Hominidae
#        V2
# 1 1773412
# 2 1706178
# 3 1595687
# 4 1390257
# 5  899853
# 6  334715
mpa$taxid = NA
#get taxid

#rename mpa style df tax name
for(i in 2:nrow(mpa)){
  t_names = mpa[i,1] %>% as.character() %>% 
  #"k__Eukaryota|k__Metazoa|p__Chordata|c__Mammalia|o__Primates"
    strsplit('\\|') %>% 
    unlist() %>% 
    str_remove('.*__') %>% 
    str_replace_all('[^[:alnum:]]', '_') 
  #map taxid from krak report
  mpa$taxid[i] = paste0('*', paste(kr$V7[match(t_names, kr$V8)], collapse = '*'), '*')
}

microbiome_output_file = read.delim(opt$microbiome_output_file, header = F)
#   V1                   V2                                V3  V4
# 1  C SRR9713114.103190677      Pentapetalae (taxid 1437201) 150
# 2  C  SRR9713114.35509535     Raphanus sativus (taxid 3726) 150
# 3  C  SRR9713114.79016790    Vigna unguiculata (taxid 3917) 150
# 4  C SRR9713114.102607713 Hemileia vastatrix (taxid 203904) 150
# 5  C  SRR9713114.22922703 Nicotiana attenuata (taxid 49451) 150
# 6  C  SRR9713114.54362314    Solanum tuberosum (taxid 4113) 150
#                                                      V5
# 1               0:54 131567:1 0:52 8289:3 1437201:5 0:1
# 2 1:13 0:9 1:1 0:5 1:5 0:2 1:2 0:1 1:5 0:39 3726:1 0:33
# 3                               9526:1 0:65 3917:3 0:47
# 4                       0:34 37293:1 0:37 203904:5 0:39
# 5                          0:20 49451:4 0:8 9347:2 0:82
# 6                        0:35 207598:2 0:22 4113:5 0:52

#remove V1
microbiome_output_file = microbiome_output_file %>% 
  select(-V1) %>% 
  separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
  mutate(taxid = str_remove(taxid, '\\)') %>% trimws(),
         name = trimws(name))
#                     V2                name   taxid  V4
# 1 SRR9713114.103190677        Pentapetalae 1437201 150
# 2  SRR9713114.35509535    Raphanus sativus    3726 150
# 3  SRR9713114.79016790   Vigna unguiculata    3917 150
# 4 SRR9713114.102607713  Hemileia vastatrix  203904 150
# 5  SRR9713114.22922703 Nicotiana attenuata   49451 150
# 6  SRR9713114.54362314   Solanum tuberosum    4113 150
#                                                      V5
# 1               0:54 131567:1 0:52 8289:3 1437201:5 0:1
# 2 1:13 0:9 1:1 0:5 1:5 0:2 1:2 0:1 1:5 0:39 3726:1 0:33
# 3                               9526:1 0:65 3917:3 0:47
# 4                       0:34 37293:1 0:37 203904:5 0:39
# 5                          0:20 49451:4 0:8 9347:2 0:82
# 6                        0:35 207598:2 0:22 4113:5 0:52
tx = kr$V7[kr$V6 %in% opt$ranks] %>% setdiff(opt$host)
tx = microbiome_output_file$taxid[microbiome_output_file$taxid %in% tx] %>% unique()
#taxid in G and S

barcode_kmer = list()
counter = 0 
for(taxa in tx){
  counter = counter + 1
  # cat(paste('\r', 'taxa processed:', round(counter/length(tx)*100, 3), '%'))
  
  lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
  #"*2759*33090*35493*3398*3699*3700*3725*3726*"
    str_extract(paste0('\\*', taxa, '\\*.*')) %>%
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  
  full.lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  

  #                          V2             name taxid  V4
  # 2     SRR9713114.35509535 Raphanus sativus  3726 150
  # 406  SRR9713114.101942223 Raphanus sativus  3726 150
  # 445   SRR9713114.58434838 Raphanus sativus  3726 150
  # 544   SRR9713114.65313981 Raphanus sativus  3726 150
  # 1229   SRR9713114.1591813 Raphanus sativus  3726 150
  # 1708  SRR9713114.84769329 Raphanus sativus  3726 150
  #                                                         V5
  # 2    1:13 0:9 1:1 0:5 1:5 0:2 1:2 0:1 1:5 0:39 3726:1 0:33
  # 406                        0:29 1:3 9639:1 0:73 3726:5 0:5
  # 445                           0:8 3726:5 0:2 33213:1 0:100
  # 544                       0:1 1:3 3726:5 0:75 10117:2 0:30
  # 1229                               3726:3 0:92 9606:2 0:19
  # 1708                                  1:3 0:20 3726:1 0:92
  
  out = subset(microbiome_output_file, taxid %in% lin) %>% separate(V5, into = c('r1', 'r2'), sep = '\\|\\:\\|') 
  out$r1[str_which(out$r1, paste0(' ', opt$host, ':'))] = NA
  out$r2[str_which(out$r2, paste0(' ', opt$host, ':'))] = NA
  out$r1[out$r1 == ''] = NA; out$r2[out$r2 == ''] = NA
  out = subset(out, !is.na(r1) | !is.na(r2))
  out$r1 = trimws(out$r1)
  out$r2 = trimws(out$r2)
  
  
  if(nrow(out) == 0){next}
#                          V2             name taxid  V4
# 2     SRR9713114.35509535 Raphanus sativus  3726 150
# 406  SRR9713114.101942223 Raphanus sativus  3726 150
# 445   SRR9713114.58434838 Raphanus sativus  3726 150
# 544   SRR9713114.65313981 Raphanus sativus  3726 150
# 1708  SRR9713114.84769329 Raphanus sativus  3726 150
# 1774  SRR9713114.20498651 Raphanus sativus  3726 150
#                                                         r1   r2
# 2    1:13 0:9 1:1 0:5 1:5 0:2 1:2 0:1 1:5 0:39 3726:1 0:33 <NA>
# 406                        0:29 1:3 9639:1 0:73 3726:5 0:5 <NA>
# 445                           0:8 3726:5 0:2 33213:1 0:100 <NA>
# 544                       0:1 1:3 3726:5 0:75 10117:2 0:30 <NA>
# 1708                                  1:3 0:20 3726:1 0:92 <NA>
# 1774                                  1:4 0:27 3726:5 0:80 <NA>
  #get fasta idx
  i = which(id %in% out$V2)
  seq = data.frame(r1 = sequences1[i] %>% as.character())
#    r1
# 1 GCCAGCGCTGGGCCCCCGTCCTGACCTGAGCGGTTACCACCAGCCCCAGGCCTGCGGAGGCGCTAGTCCACCAGAGCCCCCCCCCCCCCCCCTCCCCGCCCCCCCCCCACCCCCCCCCCCCCCCCCTCACCCCCCCCCCCCCCCCCTCAC
# 2 NTTTTTTTAATAGTCAGTCCAAATATGAGATGCGTTGTTACAGGAAGTCCCTTGTCATCCTAAAAGGCACCCCAGGTATCTCTAGGGGGAGTGGCCCAGTCTTCTGTCAAGTCCTCACAGGGGAGGTGCTAGCATTGGCTTGGGGTAAAT
# 3 TTTTTGTTTTCGGAAATACGATTTTATTTTATTTAATATTATAATTTCAGTTAATCTAAGTTTTGAGTTCTGTACATGGGATGTTTCACTGTGCCATGGACTGTGATGGAAGTTCTGCGTCGTATGTGTGTCTCTACACTGGAGCTCACG
# 4 AAGACTGCTAACCAGGCCTCTGACACCTTCTCTGGTGTCGGGAAAAAAGTCGGCCTCCTGAAATGACAGCGGGGAGACCGGGGTCGGCCTCCTGAAATGATGGCAGTGGAACGGGGGGGGCCGCCGGTCCACGCGCCGTCTGGGCGAGCG
# 5 AAGCAGTGGTATCAACGCAGAGTACCTGGGGGAGCAGGAGCAGGAGCGGGAGCGGGATCCCCAGCAGGAGCAGGAGCCGGGGCGGCCCGCGGACCGGGACCCCCCGGACCGGGGGCCGGGGGGGCGCGCCCCGGGCCCGGGCCCGCGGCC
# 6 GCTTTCACTTTCTCTTTTGTTTTAAATGACTCATAGGTCCCTGACATTTTTTTGATTATTTTTTGCTACAGGGCTGGTAGACTCTGAGTGTATATAGATTAAGGCGTGGTGGTGTCAGGAGGCAGGCTGGGAAGGCTATTGTGGTGGGCT
  barcode.x = barcode[i]
  umi.x = umi[i]
  
  tax.df = list()
  if(nrow(out) > opt$nsample){n = sample(nrow(out), opt$nsample)} else {n=1:nrow(out)}
  counter2 = 0
  
  for(i in n){
    # cat(paste('\r', 'barcodes processed:', round(counter2/length(n)*100, 3), '%   '))
    for(mate in c('r1', 'r2')){
      #get pos
            # pos
# 1    1:13
# 2     0:9
# 3     1:1
# 4     0:5
# 5     1:5
# 6     0:2
# 7     1:2
# 8     0:1
# 9     1:5
# 10   0:39
# 11 3726:1
# 12   0:33
      r = data.frame(pos = out[[mate]][i] %>% strsplit('\\s') %>% unlist()) %>% 
        separate(pos, into = c('taxid', 'nkmer'), sep = ':', convert = T) 
      # if(any(r$taxid %in% c(0, full.lin) == F)){next}
#          taxid nkmer
# 1      1    13
# 2      0     9
# 3      1     1
# 4      0     5
# 5      1     5
# 6      0     2
# 7      1     2
# 8      0     1
# 9      1     5
# 10     0    39
# 11  3726     1
# 12     0    33
      r = r %>% 
        mutate(fkmer = nkmer/sum(nkmer),
               nt_start = cumsum(nkmer) - nkmer + 1,
               nt_end = cumsum(nkmer) + opt$kmer_len - 1) %>% 
        mutate(nt_len = nt_end - nt_start + 1) %>% 
        ungroup() %>% 
        # subset(taxid == taxa)
        subset(taxid %in% c(0,full.lin))
          # taxid nkmer      fkmer      nt_start nt_end nt_len
          # 2      0     9 0.07758621       14     56     43
          # 4      0     5 0.04310345       24     62     39
          # 6      0     2 0.01724138       34     69     36
          # 8      0     1 0.00862069       38     72     35
          # 10     0    39 0.33620690       44    116     73
          # 11  3726     1 0.00862069       83    117     35
          # 12     0    33 0.28448276       84    150     67
      if(sum(r$fkmer) < opt$min_frac){next}
      #i seq index
      counter2 = counter2 + 1

      if(nrow(r) > 0){
        kmer = c()
        for(k in 1:nrow(r)){
          for(m in 1:r$nkmer[k]){
            #r$nkmer[k]=9
            # seq[[mate]][i]
            # [1] "GCCAGCGCTGGGCCCCCGTCCTGACCTGAGCGGTTACCACCAGCCCCAGGCCTGCGGAGGCGCTAGTCCACCAGAGCCCCCCCCCCCCCCCCTCCCCGCCCCCCCCCCACCCCCCCCCCCCCCCCCTCACCCCCCCCCCCCCCCCCTCAC"
            kmer = c(kmer, substr(seq[[mate]][i], r$nt_start[k] + m - 1, r$nt_start[k] + m + opt$kmer_len - 2))
          }
        }
        tax.df[[counter2]]=data.frame(barcode = barcode.x[i], taxid = taxa, k = kmer, n = sum(r$nt_len[r$taxid %in% lin]))
      } else {tax.df[[counter2]] = data.frame(barcode = barcode.x[i], taxid = taxa, k = NA, n = NA)}
    }
  }
  
  if(length(tax.df) == 0){next}
  #actually it calcaute different uniq!!!!!
  tax.df = bind_rows(tax.df) %>%
    tibble() %>% 
    subset(!is.na(k)) %>%
    group_by(barcode, taxid) %>%
    summarize(kmer = length(k),
              uniq = length(unique(k)),
              .groups = 'keep')
  
  barcode_kmer[[counter]] = tax.df
  # cat('\n')
}

barcode_kmer = rbindlist(barcode_kmer)

# c = barcode_kmer %>%
#   subset(kmer > 1) %>%
#   group_by(taxid) %>%
#   mutate(nn = n()) %>%
#   subset(nn > 2) %>%
#   group_by(taxid) %>%
#   summarize(r = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$estimate,
#             p = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$p.value) %>%
#   mutate(p = p.adjust(p)) 

write.table(barcode_kmer, file = opt$kmer_file, quote = F, row.names = F,sep="\t")


########## kmer test
# kraken report
report = read.delim(opt$kraken_report, header = F)
report$V8 = trimws(report$V8)
report[report$V8 %in% c('Homo sapiens', 'Bacteria', 'Fungi', 'Viruses'), ]
# sckmer data
length(unique(report$V8[report$V6 %in% c('G', 'S')]))

# barcode k-mer correlation tests on taxonomy IDs detected on >3 barcodes and with >1 k-mer
# taxid = NCBI taxonomic ID, r = Spearman correlation coefficient, p = adjusted p-value
kmer_test = barcode_kmer %>% 
  subset(kmer > 1) %>%
  group_by(taxid) %>%
  mutate(nn = n()) %>%
  subset(nn > 2) %>% 
  group_by(taxid) %>%
  summarize(rho = cor.test(kmer, uniq, method = 'spearman')$estimate,
            p_value = cor.test(kmer, uniq, method = 'spearman')$p.value,
            .groups='keep') %>%
  mutate(p_adjust = p.adjust(p_value))

kmer_test$name = report$V8[match(kmer_test$taxid, report$V7)] # add taxa names 

write.table(kmer_test, file = opt$kmer_test_file, quote = F, row.names = F,sep="\t")

paste('Done')


