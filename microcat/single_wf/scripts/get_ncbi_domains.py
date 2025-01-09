#!/usr/bin/env python3
from collections import defaultdict
import os
import pandas as pd
import argparse
import sys
from Bio import SeqIO
import gzip
from multiprocessing import Pool
from multiprocessing import freeze_support
# prefix = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 's1__']

# for acc in assembly_summaries.index.values:
#     taxid = assembly_summaries.loc[acc, 'taxid']
#     if taxid in full_lineage.index.values:
#         taxonomy = list(full_lineage.loc[taxid, :].values)
#         tax_string = ''
#         taxonomy.reverse()
#         count = 0
#         new_tax = []
#         for tax in taxonomy:
#             if isinstance(tax, str):
#                 new_tax.append(tax.replace('\t', ''))
#         taxonomy = new_tax
#         for t in range(len(taxonomy)):
#             if t != 0: tax_string += ';'
#             if taxonomy[t] == '':
#                 if t == 1 and assembly_summaries.loc[acc, 'Domain'] == 'protozoa': 
#                     taxonomy[t] = 'Protista'
#                 elif t < 7:
#                     taxonomy[t] = previous
#                 else:
#                     taxonomy[t] = assembly_summaries.loc[acc, 'organism_name']
#             tax_string += prefix[t]+taxonomy[t]
#             previous = taxonomy[t]
#         assembly_summaries.loc[acc, 'Taxonomy'] = tax_string

# assembly_summaries = assembly_summaries.rename(columns={'taxid':'NCBI_taxid', 'species_taxid':'NCBI_species_taxid'})
# assembly_summaries.to_csv('summary_to_download.csv')
# print('Got phylogeny and file paths for all domains to include. This is saved as summary_to_download.csv\n')



def write_log(message):
    with open(log_file, 'w+') as f:
        f.write(message+'\n')



def download_genomes(acc):
    if complete:
            # if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome':
            if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome' and assembly_summaries.loc[acc,'assembly_level'] != 'Chromosome':
                write_log("Didn't get "+acc+" because it wasn't complete or Chromosome")
                return
    try:
            ftp_path = assembly_summaries.loc[acc, 'ftp_path']
            aname = ftp_path.split('/')[-1]+'_genomic.fna.gz'
            if not os.path.exists(aname):
                ftp_path = ftp_path+'/'+aname
                cmd = 'wget -q '+ftp_path
                os.system(cmd)
                if os.path.exists(aname):
                    if acc not in existing_sequences:
                        fna_file = gzip.open(aname,'rt')
                        with open(library_fna, 'a') as library_file:
                            for record in SeqIO.parse(fna_file,"fasta"):
                                SeqIO.write(record,library_file, "fasta")
                    taxonomy[0].append(acc), taxonomy[1].append(assembly_summaries.loc[acc, 'Taxonomy'])
                    write_log("Write "+acc)
                else:
                    write_log("Couldn't get "+acc)
            else:
                if acc not in existing_sequences:
                    fna_file = gzip.open(aname,'rt')
                    with open(library_fna, 'a') as library_file:
                        for record in SeqIO.parse(fna_file,"fasta"):
                            SeqIO.write(record,library_file, "fasta")
                taxonomy[0].append(acc), taxonomy[1].append(assembly_summaries.loc[acc, 'Taxonomy'])
                write_log("Already had "+acc+" as file "+aname+" so didn't download it again")
    except:
            write_log("Couldn't get "+acc)
    return



def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

# def run_multiprocessing(func, existing_sequences, index_values,library_fna, n_processors):
#     args = zip(existing_sequences, index_values,library_fna)
#     with Pool(processes=n_processors) as pool:
#         return pool.starmap(func, args)
parser = argparse.ArgumentParser(description='This script is to download all sequences from a certain domain that are present in NCBI refseq (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/).\n\
                                This requires the additional packages pandas and Biopython\nAs long as the script is able to download the assembly summary file, then it will create a log_file that tells you about whether each sequence was downloaded or not\n\
                                Re-running it with additional domains will by default just add these to what you already have')
parser.add_argument('--complete', dest='complete', default=True, 
                    help="choose whether to only download complete genomes, or all genomes. Default is False, meaning all genomes are downloaded")
parser.add_argument('--candidate', dest='candidate', 
                    help="candidate species")
parser.add_argument('--library_report', dest='library_report', 
                    help="library_report")
parser.add_argument('--seqid2taxid', dest='seqid2taxid', 
                    help="seqid2taxid")
parser.add_argument('--library_fna', dest='library_fna', 
                    help="library_fna")
parser.add_argument('--interest_fna', dest='interest_fna', 
                    help="interest_fna")
parser.add_argument('--folder', dest='folder',
                    help="name of the folder to download the genomes to. If this already exists, the genomes will be added to it. By default this is ncbi_genomes")
parser.add_argument('--log_file', dest='log_file', default='logfile_download_genomes.txt',
                    help="File to write the log to")
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to rename genome files")

args = parser.parse_args()

complete = args.complete
folder = args.folder
log_file = args.log_file
n_processors =args.proc
library_fna = args.library_fna
print('Starting processing')
domains = ['bacteria','archaea','fungi','viral']

taxonomy_folder = os.path.join(folder, "taxonomy")

# 检查 taxonomy 文件夹是否存在
if not os.path.exists(taxonomy_folder):
    # 如果不存在，创建它
    os.makedirs(taxonomy_folder)

assembly_summaries = []
os.chdir(taxonomy_folder)
for domain in domains:
    try:
        if not os.path.exists(str(domain)+"_assembly_summary.txt"):
            os.system("wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+str(domain)+"/assembly_summary.txt")
            # os.system("wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+str(domain)+"/assembly_summary.txt")
            os.system("mv assembly_summary.txt "+str(domain)+"_assembly_summary.txt")
        else:
            print("Already got "+str(domain)+"_assembly_summary.txt")
        summary = pd.read_csv(str(domain)+"_assembly_summary.txt", sep='\t', header=1, index_col=0)
        summary = summary.loc[:, ['taxid', 'species_taxid', 'organism_name', 'assembly_level', 'ftp_path','refseq_category']]
        summary['Domain'] = str(domain)
        assembly_summaries.append(summary)
    except:
        print("Unable to download assembly_summary.txt for "+domain)
        
assembly_summaries = pd.concat(assembly_summaries)
print('Joined the assembly summaries\n')

# summary = pd.read_csv("/data/project/host-microbiome/microcat_bowtie2/database/taxonomy/bacteria_assembly_summary.txt", sep='\t', header=1, index_col=0)
# summary = summary.loc[:, ['taxid', 'species_taxid', 'organism_name', 'assembly_level', 'ftp_path','refseq_category']]
try:
    if not os.path.exists('rankedlineage.dmp'):
        # subprocess.call('wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz', shell=True)
        os.system('wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz')
        os.system('tar -xf new_taxdump.tar.gz')
    full_lineage = pd.read_csv('rankedlineage.dmp', sep='|', header=None, index_col=0)
    print('Got the full lineage from the current NCBI taxdump\n')
    other_files = ['new_taxdump.tar.gz']
    for f in other_files:
        if os.path.exists(f):
            os.system('rm '+f)
    print('Removed the other files from the taxdump folder')
except:
    print("Couldn't get the full lineage from the current NCBI taxdump\n")
    sys.exit()
assembly_summaries['Taxonomy'] = ''
download_folder = os.path.join(folder, "download")

# 检查 taxonomy 文件夹是否存在
if not os.path.exists(download_folder):
    # 如果不存在，创建它
    os.makedirs(download_folder)
os.chdir(download_folder)

library_folder = os.path.join(folder, "library")

# 检查 taxonomy 文件夹是否存在
if not os.path.exists(library_folder):
    # 如果不存在，创建它
    os.makedirs(library_folder)

candidate_species = pd.read_csv(args.candidate,sep="\t")
# candidate_species = pd.read_csv("/data/scRNA_analysis/benchmark/Galeano2022_HT29/results/03.classifier/rmhost_classified_qc/4_HT29Cells_BMixB_Fn_Ec_16S_S4/4_HT29Cells_BMixB_Fn_Ec_16S_S4_krak_study_denosing.txt",sep="\t")
desired_taxid_list = set(candidate_species["species_level_taxid"].unique())

library_report = pd.read_csv(args.library_report,sep="\t")
# library_report = pd.read_csv("/data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/library_report.tsv",sep="\t")
# 处理第二列
library_report['sci_name'] = library_report['Sequence Name'].str.split(',').str[0].str.split('.1 ').str[1]
library_report['species_name'] = library_report['sci_name'].str.split(' ').str[0] + " "+library_report['sci_name'].str.split(' ').str[1]
library_report['seqid'] = library_report['Sequence Name'].str.split(' ').str[0].str.replace('>', '')
library_report = library_report[library_report['#Library'].isin(["bacteria","viral","fungi","archaea"])]
library_report['acc'] = library_report['URL'].str.split('/').str[-1].str.extract(r'(GCF_\d+\.?\d*)')[0]
seqid_interest =  set(library_report['seqid'].unique())


taxid2genomeid = defaultdict(list)
with open(os.path.join(args.seqid2taxid), "r") as mapping:
    for line in mapping:
        info, taxid = line.strip().split("\t")
        info_tuple = info.split("|")
        if len(info_tuple) > 2:
            genome = info.split("|")[2]
        else:
            genome = info
        # only keep "bacteria","viral","fungi","archaea"
        if genome in seqid_interest:
            taxid2genomeid[taxid].append(genome)
# with open(os.path.join("/data/database/kraken2uniq_database/k2_pluspf_16gb_20231009/seqid2taxid.map"), "r") as mapping:
#     for line in mapping:
#         info, taxid = line.strip().split("\t")
#         info_tuple = info.split("|")
#         if len(info_tuple) > 2:
#             genome = info.split("|")[2]
#         else:
#             genome = info
#         # only keep "bacteria","viral","fungi","archaea"
#         if genome in seqid_interest:
#             taxid2genomeid[taxid].append(genome)

genome_list = defaultdict(list)
genome2taxid = dict()

for taxid in desired_taxid_list:
    taxid = str(taxid)
    genome_list[taxid] = taxid2genomeid[taxid]
    if len(genome_list[taxid]) > 0:
        for genome in genome_list[taxid]:
            genome2taxid[genome] = taxid

filtered_library_report = library_report[library_report['seqid'].isin(genome2taxid.keys())]

# Filter the 'summary' DataFrame to only include rows where 'acc' matches 'library_report['acc']'
assembly_summaries = assembly_summaries[assembly_summaries.index.isin(filtered_library_report['acc'])]

ref_assembly_summaries = assembly_summaries[assembly_summaries['refseq_category'].isin(['representative genome','reference genome'])]

# 根据taxid去重，只保留第一次出现的行
ref_assembly_summaries = ref_assembly_summaries.drop_duplicates(subset='taxid')
# 使用 ~ 操作符来反转条件，获取不包含在 ref_assembly_summaries 中的行
noref_assembly_summaries= assembly_summaries[~assembly_summaries['taxid'].isin(ref_assembly_summaries["taxid"])]

# 首先筛选出Complete Genome的行
complete_noref_assembly_summaries = noref_assembly_summaries[noref_assembly_summaries['assembly_level'] == 'Complete Genome']

# 从 noref_assembly_summaries 中移除 Complete Genome 的行
complete_noref_assembly_summaries = complete_noref_assembly_summaries.drop_duplicates(subset='taxid')

nocomplete_noref_assembly_summaries = noref_assembly_summaries[~noref_assembly_summaries['taxid'].isin(complete_noref_assembly_summaries["taxid"])]

nocomplete_noref_assembly_summaries  = nocomplete_noref_assembly_summaries.drop_duplicates(subset='taxid')

# 将 Complete Genome 的行重新添加到 noref_assembly_summaries
assembly_summaries = pd.concat([ref_assembly_summaries, complete_noref_assembly_summaries])
assembly_summaries = pd.concat([assembly_summaries, nocomplete_noref_assembly_summaries])

if not download_genomes:
    print("Finished running everything. The genomes haven't been downloaded because you didn't want to download them.")
    sys.exit()

existing_sequences = set()
# 如果library.fna文件存在，将其中的序列添加到集合中
if os.path.exists(args.library_fna):
    with open(args.library_fna, 'r') as library_file:
        for record in SeqIO.parse(library_file, "fasta"):
            existing_sequences.add(record.id)
else:
        library_file = open(args.library_fna, 'w')
        library_file.close()

taxonomy = [[], []]


def main():
    run_multiprocessing(download_genomes, assembly_summaries.index.values, int(n_processors))
    # run_multiprocessing(download_genomes, existing_sequences,assembly_summaries.index.values,args.library_fna, int(n_processors))


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()

    taxonomy = pd.DataFrame(taxonomy).transpose()
    taxonomy.to_csv(os.path.join(folder, "taxonomy",'download_genomes.tsv'), sep='\t', index=False, header=False)
    with open(args.library_fna, 'r') as library_file:
        with open(args.interest_fna, 'w') as interest_file:
            for record in SeqIO.parse(library_file, "fasta"):
                if record.id in genome2taxid:
                    SeqIO.write(record,interest_file, "fasta")
    print("Finished running everything. The genomes should be downloaded in "+download_folder+" and the list of these and their taxonomy is in genomes.tsv \nA log of any genomes that couldn't be downloaded is in logfile.txt")