"""
This script is used to download sequences from a certain domain that are present in NCBI refseq (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/).
It requires the additional packages pandas and Biopython.

The script downloads the assembly summary file and creates a log file that indicates whether each sequence was downloaded or not.

Command line arguments:
    --complete: Choose whether to only download complete genomes or all genomes. Default is True, meaning only complete genomes are downloaded.
    --candidate: Candidate species.
    --library_report: Library report.
    --seqid2taxid: Seqid to taxid.
    --library_fna: Library FNA file.
    --project: Project name.
    --interest_fna: Interest FNA file.
    --acc2tax: Accession to TAXID databases.
    --folder: Name of the folder to download the genomes to. If the folder already exists, the genomes will be added to it. By default, this is "ncbi_genomes".
    --log_file: File to write the log to. Default is "logfile_download_genomes.txt".
    --verbose: Enable detailed print.
    --processors: Number of processors to use for renaming genome files. Default is 1.

"""

#!/usr/bin/env python3
from collections import defaultdict
import os
import pandas as pd
import argparse
import sys
import logging
from Bio import SeqIO
import gzip
import json
import subprocess
import requests
import json
from multiprocessing import Pool
from multiprocessing import freeze_support
# Create a logger object
logger = logging.getLogger('my_logger')

# Create a formatter object with the desired log format
log_format = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

# Create a handler and add the formatter to it
console_handler = logging.StreamHandler()  # Output logs to the console
console_handler.setFormatter(log_format)

# Add the handler to the logger object
logger.addHandler(console_handler)

# Customize logger.info function to include status
def custom_log(level, msg, *args, status=None):
    if status:
        msg = f'({status}) {msg}'  # Concatenate the message and status
    logger.log(level, msg, *args)

# Bind the custom_log function to the logger object for different log levels
logger.info = lambda msg, *args, status=None: custom_log(logging.INFO, msg, *args, status=status)
logger.warning = lambda msg, *args, status=None: custom_log(logging.WARNING, msg, *args, status=status)
logger.error = lambda msg, *args, status=None: custom_log(logging.ERROR, msg, *args, status=status)
logger.debug = lambda msg, *args, status=None: custom_log(logging.DEBUG, msg, *args, status=status)


def select_genomes(group):
    # First get genome from subspecies
    filtered_genomes = group[group["acc"].isin(specifice_genome_list)].head(5)

    priority_categories = ['representative genome', 'reference genome']
    
    remaining_needed = 5 - len(filtered_genomes)
    if len(filtered_genomes) >= 3:
        represent_genomes = group[group['refseq_category'].isin(priority_categories)].head(2)
    else:
        represent_genomes = group[group['refseq_category'].isin(priority_categories)].head(5)
        
    filtered_genomes = pd.concat([filtered_genomes, represent_genomes])
    remaining_needed = 5 - len(filtered_genomes)

    if remaining_needed > 0:
        # 如果找到的基因组不足5个，先从包含 'FDAARGOS' 的基因组中选择
        contains_fdaargos = group[group['infraspecific_name'].str.contains('FDAARGOS', na=False)]
        fdaargos_filtered = contains_fdaargos[~contains_fdaargos.index.isin(filtered_genomes.index)]
        filtered_genomes = pd.concat([filtered_genomes, fdaargos_filtered.head(remaining_needed)])

        remaining_needed = 5 - len(filtered_genomes)

        if remaining_needed > 0:
            # 如果包含 'FDAARGOS' 的基因组仍不足5个，从剩余基因组中继续选择，按照指定的排序条件
            remaining_genomes = group[~group.index.isin(filtered_genomes.index)]
            sorted_remaining_genomes = remaining_genomes.sort_values(by=['assembly_level', '#assembly_accession']).head(remaining_needed)
            filtered_genomes = pd.concat([filtered_genomes, sorted_remaining_genomes])

    return filtered_genomes
    
def download_genomes(acc):
    download =False
    if complete:
            # if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome':
            if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome' and assembly_summaries.loc[acc,'assembly_level'] != 'Chromosome':
                logger.error(f"Didn get {acc} because it wasn't complete or Chromosom") 
                return
            
    if acc not in existing_sequences:
        try:
            ftp_path = assembly_summaries.loc[acc, 'ftp_path']
            aname = ftp_path.split('/')[-1]+'_genomic.fna.gz'
            file_path = os.path.join(download_folder,aname)
            ftp_path = ftp_path + "/" + aname
            attempts = 0
            
            if not os.path.exists(file_path) or os.stat(file_path).st_size == 0:
                download = True
            else:
                logger.info(f"Already had {acc} as file {aname} so didn't download it again",status='complete')
                attempts = 4
                try:
                    fna_file = gzip.open(file_path,'rt')
                except:
                    logger.error(f"Error in open {acc} :", e)
                    download = True
                    # remove the file
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        attempts = 0
            if download:
                while attempts < 3:
                    try:
                        # 构建下载命令
                        command = ["wget", "-q", ftp_path, "-O", os.path.join(download_folder,aname)]
                        # 执行下载命令，并检查返回码
                        subprocess.run(command, check=True, capture_output=True, text=True)
                        logger.info(f"Successfully downloaded {acc}", status='complete')
                        break
                    except subprocess.CalledProcessError as e:
                        # 处理调用命令时的异常
                        logger.error(f"Error executing download {acc} :", e)
                        attempts += 1
                        # 如果尝试了3次还是失败，抛出异常
                        if attempts == 3:
                            raise Exception(f"Failed to download {acc} after 3 attempts")
                    except Exception as e:
                        # 处理其他异常
                        logger.error("Unexpected error:", e)
                        attempts += 1
                        # 如果尝试了3次还是失败，抛出异常
                        if attempts == 3:
                            raise Exception(f"Failed to download {acc} after 3 attempts")

            taxonomy[0].append(acc), taxonomy[1].append(assembly_summaries.loc[acc, 'Taxonomy'])
            logger.info(f'Save with {acc}', status='complete')
            fna_file = gzip.open(file_path,'rt')
            with open(library_fna, 'a') as library_file:
                for record in SeqIO.parse(fna_file,"fasta"):
                    SeqIO.write(record,library_file, "fasta")

        except Exception as e:
            logger.error(f"Error downloading {acc} : {e}")


    else:
        logger.info(f"Already had {acc} as file {aname} so didn't download it again",status='complete')
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
parser.add_argument('--project', dest='project', 
                    help="project name")
parser.add_argument('--interest_fna', dest='interest_fna', 
                    help="interest_fna")
parser.add_argument('--acc2tax',
                    help="accession to TAXID databases")
parser.add_argument('--folder', dest='folder',
                    help="name of the folder to download the genomes to. If this already exists, the genomes will be added to it. By default this is ncbi_genomes")
parser.add_argument('--log_file', dest='log_file', default='logfile_download_genomes.txt',
                    help="File to write the log to")
parser.add_argument('--verbose', action='store_true', help='Detailed print')
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to rename genome files")

args = parser.parse_args()
# Set log level based on command line arguments
if args.verbose:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)

# Create a file handler and add the formatter to it
file_handler = logging.FileHandler(args.log_file)  # Output logs to the specified file
file_handler.setFormatter(log_format)
logger.addHandler(file_handler)

complete = args.complete
folder = args.folder
log_file = args.log_file
n_processors =args.proc
library_fna = args.library_fna
# print('Starting processing')
domains = ['bacteria','archaea','fungi','viral']
# 获取当前运行时的工作目录路径
current_working_directory = os.getcwd()+"/"
# print(f"workdir：{current_working_directory}")

logger.info('Downloading assembly summary', status='run')

taxonomy_folder = os.path.join(folder, "taxonomy")

# 检查 taxonomy 文件夹是否存在
if not os.path.exists(taxonomy_folder):
    # 如果不存在，创建它
    os.makedirs(taxonomy_folder)

assembly_summaries = []
# os.chdir(taxonomy_folder)

for domain in domains:
    try:
        summary = pd.read_csv(os.path.join(taxonomy_folder,str(domain)+"_assembly_summary.txt"), sep='\t', header=1, index_col=0)
        summary = summary.loc[:, ['taxid', 'species_taxid', 'organism_name', 'assembly_level','ftp_path','infraspecific_name','refseq_category']]
        summary['Domain'] = str(domain)
        assembly_summaries.append(summary)
    except (ValueError, KeyError) as e:
        logger.error(f"Unable to read {domain}_assembly_summary.txt with {e}")
        sys.exit()

logger.info('Finished downloading interest assembly summary', status='complete')        
assembly_summaries = pd.concat(assembly_summaries)
logger.info('Finished Joining the assembly summaries', status='complete')   

# taxid_change_log = pd.read_csv(os.path.join(taxonomy_folder,"merged.dmp"),sep="|",header=None,names=["origin","change","empty"])
# taxid_change_log = pd.read_csv("/data/scRNA_analysis/benchmark/Galeano2022_16s/database/bowtie2/taxonomy/merged.dmp",sep="|",header=None,names=["origin","change","empty"])

# summary = pd.read_csv("/data/project/host-microbiome/microcat_bowtie2/database/taxonomy/bacteria_assembly_summary.txt", sep='\t', header=1, index_col=0)
# summary = summary.loc[:, ['taxid', 'species_taxid', 'organism_name', 'assembly_level', 'ftp_path','refseq_category']]

logger.info('Downloading tax dump', status='run')  


# if not os.path.exists(os.path.join(taxonomy_folder,'rankedlineage.dmp')) :

#     try:
#         # 下载新的 taxdump 文件
#         download_command = ["wget", "https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz", "-O", os.path.join(taxonomy_folder, "new_taxdump.tar.gz")]
#         subprocess.run(download_command, check=True)

#         # 解压缩 taxdump 文件
#         extract_command = ["tar", "-xf", os.path.join(taxonomy_folder, "new_taxdump.tar.gz"),"-C", taxonomy_folder]
#         subprocess.run(extract_command, check=True)

#         print("Download and extraction completed successfully.")

#     except subprocess.CalledProcessError as e:
#         # 处理调用命令时的异常
#         logger.error(f"Error executing download tax dump: {e}")
#     except Exception as e:
#         # 处理其他异常
#         logger.error(f"Unexpected error executing download tax dump: {e}")
# else:
#     logger.info("Already got full lineage", status='complete')

# try:
#     full_lineage = pd.read_csv(os.path.join(taxonomy_folder,'rankedlineage.dmp'), sep='|', header=None, index_col=0)
#     logger.info('Got the full lineage from the current NCBI taxdump', status='complete')
#     if os.path.exists(os.path.join(taxonomy_folder, "new_taxdump.tar.gz")): 
#         os.remove(os.path.join(taxonomy_folder, "new_taxdump.tar.gz"))
#         logger.info('Removed the other files from the taxdump folder', status='complete') 
    
# except (ValueError, KeyError) as e:
#     logger.error(f"Couldn't get the full lineage from the current NCBI taxdump {e}")
#     sys.exit()
assembly_summaries['Taxonomy'] = ''
# os.chdir(current_working_directory)

download_folder = os.path.join(folder, "download")

# 检查 taxonomy 文件夹是否存在
if not os.path.exists(download_folder):
    # 如果不存在，创建它
    os.makedirs(download_folder)
# os.chdir(download_folder)

library_folder = os.path.join(folder, "library")

# 检查 library 文件夹是否存在
if not os.path.exists(library_folder):
    # 如果不存在，创建它
    os.makedirs(library_folder)

candidate_species = pd.read_csv(args.candidate,sep="\t")
candidate_species = candidate_species.sort_values(by="reads",ascending=False)
candidate_species = candidate_species[~candidate_species.duplicated(subset=['ncbi_taxa'],  keep='first')]


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

# genome_list = defaultdict(list)
# genome2taxid = dict()
# # some taxid may only have subspecies taxid
# failed_taxids = []
# for taxid in desired_taxid_list:
#     taxid = str(taxid)
#     genome_list[taxid] = taxid2genomeid[taxid]
#     if len(genome_list[taxid]) > 0:
#         for genome in genome_list[taxid]:
#             genome2taxid[genome] = taxid
#     if len(genome_list[taxid]) == 0 :
#         failed_taxids.append(taxid)
# 构建defaultdict
desired_taxid_dict = defaultdict(set)
accession2taxid = defaultdict(dict)  # Use defaultdict to automatically create a set for new keys
genome_list = defaultdict(list)
genome2taxid = dict()
specifice_genome_list = set()
grouped = candidate_species.groupby('main_level_taxid') 
count = 0
for name, group in grouped:
    count += 1
    # if count > 3:
    #     break    
    if len(group) >= 1:
        for index, row in group.head(5).iterrows():
            main_level_taxid = str(row['main_level_taxid'])
            rank = row['rank']
            taxid = str(row['ncbi_taxa'])
            # if rank != "S":
            genome_list = taxid2genomeid[taxid]
            if len(genome_list) > 5:
                acc_list = library_report[library_report['seqid'].isin(genome_list)]['acc'].tolist()
                for acc in acc_list:
                    accession2taxid[acc]['main_level_taxid'] = main_level_taxid
                    accession2taxid[acc]['rank'] = rank
            elif len(genome_list) > 0:
                acc_list = library_report[library_report['seqid'].isin(genome_list)]['acc'].tolist()
                for acc in acc_list:
                    accession2taxid[acc]['main_level_taxid'] = main_level_taxid
                    accession2taxid[acc]['rank'] = rank
                    specifice_genome_list.add(acc)


# filtered_library_report = library_report[library_report['seqid'].isin(genome2taxid.keys())]
# assembly_summaries = assembly_summaries[assembly_summaries.index.isin(filtered_library_report['acc'])]
# selected_assembly_summaries  = assembly_summaries.copy()


# Filter the 'summary' DataFrame to only include rows where 'acc' matches 'library_report['acc']'
selected_assembly_summaries = assembly_summaries.copy()[assembly_summaries.index.isin(accession2taxid.keys())]
selected_assembly_summaries['acc'] = selected_assembly_summaries.index.tolist()
selected_assembly_summaries = selected_assembly_summaries.groupby('species_taxid').apply(select_genomes).reset_index(drop=True)
selected_assembly_summaries = selected_assembly_summaries.drop_duplicates(subset='acc')

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
# print(selected_assembly_summaries.columns.tolist())
# 将 'acc' 列的数据类型都转换为字符串
selected_assembly_summaries['acc'] = selected_assembly_summaries['acc'].astype(str)
# filtered_library_report['acc'] = filtered_library_report['acc'].astype(str)

# 然后再进行合并
selected_assembly_summaries_map = pd.merge(selected_assembly_summaries, library_report[['acc', 'seqid']], on='acc', how='left')

# selected_assembly_summaries_map = pd.merge(selected_assembly_summaries, filtered_library_report[['acc', 'seqid']], on='acc', how='inner')
selected_assembly_summaries_map = selected_assembly_summaries_map[["seqid","taxid"]]
selected_assembly_summaries_map['accession'] = selected_assembly_summaries_map['seqid'].str.replace(r'\.\d+', '', regex=True)
selected_assembly_summaries_map['accession.version'] = selected_assembly_summaries_map["seqid"]
selected_assembly_summaries_map['taxid'] = selected_assembly_summaries_map["taxid"]
selected_assembly_summaries_map = selected_assembly_summaries_map[['accession',"accession.version","taxid"]]
genome_list = set(library_report[library_report['acc'].isin(selected_assembly_summaries['acc'])]['seqid'].unique().tolist())
selected_assembly_summaries_map.to_csv(args.acc2tax,sep="\t",index=False)
selected_assembly_summaries.set_index('acc', inplace=True)


def main():
    run_multiprocessing(download_genomes, selected_assembly_summaries.index.values, int(n_processors))
    # run_multiprocessing(download_genomes, existing_sequences,assembly_summaries.index.values,args.library_fna, int(n_processors))


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()

    # taxonomy = pd.DataFrame(taxonomy).transpose()
    # taxonomy.to_csv(os.path.join(folder, "taxonomy",'download_genomes.tsv'), sep='\t', index=False, header=False)
    with open(args.library_fna, 'r') as library_file:
        with open(args.interest_fna, 'w') as interest_file:
            for record in SeqIO.parse(library_file, "fasta"):
                if record.id in genome_list:
                    SeqIO.write(record,interest_file, "fasta")
    logger.info('Finished downloaded', status='complete') 
    # print("Finished running everything. The genomes should be downloaded in "+download_folder+" and the list of these and their taxonomy is in genomes.tsv \nA log of any genomes that couldn't be downloaded is in logfile.txt")