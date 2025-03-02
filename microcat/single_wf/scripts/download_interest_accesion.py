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
from multiprocessing import Pool, Manager
import urllib.request
from urllib.parse import urlparse
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



def url_for_accession(accession, *, verbose=False, quiet=False):
    accsplit = accession.strip().split("_")
    assert len(accsplit) == 2, f"ERROR: '{accession}' should have precisely one underscore!"

    db, acc = accsplit
    if '.' in acc:
        number, version = acc.split(".")
    else:
        number, version = acc, '1'
    number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"
    if verbose:
        print(f"opening directory: {url}", file=sys.stderr)

    with urllib.request.urlopen(url) as response:
        all_names = response.read()

    if verbose:
        print("done!", file=sys.stderr)

    all_names = all_names.decode("utf-8")

    full_name = None
    for line in all_names.splitlines():
        if line.startswith(f'<a href='):
            name=line.split('"')[1][:-1]
            db_, acc_, *_ = name.split("_")
            if db_ == db and acc_.startswith(acc):
                full_name = name
                break

    if full_name is None:
        return None
    else:
        url = "htt" + url[3:]
        return (
            f"{url}/{full_name}/{full_name}_genomic.fna.gz"
        )

def download_genomes(acc):
    download =False
    acc2seqid_local = {}

    # if complete:
    #         # if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome':
    #         if assembly_summaries.loc[acc,'assembly_level'] != 'Complete Genome' and assembly_summaries.loc[acc,'assembly_level'] != 'Chromosome':
    #             logger.error(f"Didn get {acc} because it wasn't complete or Chromosom") 
    #             return

    try:

        ftp_path = url_for_accession(acc)
        aname = ftp_path.split('/')[-1]
        file_path = os.path.join(download_folder,aname)
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

        logger.info(f'Save with {acc}', status='complete')
        fna_file = gzip.open(file_path,'rt')
        with open(interest_fna, 'a') as library_file:
            for record in SeqIO.parse(fna_file,"fasta"):
                acc2seqid_local[record.id] = acc
                SeqIO.write(record,library_file, "fasta")
    except Exception as e:
        logger.error(f"Error downloading {acc} : {e}")

    return acc2seqid_local  # Return empty or partial dict on error

def run_multiprocessing(func, acc_list, n_processors):
    acc2seqid = {}
    with Pool(processes=n_processors) as pool:
        results = pool.map(func, acc_list)  # Collect results from each process
    
    for result in results:
        acc2seqid.update(result)  # Merge each result dictionary into acc2seqid
    
    return acc2seqid

# def run_multiprocessing(func, i, n_processors):
#     with Pool(processes=n_processors) as pool:
#         return pool.map(func, i)
    # with Manager() as manager:
    #     acc2seqid = manager.dict()
    #     with Pool(processes=n_processors) as pool:
    #         pool.starmap(func, [(acc, acc2seqid) for acc in acc_list])
    #     return dict(acc2seqid)

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
parser.add_argument('--project', dest='project', 
                    help="project name")
parser.add_argument('--interest_fna', dest='interest_fna', 
                    help="interest_fna")
parser.add_argument('--folder', dest='folder',
                    help="name of the folder to download the genomes to. If this already exists, the genomes will be added to it. By default this is ncbi_genomes")
parser.add_argument('--acc2tax', dest='acc2tax',
                    help="Accession to TAXID databases")
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
interest_fna = args.interest_fna

logger.info('Downloading assembly summary', status='run')
# print('Starting processing')
domains = ['bacteria','archaea','fungi','viral']
# 获取当前运行时的工作目录路径
current_working_directory = os.getcwd()+"/"
# print(f"workdir：{current_working_directory}")


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

candidate_species = pd.read_csv(args.candidate,sep="\t",names = ['accesion'])
if not download_genomes:
    print("Finished running everything. The genomes haven't been downloaded because you didn't want to download them.")
    sys.exit()


# global acc2seqid
candidate_species.set_index('accesion', inplace=True)
selected_assembly_summaries_map = assembly_summaries[assembly_summaries.index.isin(candidate_species.index)]

origin_assembly_summaries = pd.read_csv("/data/comics-sucx/database/microcat/assembly_summaries.csv", sep=',',index_col=0)

def main():
    # Collect and update acc2seqid from multiprocessing results
    acc2seqid = run_multiprocessing(download_genomes, candidate_species.index.values, int(n_processors))
    # run_multiprocessing(download_genomes, existing_sequences,assembly_summaries.index.values,args.library_fna, int(n_processors))
    with open(args.acc2tax, 'w') as f:
        for seqid ,acc in acc2seqid.items():
            try:
                taxid = origin_assembly_summaries.loc[acc,'taxid']
                f.write(f"{seqid}\t{taxid}\n")
            except:
                logger.error(f"Couldn't find {acc} in the assembly summaries")
                f.write(f"{seqid}\t{acc}\n")
                continue


if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()

    # taxonomy = pd.DataFrame(taxonomy).transpose()
    # taxonomy.to_csv(os.path.join(folder, "taxonomy",'download_genomes.tsv'), sep='\t', index=False, header=False)
    
    # since acc2seqid[acc].append(record.id)
    

    logger.info('Finished downloaded', status='complete') 
    # print("Finished running everything. The genomes should be downloaded in "+download_folder+" and the list of these and their taxonomy is in genomes.tsv \nA log of any genomes that couldn't be downloaded is in logfile.txt")