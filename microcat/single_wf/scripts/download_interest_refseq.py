"""
This script is used to download sequences from a certain domain that are present in NCBI refseq.
It downloads sequences based on a kraken2 library report and candidate species list.
It requires the additional packages pandas and Biopython.

Command line arguments:
    --library_report: Library report file path.
    --candidate_species: Candidate species file path.
    --seqid2taxid: Seqid to taxid map file path.
    --library_fna: Library FNA file path.
    --project: Project name.
    --interest_fna: Interest FNA file path.
    --acc2tax: Accession to TAXID database output path.
    --folder: Name of the folder to download the genomes to. Default is "ncbi_genomes".
    --log_file: File to write the log to. Default is "logfile_download_genomes.txt".
    --verbose: Enable detailed print.
    --processors: Number of processors to use for downloading. Default is 1.
    --max_urls_per_taxid: Maximum number of URLs per taxid. Default is 10.
    --download_gff: Download GFF files along with FNA files. Default is False.
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
import subprocess
import requests
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

def download_with_requests(url, local_file):
    """Download a file using requests with streaming to save memory.
    
    Parameters
    ----------
    url : str
        URL to download from
    local_file : str
        Path to save the file to
        
    Returns
    -------
    bool
        True if download was successful, False otherwise
    """
    # TODO: add md5 check
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for 4XX/5XX responses
        
        with open(local_file, 'wb') as handle:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:  # Filter out keep-alive chunks
                    handle.write(chunk)
        return True
    except Exception as e:
        logger.error(f"Error downloading {url}: {e}")
        return False
    
def download_genome(row):
    """Only responsible for downloading genomes, not handling writing"""
    download = False
    seqid = row['seqid']
    taxid = row['taxid']
    accession = row['accession']
    fna_ftp_path = row['URL']
    # here we use the http path to download the fna file
    fna_http_path = fna_ftp_path.replace('ftp://', 'http://')
    try:
        fna_name = fna_http_path.split('/')[-1]
        fna_file_path = os.path.join(download_folder, fna_name)
        
        # Download fna file
        if not os.path.exists(fna_file_path) or os.stat(fna_file_path).st_size == 0:
            download = True
        else:
            logger.info(f"Already had {accession} as file {fna_name}", status='complete')
        
        if download:
            logger.info(f"Downloading {accession} fna file", status='run')
            attempts = 0
            while attempts < 3:
                try:
                    # Replace wget with requests
                    if download_with_requests(fna_http_path, fna_file_path):
                        logger.info(f"Successfully downloaded {accession} fna file", status='complete')
                        break
                    else:
                        raise Exception("Download failed")
                except Exception as e:
                    logger.error(f"Error downloading {accession}: {e}")
                    attempts += 1
                    # if have tried 3 times, raise an exception
                    if attempts == 3:
                        raise Exception(f"Failed to download {accession} after 3 attempts")

        # Download gff file only if requested
        if args.download_gff:
            gff_http_path = fna_http_path.replace('genomic.fna.gz', 'genomic.gff.gz')
            # gff_http_path = gff_ftp_path.replace('ftp://', 'http://')
            gff_name = gff_http_path.split('/')[-1]
            gff_file_path = os.path.join(download_folder, gff_name)

            if not os.path.exists(gff_file_path) or os.stat(gff_file_path).st_size == 0:
                attempts = 0
                while attempts < 3:
                    try:
                        # Download gff file
                        if download_with_requests(gff_http_path, gff_file_path):
                            logger.info(f"Successfully downloaded {accession} gff file", status='complete')
                            break
                        else:
                            raise Exception("Download failed")
                    except Exception as e:
                        logger.error(f"Error downloading {accession} gff: {e}")
                        attempts += 1
                        # if have tried 3 times, raise an exception
                        if attempts == 3:
                            raise Exception(f"Failed to download {accession} after 3 attempts")
            else:
                logger.info(f"Already had {accession} as file {gff_name}", status='complete')
            
        return {'success': True, 'path': fna_file_path, 'seqid': seqid, 'taxid': taxid, 'accession': accession}
    
    except Exception as e:
        logger.error(f"Error processing {accession}: {str(e)}")
        return {'success': False, 'seqid': seqid, 'taxid': taxid, 'accession': accession}

def run_multiprocessing(func, i, n_processors):
    with Pool(processes=n_processors) as pool:
        return pool.map(func, i)

parser = argparse.ArgumentParser(description='This script is to download all sequences from a certain domain that are present in NCBI refseq (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/).\n\
                                This requires the additional packages pandas and Biopython\nAs long as the script is able to download the assembly summary file, then it will create a log_file that tells you about whether each sequence was downloaded or not\n\
                                Re-running it with additional domains will by default just add these to what you already have')
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
parser.add_argument('--max_urls_per_taxid', dest='max_urls_per_taxid', default=5,
                    help="Maximum number of URLs per taxid")
parser.add_argument('--download_gff', action='store_true', default=False,
                    help="Download GFF files along with FNA files. Default is False.")

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

# Set up basic variables
folder = args.folder
log_file = args.log_file
n_processors =args.proc
library_fna = args.library_fna
max_urls_per_taxid = args.max_urls_per_taxid

# Create necessary directories
download_folder = os.path.join(folder, "download")
if not os.path.exists(download_folder):
    logger.info(f'Creating download folder {download_folder}', status='run')
    os.makedirs(download_folder)

# Load candidate species
logger.info('Loading candidate species', status='run')
candidate_species = pd.read_csv(args.candidate,sep="\t",dtype={'ncbi_taxa': str, 'main_level_taxid': str})
candidate_species = candidate_species.sort_values(by="reads",ascending=False)
# count the prevelence of each taxid
candidate_species['sample_count'] = candidate_species.groupby('ncbi_taxa')['ncbi_taxa'].transform('size')
# then remove the duplicates
candidate_species = candidate_species[~candidate_species.duplicated(subset=['ncbi_taxa'],  keep='first')]
logger.info(f'Loaded {len(candidate_species)} candidate species', status='complete')


# Load library report
logger.info('Loading library report', status='run')
library_report = pd.read_csv(args.library_report, sep="\t")

# Extract scientific name and seqid from sequence name
library_report['sci_name'] = library_report['Sequence Name'].str.split(',').str[0].str.split('.1 ').str[1]
library_report['seqid'] = library_report['Sequence Name'].str.split(' ').str[0].str.replace('>', '')
# Filter out non-bacterial, viral, fungal, or archaeal sequences
library_report = library_report[library_report['#Library'].isin(["bacteria","viral","fungi","archaea"])]
# Extract accession from URL
library_report['accession'] = library_report['URL'].str.split('/').str[-1].str.extract(r'(GCF_\d+\.?\d*)')[0]
# Get unique seqids of interest
seqid_interest =  set(library_report['seqid'].unique())

# Load seqid to taxid mapping if provided
seqid2taxid = {}
if args.seqid2taxid:
    logger.info('Loading kraken2 seqid to taxid mapping', status='run')
    with open(args.seqid2taxid, "r") as mapping:
        for line in mapping:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                info, taxid = parts
                info_parts = info.split("|")
                if len(info_parts) > 2:
                    seqid = info.split("|")[2]
                else:
                    seqid = info
                seqid2taxid[seqid] = taxid
    logger.info(f'Loaded {len(seqid2taxid)} seqid to taxid mappings', status='complete')
else:
    logger.error('No seqid to taxid mapping provided', status='error')
    sys.exit(1)


# Add taxids to library report if possible
if seqid2taxid:
    logger.info('Adding taxids to library report', status='run')
    library_report['taxid'] = library_report['seqid'].map(seqid2taxid)
    # check if there are any missing taxids
    missing_taxids = library_report[library_report['taxid'].isna()]
    if len(missing_taxids) > 0:
        logger.warning(f'Missing taxids for {len(missing_taxids)} sequences', status='check')
    logger.info(f'Added {len(library_report)} taxids to library report', status='complete')


# Get existing sequences to skip download
existing_sequences = set()
if os.path.exists(args.library_fna):
    with open(args.library_fna, 'r') as library_file:
        for record in SeqIO.parse(library_file, "fasta"):
            existing_sequences.add(record.id)
else:
    with open(args.library_fna, 'w') as library_file:
        logger.warning("No existing sequences found, will create empty library file and download all sequences", status='check')
        pass  # Create empty file


# Extract taxids from candidate species to use for filtering
taxids_to_download = set(candidate_species['ncbi_taxa'].astype(str).tolist())
main_level_taxids = set(candidate_species['main_level_taxid'].astype(str).tolist())

# Match library_report entries to candidate species
# First try direct taxid matching
matched_entries = []

# Group by main_level_taxid
grouped_candidates = candidate_species.groupby('main_level_taxid')

# Process each main_level_taxid group
## Init the list of taxids to download
taxids_to_download = defaultdict(set)
all_selected_entries = set()
# select the microbiome library report for genome selection
microbiome_library_report = library_report[library_report['#Library'].isin(["bacteria","viral","fungi","archaea"])]
for main_level_taxid, group in grouped_candidates:
    taxids_to_download[main_level_taxid] = set()

    # Here we sort the group by rank, sample_count, and reads
    # This is to prioritize higher rank taxa, taxa that are more prevalent, and taxa that have more reads
    sorted_group = group[group['ncbi_taxa'] != main_level_taxid].sort_values(['rank','sample_count','reads'], ascending=False)
    
    # Find entries in library report matching this main_level_taxid
    matching_entries = microbiome_library_report[microbiome_library_report['taxid'].isin(sorted_group['ncbi_taxa'])]

    # sometime only main_level_taxid is in the library report
    if matching_entries.empty:
        try:
            # only use the main_level_taxid
            matching_entries_accession = microbiome_library_report[microbiome_library_report['taxid'] == main_level_taxid]['URL'].unique()[:max_urls_per_taxid]
            taxids_to_download[main_level_taxid].update(matching_entries_accession)
        except Exception as e:
            logger.error(f"Error in matching entries for {main_level_taxid}: {e}")
            sys.exit(1)
    else:

        # First, ensure each ncbi_taxa gets one URL
        for taxid in sorted_group['ncbi_taxa'].astype(str):
            if len(taxids_to_download[main_level_taxid]) >= max_urls_per_taxid:
                break
                
            taxid_entries = matching_entries[matching_entries['taxid'] == taxid]
            if not taxid_entries.empty:
                # Take the first URL for this taxid
                taxid_url = taxid_entries['URL'].unique()[0]
                
                if taxid_url not in taxids_to_download[main_level_taxid]:
                    taxids_to_download[main_level_taxid].add(taxid_url)
        
        # If we still have room for more URLs, add additional ones in priority order
        remaining_slots = max_urls_per_taxid - len(taxids_to_download[main_level_taxid])
        if remaining_slots > 0:
            for taxid in sorted_group['ncbi_taxa'].astype(str):
                taxid_entries = matching_entries[matching_entries['taxid'] == taxid]
                if not taxid_entries.empty:
                    # Take the first URL for this taxid
                    taxid_url_list = taxid_entries['URL'].unique()[:remaining_slots]
                    # if the taxid_url is not in the selected_urls, add it
                    for taxid_url in taxid_url_list:
                        if taxid_url not in taxids_to_download[main_level_taxid]:
                            taxids_to_download[main_level_taxid].add(taxid_url)
    
    all_selected_entries.update(taxids_to_download[main_level_taxid])

logger.info(f"Selected {len(all_selected_entries)} entries", status='complete')
# use all_selected_entries to filter library_report
selected_library_report = library_report[library_report['URL'].isin(all_selected_entries)]
# get the interest seqids
selected_genome_set = set(selected_library_report['seqid'])
# for each URL, keep only the first row
selected_library_report_to_download = selected_library_report.drop_duplicates(subset=['URL'])
num_to_download_url = len(selected_library_report_to_download)
logger.info(f"Selected {num_to_download_url} URLs to download", status='complete')
# only keep the url with the seqid dont exist in the existing_sequences
selected_library_report_to_download = selected_library_report_to_download[~selected_library_report_to_download['seqid'].isin(existing_sequences)]

def main():

    # get the number of processors
    n_processors = int(args.proc)
    
    logger.info("Beginning to download genomes", status='run')
    # Convert DataFrame rows to dictionaries for multiprocessing
    rows_to_process = selected_library_report_to_download.to_dict('records')
    # Check if rows_to_process is empty and selected_library_report_to_download is not empty
    download_genome_status = True
    if len(rows_to_process) == 0:
        if num_to_download_url > 0:
            logger.info("All sequences have already been downloaded.", status='complete')

        else:
            logger.error("No sequences to download, please check the library report and the existing sequences", status='error')
            sys.exit(1)
    else:
        # First stage: parallel download all genome files
        n_processors = int(args.proc)
        if len(rows_to_process) < n_processors:
            n_processors = len(rows_to_process)
            
        logger.info(f"Parallel downloading {len(rows_to_process)} genomes with {n_processors} processors", status='run')
        with Pool(processes=n_processors) as pool:
            download_results = pool.map(download_genome, rows_to_process)
        
        # Filter successful downloads
        successful_downloads = [r for r in download_results if r['success']]
        logger.info(f"Successfully downloaded {len(successful_downloads)} of {len(rows_to_process)} genomes", status='complete')
        
        # Second stage: serial processing of all files for writing to library_fna
        logger.info("Processing downloaded genomes and writing to library file", status='run')
        
        # Open library file for appending
        with open(args.library_fna, 'a') as library_file:
            for result in successful_downloads:
                try:
                    fna_path = result['path']
                    seqid = result['seqid']
                    is_interest = seqid in selected_genome_set
                    
                    with gzip.open(fna_path, 'rt') as fna_file:
                        for record in SeqIO.parse(fna_file, "fasta"):
                            if record.id not in existing_sequences:
                                # write to library file
                                SeqIO.write(record, library_file, "fasta")                                
                            else:
                                logger.info(f'Skipping duplicate sequence {record.id}', status='skip')
                    
                    logger.info(f'Successfully processed {result["accession"]} fna file', status='complete')
                except Exception as e:
                    logger.error(f"Error processing {result['accession']}: {str(e)}")
    
    logger.info("Processing genomes and writing to interest file", status='run')
    # Create interest_fna file (if needed)
    interest_file = None
    if args.interest_fna:
        interest_file = open(args.interest_fna, 'w')
    
    # Open library file for read
    already_added_genomes = set()
    with open(args.library_fna, 'r') as library_file:
        for record in SeqIO.parse(library_file, "fasta"):
            if record.id in selected_genome_set and record.id not in already_added_genomes:
                SeqIO.write(record, interest_file, "fasta")
                already_added_genomes.add(record.id)
    
    # close the interest file
    if interest_file:
        interest_file.close()

    # save the accession to taxid mapping as a tsv file
    selected_library_report[['seqid','taxid']].to_csv(args.acc2tax, index=False, sep="\t")
    
    # TODO: the error cutoff select
    if len(already_added_genomes)/len(selected_genome_set) != 1:
        logger.warning(f"Only {len(already_added_genomes)/len(selected_genome_set) * 100}% of the selected genomes were added to the library file", status='check')
        sys.exit(1)
    else:
        logger.info(f"Successfully processed {len(already_added_genomes)} genomes, completed ratio {len(already_added_genomes)/len(selected_genome_set) * 100}%", status='complete')

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    main()