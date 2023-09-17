import pysam
import logging
import os
import argparse
import pandas as pd 
import re
import multiprocessing as mp
from collections import defaultdict
import mmap
import collections
import numpy as np
from collections import Counter
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, Manager
import itertools
import logging
import sys
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


# def calc_correlation(group):
#     results = {}
#     if len(group['kmer_counts'].unique()) > 1 and len(group['global_unique_kmer_counts'].unique()) > 1:
#         corr, p_value = spearmanr(group['kmer_counts'], group['global_unique_kmer_counts'])
#         results['corr_kmer_global_uniq_counts'], results['p_value_kmer_global_uniq_counts'] = corr, p_value

#     if len(group['kmer_counts'].unique()) > 1 and len(group['unique_kmer_counts'].unique()) > 1:
#         corr, p_value = spearmanr(group['kmer_counts'], group['unique_kmer_counts'])
#         results['corr_kmer_uniq_counts'], results['p_value_kmer_uniq_counts'] = corr, p_value
        
#     return pd.DataFrame([results])


def calc_correlation(group):
    results = {}
    if len(group['kmer_counts']) > 1 and len(group['global_unique_kmer_counts']) > 1:
        if len(group['kmer_counts'].unique()) == 1 and len(group['global_unique_kmer_counts'].unique()) == 1:
            corr, p_value = 1,1e-10 # since zero will make error, use a very small number
        else:
            corr, p_value = spearmanr(group['kmer_counts'], group['global_unique_kmer_counts'])
        results['corr_kmer_global_uniq_counts'], results['p_value_kmer_global_uniq_counts'] = corr, p_value

    if len(group['kmer_counts']) > 1 and len(group['unique_kmer_counts']) > 1:
        if len(group['kmer_counts'].unique()) == 1 and len(group['unique_kmer_counts'].unique()) == 1:
            corr, p_value = 1,1e-10
        else:
            corr, p_value = spearmanr(group['kmer_counts'], group['unique_kmer_counts'])
        results['corr_kmer_uniq_counts'], results['p_value_kmer_uniq_counts'] = corr, p_value
        
    return pd.DataFrame([results])


def process_taxa_kmer_dict(taxa_name, child_parent, taxid_rank, kmer_map_dict, kmer_len, min_frac):
    species_taxid = taxid_to_desired_rank(str(taxa_name), 'species', child_parent, taxid_rank)
    genus_taxid = taxid_to_desired_rank(str(taxa_name), 'genus', child_parent, taxid_rank)
    family_taxid = taxid_to_desired_rank(str(taxa_name), 'family', child_parent, taxid_rank)
        

    matching_indices = []

    for i, kread in kmer_map_dict.items():

        main_taxid = kread['main_taxid']

        if taxa_name == main_taxid or species_taxid == main_taxid or genus_taxid == main_taxid or family_taxid == main_taxid:
            matching_indices.append(i)
    
    kmer_df_batch = []  # Temporary list to store data for this batch

    for i in matching_indices:
        kread = kmer_map_dict[i]
        query_name = kread['query_name']
        total_kmer_count = np.sum(kread['kmer_counts'])
        kread_main_taxid = kread['main_taxid']
        selected_mask = np.isin(kread['taxids'], [0, taxa_name, species_taxid, genus_taxid,family_taxid])
        selected_kmer_counts = kread['kmer_counts'][selected_mask]
        selected_kmer_count = np.sum(selected_kmer_counts)

        # Calculate the percentage of selected k-mer counts out of total k-mer counts
        selected_percentage = selected_kmer_count / total_kmer_count
        temp_kmer_df = []  # Temporary list to store data for this iteration
        start_pos = 1
        # If the selected percentage is less than min_frac, skip adding to kmer_df_batch
        if selected_percentage >= min_frac['value']:
            for j, taxid in enumerate(kread['taxids']):
                kmer_count = kread['kmer_counts'][j]
                end_pos = start_pos + kmer_len + kmer_count - 1 - 1
                if taxid in (taxa_name, species_taxid):
                    temp_kmer_df.append([query_name, taxid, start_pos, end_pos, kmer_count,species_taxid,kread_main_taxid])
                start_pos += kmer_count - 1
            # Extend the kmer_df_batch list with the data from temp_kmer_df
            kmer_df_batch.extend(temp_kmer_df)

    return kmer_df_batch

# 计算k-mer一致性
def kmer_consistency(sequence, k=6):
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    return len(kmer_counts) / len(kmers)

# 计算DUST得分
def dust_score(sequence):
    total_length = len(sequence)
    unique_chars = set(sequence)
    num_unique_chars = len(unique_chars)

    if num_unique_chars == 0:
        return 0

    frequency = {}
    for char in sequence:
        if char in frequency:
            frequency[char] += 1
        else:
            frequency[char] = 1

    dust_score = num_unique_chars / total_length
    return dust_score

# 计算信息熵
def calculate_entropy(sequence):
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    sequence_length = len(sequence)
    
    for nucleotide in sequence:
        if nucleotide in nucleotide_counts:
            nucleotide_counts[nucleotide] += 1
    
    nucleotide_probabilities = [count / sequence_length for count in nucleotide_counts.values()]
    nucleotide_probabilities = [p for p in nucleotide_probabilities if p > 0]  # 只对非零概率的核苷酸进行计算
    entropy = -np.sum(nucleotide_probabilities * np.log2(nucleotide_probabilities))
    
    return entropy

def make_dicts(nodes_file):
    with open(nodes_file, 'r') as infile:
        # make a child to parent dictionary
        # and a taxid to rank dictionary
        child_to_parent = {}
        taxid_to_rank = {}
        for line in infile:
            line=line.rstrip('\n').split('\t')
            child, parent, rank = line[0], line[2], line[4]
            child_to_parent[child] = parent
            taxid_to_rank[child] = rank
    return child_to_parent, taxid_to_rank


def taxid_to_desired_rank(taxid, desired_rank, child_parent, taxid_rank):
    # look up the specific taxid,
    # build the lineage using the dictionaries
    # stop at the desired rank and return the taxid
    lineage = [[taxid, taxid_rank[taxid]]]
    if taxid_rank[taxid] == desired_rank:
        return taxid
    child, parent = taxid, None
    if child == '0':
        return 'unclassified'
    while not parent == '1':
        # print(child, parent)
        # look up child, add to lineage
        parent = child_parent[child]
        rank = taxid_rank[parent]
        lineage.append([parent, rank])
        if rank == desired_rank:
            return parent
        child = parent # needed for recursion
    return 'error - taxid above desired rank, or not annotated at desired rank'

def testFilesCorrespondingReads(inputfile_krakenAlign, inputfile_unmappedreads,numberLinesToTest=500):
    lines_tested = 0
    kraken_query_names = set(inputfile_krakenAlign['query_name'])  # Assuming 'query_name' is the column containing read names in inputfile_krakenAlign
    
    with pysam.AlignmentFile(inputfile_unmappedreads, "rb") as bam_file:
        for sread in bam_file:
            # 检查read的query_name是否在Kraken的DataFrame中
            if sread.query_name not in kraken_query_names:
                print("ERROR: corresponding test failed for files:", inputfile_krakenAlign, "and", inputfile_unmappedreads)
                return False
            
            lines_tested += 1
            if lines_tested >= numberLinesToTest:
                break

    return True


#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--krak_report', required=True, 
        dest="krak_report_file", help='Input kraken report file for denosing')
    parser.add_argument('--krak_output', required=True,
        dest='krak_output_file', help='Input kraken output file for denosing')
    parser.add_argument('--krak_mpa_report', required=True,
        dest='krak_mpa_report_file', help='Input kraken output file for denosing')
    parser.add_argument('--bam', required=True,
        dest='bam_file', help='Input origin bam file for denosing')
    parser.add_argument('--raw_qc_output_file', required=True,
        help='Output denosed info at individual level')
    parser.add_argument('--qc_output_file', required=True,
        help='Output denosed info at individual level')
    parser.add_argument('--nodes_dump', required=True,
        help='Kraken2 database node tree file path')
    parser.add_argument('--inspect', required=True,
        dest="inspect_file", help='Kraken2 database inspect file path')
    parser.add_argument('--kmer_len', required=False,
        default=35, help='Kraken classifer kmer length [default=35]')
    parser.add_argument('--num_processes', required=False,
        default=8, help='Cores use [default=8]')
    parser.add_argument('--exclude', required=False,
        default=9606, nargs='+',
        help='Taxonomy ID[s] of reads to exclude (space-delimited)')
    parser.add_argument('--min_frac', required=False,
        default=0.5, type=float, help='minimum fraction of kmers directly assigned to taxid to use read [default=0.5]')
    parser.add_argument('--min_entropy', required=False,
        default=1.2, type=float, help='minimum entropy of sequences cutoff [default=1.2]')
    parser.add_argument('--min_dust', required=False,
        default=0.1, type=float, help='minimum dust score of sequences cutoff [default=1.2]')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
    parser.add_argument("--barcode_tag", default="CB", help="Barcode tag to use for extracting barcodes")

    args=parser.parse_args()
    
    # Set log level based on command line arguments
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Create a file handler and add the formatter to it
    file_handler = logging.FileHandler(args.log_file)  # Output logs to the specified file
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    logger.info('Parsing taxonmy full lineage infomation from NCBI nodes.dump', status='run')
    try:
        child_parent, taxid_rank = make_dicts(args.nodes_dump)
        logger.info('Successfully parsing taxonmy full lineage infomation from NCBI nodes.dump', status='complete')
    except:
        logger.error("Couldn't get the taxonmy full lineage infomation from NCBI nodes.dump")
        sys.exit()

    logger.info('Reading kraken2 classifier result infomation from report', status='run')
    krak_report = pd.read_csv(args.krak_report_file, sep="\t", names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
    # remove space
    krak_report['scientific name'] = krak_report['scientific name'].str.strip() 
    # replace space
    krak_report['scientific name'] = krak_report['scientific name'].str.replace(r' ', '_')
    logger.info('Finishing reading kraken2 classifier result infomation from report', status='complete')
    logger.info('Reading kraken2 database minimizers from inspect txt', status='run')
    krak2_inspect = pd.read_csv(args.inspect_file, sep="\t", names=['frac','minimizers_clade', 'minimizers_taxa', 'rank','ncbi_taxonomy','sci_name'])

    krak_report = krak_report.merge(krak2_inspect[['ncbi_taxonomy', 'minimizers_taxa', 'minimizers_clade']],
                                left_on='ncbi_taxa',
                                right_on='ncbi_taxonomy',
                                how='left')

    krak_report.drop(columns='ncbi_taxonomy', inplace=True)
    krak_report['cov'] = krak_report['uniqminimizers']/krak_report['minimizers_taxa']
    krak_report['dup'] = krak_report['minimizers']/krak_report['uniqminimizers']

    bacteria_mean_cov =(krak_report[krak_report['scientific name'].str.strip().isin(['Bacteria'])]['minimizers']/krak_report[krak_report['scientific name'].str.strip().isin(['Bacteria'])]['minimizers_clade']).sum()
    archaea_mean_cov =(krak_report[krak_report['scientific name'].str.strip().isin(['Archaea'])]['minimizers']/krak_report[krak_report['scientific name'].str.strip().isin(['Archaea'])]['minimizers_clade']).sum()
    microbiome_reads = krak_report[krak_report['scientific name'].isin(['Bacteria', 'Fungi', 'Viruses'])]['fragments'].sum()
    total_reads = krak_report['fragments'].iloc[0] + krak_report['fragments'].iloc[1]
    logger.info('Reading kraken2 bacteria, fungi, virus classifier rank infomation from mpa report', status='run')
    krak_mpa_report = pd.read_csv(args.krak_mpa_report_file, sep='\t', names=['mpa_taxa','reads'])
    krak_mpa_report['taxa'] = krak_mpa_report['mpa_taxa'].apply(lambda x: re.sub(r'[a-z]__', '', x.split('|')[-1]))

    # we only focus on  k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea
    keywords = ["k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea"]

    krak_mpa_report_subset = krak_mpa_report[krak_mpa_report['mpa_taxa'].str.contains('|'.join(keywords))]
    df = pd.merge(krak_mpa_report_subset, krak_report, left_on='taxa', right_on='scientific name')
    # filter kraken_file to species and subspecies only
    desired_krak_report = krak_report.copy()[krak_report['classification_rank'].str.startswith(('S'), na=False)]
    desired_krak_report['species_level_taxid'] = desired_krak_report.apply(lambda x: taxid_to_desired_rank(str(x['ncbi_taxa']), 'species', child_parent, taxid_rank), axis=1)

    desired_krak_report = desired_krak_report[desired_krak_report['scientific name'].isin(df['scientific name'])] 

    logger.info('Finished processing kraken2 classifier result', status='complete')
    del df
    # Reading kraken2 classifier output information
    logger.info('Reading kraken2 classifier output information', status='run')
    krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type','query_name', 'taxid_info', 'len','kmer_position'])
    krak2_output[['taxa', 'taxid']] = krak2_output['taxid_info'].str.extract(r'(.*) \(taxid (\d+)\)')

    # Remove the ')' and leading/trailing whitespace from the 'taxid' and 'name' columns
    krak2_output['taxid'] = krak2_output['taxid'].str.replace(r'\)', '').str.strip()
    krak2_output['r1_len'], krak2_output['r2_len'] = zip(*krak2_output['len'].apply(lambda x: x.split('|')))
    krak2_output['r1_kmer_position'], krak2_output['r2_kmer_position']  = zip(*krak2_output['kmer_position'].apply(lambda x: x.split('|:|')))
    r1_df = krak2_output[['query_name', 'taxa', 'taxid', 'r1_len', 'r1_kmer_position']] 
    r2_df = krak2_output[['query_name', 'taxa', 'taxid', 'r2_len', 'r2_kmer_position']]
    # since pair end read maybe less than k-mer,so kraken could not classified it, filter less than kmer
    r1_df = r1_df[r1_df['r1_len'] >= str(args.kmer_len)]
    r2_df = r2_df[r2_df['r2_len'] >= str(args.kmer_len)]

    # Remove krak2_output
    del krak2_output 
    logger.info('Finishing reading kraken2 classifier output information', status='complete')
    logger.info('Pre-calculate and store taxid, kmer_count, and taxid index information', status='run')
    
    # Pre-calculate and store taxid, kmer_count, and taxid index information in a dictionary
    kmer_map_dict_r1 = {}
    for i, row in r1_df.iterrows():
        kmer_position = row['r1_kmer_position'].strip()
        kmer_info = kmer_position.split(" ")
        kmer_map_dict_r1[i] = {
            'taxids': np.array([info.split(":")[0] for info in kmer_info]),
            'kmer_counts': np.array([int(info.split(":")[1]) for info in kmer_info]),
            "main_taxid": row['taxid'],
            'query_name': row['query_name'],
        }
    # Pre-calculate and store taxid, kmer_count, and taxid index information in a dictionary
    kmer_map_dict_r2 = {}
    for i, row in r2_df.iterrows():
        kmer_position = row['r2_kmer_position'].strip()
        kmer_info = kmer_position.split(" ")
        kmer_map_dict_r2[i] = {
            'taxids': np.array([info.split(":")[0] for info in kmer_info]),
            'kmer_counts': np.array([int(info.split(":")[1]) for info in kmer_info]),
            "main_taxid": row['taxid'],
            'query_name': row['query_name'],
        }
    logger.info(f'Parsing taxid, kmer_count, and taxid index information with {args.num_processes} processes', status='run')


    # desired_krak_report['ncbi_taxa'] = desired_krak_report['ncbi_taxa'].astype(str)
    desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    desired_taxid_list = set(desired_krak_report['ncbi_taxa'].unique())
    # Specify the number of processes to use (e.g., 4)
    num_processes = int(args.num_processes)
    manager = Manager()
    min_frac = manager.dict()
    min_frac['value'] = float(args.min_frac)
    # Use multiprocessing to parallelize the execution for each taxa_name
    with Pool(processes=num_processes) as pool:
        # Pass additional arguments to process_taxa_kmer_dict
        kmer_df_batches_r1 = pool.starmap(process_taxa_kmer_dict, [(taxa_name, child_parent, taxid_rank, kmer_map_dict_r1, args.kmer_len, min_frac) for taxa_name in desired_taxid_list])
    # Flatten the list of batches
    kmer_df_r1 = [item for sublist in kmer_df_batches_r1 for item in sublist]

    kraken_df_r1 = pd.DataFrame(kmer_df_r1, columns=[
        "query_name", "taxid", "start_pos", "end_pos", "kmer_count", "species_level_taxid", "sequence_taxid"
    ])

    kraken_df_r1 = kraken_df_r1.drop_duplicates(subset=['query_name', 'start_pos', 'end_pos'], keep='first')

    # 将Kraken的DataFrame的query_name列转换为一个集合
    kraken_r1_query_names = set(kraken_df_r1["query_name"])
    # Use multiprocessing to parallelize the execution for each taxa_name
    with Pool(processes=num_processes) as pool:
        # Pass additional arguments to process_taxa_kmer_dict
        kmer_df_batches_r2 = pool.starmap(process_taxa_kmer_dict, [(taxa_name, child_parent, taxid_rank, kmer_map_dict_r2, args.kmer_len, min_frac) for taxa_name in desired_taxid_list])

    # Flatten the list of batches
    kmer_df_r2 = [item for sublist in kmer_df_batches_r2 for item in sublist]

    kraken_df_r2 = pd.DataFrame(kmer_df_r2, columns=[
        "query_name", "taxid", "start_pos", "end_pos", "kmer_count", "species_level_taxid", "sequence_taxid"
    ])

    kraken_df_r2 = kraken_df_r2.drop_duplicates(subset=['query_name', 'start_pos', 'end_pos'], keep='first')
    
    # 将Kraken的DataFrame的query_name列转换为一个集合
    kraken_r2_query_names = set(kraken_df_r2["query_name"])

    logger.info(f'Finished parsing taxid, kmer_count, and taxid index information', status='complete')

    combined_query_names = kraken_r1_query_names | kraken_r2_query_names
    # Get species level taxid
    # kraken_df['species_level_taxid'] = kraken_df.apply(lambda x: taxid_to_desired_rank(str(x['taxid']), 'species', child_parent, taxid_rank), axis=1)
    # kraken_df['species_level_taxid'] = kraken_df.apply(lambda x: taxid_to_desired_rank(str(x['taxid']), 'species', child_parent, taxid_rank)['species_taxid'], axis=1)
    num_unique_species = len(kraken_df_r1['species_level_taxid'].unique())

    logger.info(f'Found {num_unique_species} unique species level taxids', status='summary')

    logger.info(f'Get the raw classified reads from bam file', status='run')

    # 将DataFrame转换为字典
    kraken_dict_r1 = kraken_df_r1.groupby('query_name').apply(lambda x: x.to_dict('records')).to_dict()
    # 将DataFrame转换为字典
    kraken_dict_r2 = kraken_df_r2.groupby('query_name').apply(lambda x: x.to_dict('records')).to_dict()
    # Since we have kraken_dict , we can del kraken_df to release memory
    del kraken_df_r1
    del kraken_df_r2

    skipped = 0
    # Move constant values outside the loop
    kmer_key = "kmer_count"
    start_pos_key = "start_pos"
    end_pos_key = "end_pos"
    species_level_taxid_key = "species_level_taxid"
    # Init bam count
    read_count = 0
    use_count = 0
    krak_count = 0
    # # Create a defaultdict to map species_level_taxid to its set of all kmers
    # taxid_to_kmers = defaultdict(set)
    # Create a dictionary to map CB and taxid to its set of all UB and kmers
    cb_taxid_to_ub_kmers = defaultdict(lambda: {"kmers": []})  # Using a nested defaultdict
    cb_metrics_dict = defaultdict(list)
    logger.info(f'Parsing the raw classified reads from bam file', status='run')
    # 遍历bam文件的每一个read
    for sread in pysam.AlignmentFile(args.bam_file, "rb"):
        read_count += 1
        # 检查read的query_name是否在Kraken的DataFrame中
        if sread.query_name not in combined_query_names:
            skipped += 1
            continue

        krak_count += 1
        # 获取CB和UB
        try:
            sread_CB = sread.get_tag(args.barcode_tag)
            # sread_UB = sread.get_tag('UB')
            use_count += 1
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

        if sread_CB == '-':
            skipped += 1 
            continue
        
        # 获取对应的Kraken信息
        if sread.is_read1:
            try:
                kraken_info = kraken_dict_r1.get(sread.query_name, [])
            except:
                continue
        # 获取对应的Kraken信息
        if sread.is_read2:
            try:
                kraken_info = kraken_dict_r2.get(sread.query_name, [])
            except:
                continue
        
        # 遍历每一个Kraken信息
        for row in kraken_info:
            # 获取序列的子串
            sequence = sread.seq[row[start_pos_key] - 1:row[end_pos_key]]

            # Check sequence length and skip if necessary
            if len(sequence) < args.kmer_len:
                print(f"Read: {sread.query_name}, start_pos: {row[start_pos_key]}, end_pos: {row[end_pos_key]}, sequence: {sequence}")
                print("Sequence length less than kmer length")
                print(len(sread.seq))
                continue

            # 计算指标
            seq_kmer_consistency = kmer_consistency(sequence)
            seq_entropy = calculate_entropy(sequence)
            seq_dust_score = dust_score(sequence)
            seq_length = len(sequence)

            # Use a sliding window of 1 to get all kmers
            kmers = (sequence[i:i+args.kmer_len] for i in range(seq_length - args.kmer_len + 1))
            # Limit the number of kmers to kmer_count
            kmers = list(itertools.islice(kmers, row[kmer_key]))


            # Store the result in the dictionary
            cb_metrics_dict[sread.query_name].append([
                sread_CB, row[species_level_taxid_key], seq_kmer_consistency,
                seq_entropy, seq_dust_score, seq_length
            ])
            
            key = (sread_CB, row[species_level_taxid_key])
            # cb_taxid_to_ub_kmers[key]["UB"].append(sread_UB)
            cb_taxid_to_ub_kmers[key]["kmers"].extend(kmers)

    logger.info(f'Finished parsing the raw classified reads from bam file', status='run')
    logger.info(f'Total unmapped reads: {read_count}', status='summary')
    logger.info(f'Total unmapped reads classified by Kraken: {krak_count}', status='summary')
    logger.info(f'Total classified Reads with CB : {use_count}', status='summary')
    logger.info(f'Skipped reads: {skipped}', status='summary')
    logger.info(f'Finishing getting the raw classified reads from bam file', status='complete')
    logger.info(f'Calculating quality control indicators', status='run')
    # 将 cb_metrics_dict 转换为 DataFrame
    cb_metrics_list = []
    for cb_data_list in cb_metrics_dict.values():
        cb_metrics_list.extend(cb_data_list)

    metrics_columns = [
        'CB', 'species_level_taxid', 'seq_kmer_consistency',
        'seq_entropy', 'seq_dust_score', 'seq_length'
    ]
    cb_metrics_df = pd.DataFrame(cb_metrics_list, columns=metrics_columns)
    del cb_metrics_dict

    # 根据 species 计算统计值
    species_seq_metrics = cb_metrics_df.groupby('species_level_taxid').agg({
        'seq_kmer_consistency': 'mean',
        'seq_entropy': 'mean',
        'seq_dust_score': 'mean',
        'seq_length': ['max', 'mean']
    }).reset_index()

    # 重命名统计值的列名
    species_seq_metrics.columns = [
        'species_level_taxid', 'average_kmer_consistency',
        'average_seq_entropy', 'average_seq_dust_score',
        'max_seq_length', 'mean_seq_length'
    ]

    # Create a list of dictionaries for
    data = [{"CB": cb, "species_level_taxid": species_level_taxid, "kmers": kmers["kmers"]} 
            for (cb, species_level_taxid), kmers in cb_taxid_to_ub_kmers.items()]

    # Create the DataFrame from the list of dictionaries
    cb_taxid_ub_kmer_count_df = pd.DataFrame(data)
    # Del data
    del data
    del cb_taxid_to_ub_kmers

    num_unique_CB = len(cb_taxid_ub_kmer_count_df['CB'].unique())
    if num_unique_CB > 1:
        # Convert the DataFrame to long format, each row contains a kmer
        cb_taxid_ub_kmer_count_df = cb_taxid_ub_kmer_count_df.explode('kmers')

        # Calculate total kmer counts for each CB and species_level_taxid combination
        total_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'species_level_taxid']).size().reset_index(name='kmer_counts')

        # Calculate number of unique kmers for each CB and species_level_taxid combination 
        unique_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'species_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'unique_kmer_counts'})


        # 标识重复的 kmers，并获取不重复的行
        cb_taxid_ub_global_unique_count_df = cb_taxid_ub_kmer_count_df[~cb_taxid_ub_kmer_count_df.duplicated(subset=['kmers'], keep=False)]

        global_unique_kmer_counts =cb_taxid_ub_global_unique_count_df.groupby(['CB', 'species_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'global_unique_kmer_counts'})

        cb_taxid_kmer_count_df = pd.merge(total_kmer_counts, unique_kmer_counts, on=['CB', 'species_level_taxid'])
        cb_taxid_kmer_count_df = pd.merge(cb_taxid_kmer_count_df,global_unique_kmer_counts, on=['CB', 'species_level_taxid'])
        del cb_taxid_ub_global_unique_count_df
        del cb_taxid_ub_kmer_count_df
        cb_taxid_kmer_corr_df = cb_taxid_kmer_count_df.groupby('species_level_taxid').apply(calc_correlation).reset_index()

        p_val_cols = ['p_value_kmer_global_uniq_counts', 'p_value_kmer_uniq_counts']

        non_na_rows = cb_taxid_kmer_corr_df[p_val_cols].notna().any(axis=1)

        # Calculate ntests using non-na rows
        ntests = non_na_rows.sum()
        # Perform multiple testing correction only if ntests is non-zero
        if ntests > 5:
            # Filter out NaN p-values and adjust them for both sets
            for col in p_val_cols:
                if col in cb_taxid_kmer_corr_df.columns:
                    non_nan = cb_taxid_kmer_corr_df[col].notna()
                    cb_taxid_kmer_corr_df.loc[non_nan, col] = multipletests(cb_taxid_kmer_corr_df.loc[non_nan, col], method='fdr_bh')[1]
        else:
            pass
    else:
        ntests = 0
        pass

    final_desired_krak_report = desired_krak_report.copy()
    # Convert 'ncbi_taxa' column to string data type
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    # final_desired_krak_report.drop('fraction', axis=1, inplace=True)
    final_desired_krak_report['cov'].replace([float('inf'), float('-inf')], float('nan'), inplace=True)
    final_desired_krak_report['max_cov'] = final_desired_krak_report.groupby('species_level_taxid')['cov'].transform('max')
    final_desired_krak_report['max_minimizers'] = final_desired_krak_report.groupby('species_level_taxid')['minimizers'].transform('max')
    final_desired_krak_report['max_uniqminimizers'] = final_desired_krak_report.groupby('species_level_taxid')['uniqminimizers'].transform('max')
    final_desired_krak_report["rpm"] =final_desired_krak_report['fragments']/ total_reads * 1e6
    final_desired_krak_report["rpmm"] =final_desired_krak_report['fragments']/ microbiome_reads * 1e6
    final_desired_krak_report = final_desired_krak_report.merge(species_seq_metrics,left_on='species_level_taxid', right_on='species_level_taxid')
    # final_desired_krak_report.drop('taxid', axis=1, inplace=True)
    if num_unique_CB > 1:
        final_desired_krak_report = final_desired_krak_report.merge(cb_taxid_kmer_corr_df,left_on='species_level_taxid', right_on='species_level_taxid')
    # final_desired_krak_report.drop('species_level_taxid', axis=1, inplace=True)
    else:
        pass
    
    final_desired_krak_report['superkingdom'] = final_desired_krak_report.apply(lambda x: taxid_to_desired_rank(str(x['ncbi_taxa']),'superkingdom', child_parent, taxid_rank), axis=1)
    final_desired_krak_report['bacteria_mean_cov_cutoff'] = bacteria_mean_cov/50
    final_desired_krak_report['archaea_mean_cov_cutoff'] = archaea_mean_cov/50

    logger.info(f'Finishging calculating quality control indicators', status='complete')

    num_unique_species = len(final_desired_krak_report['ncbi_taxa'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids having qc indictor', status='summary')

    # Save data
    logger.info(f'Saving the raw result', status='run')
    final_desired_krak_report.to_csv(args.raw_qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')

    logger.info(f'Filtering taxa with quality control indicators', status='run')
    final_desired_krak_report['superkingdom'] = final_desired_krak_report['superkingdom'].astype(str)
    ## For many corr
    if ntests > 5 and num_unique_CB > 1:
        filter_desired_krak_report = final_desired_krak_report.copy()[
            (final_desired_krak_report['average_seq_entropy'] > args.min_entropy) &
            (final_desired_krak_report['max_minimizers'] > 5) &
            (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
            (
                (
                    ((final_desired_krak_report['superkingdom'] == '2') & (final_desired_krak_report['max_cov']*50 >=  bacteria_mean_cov)) |
                    ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov']*50 >=  archaea_mean_cov)) |
                    ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >  0)) |
                    ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0)) 
                )
            ) &
            (
                (
                    (final_desired_krak_report['corr_kmer_uniq_counts'] > 0.5) &
                ((final_desired_krak_report['p_value_kmer_uniq_counts'] < float(0.05)) | (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()))
                ) |
                (
                    (final_desired_krak_report['corr_kmer_uniq_counts'].isna()) &
                    (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()) &
                    (final_desired_krak_report['average_seq_entropy'] > 1.5*args.min_entropy) &
                    (final_desired_krak_report['max_seq_length'] > float(38))
                )
            )
        ]
    else:
        filter_desired_krak_report = final_desired_krak_report.copy()[
            (final_desired_krak_report['average_seq_entropy'] > args.min_entropy) &
            (final_desired_krak_report['max_minimizers'] > 5) &
            (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
            (
                (
                    ((final_desired_krak_report['superkingdom'] == '2') & (final_desired_krak_report['max_cov']*50 >=  bacteria_mean_cov)) |
                    ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov']*50 >=  archaea_mean_cov)) |
                    ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >  0)) |
                    ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0)) 
                )
            ) 
        ]
    # filter_desired_krak_report.drop(['frac','classification_rank','fraction','minimizers_clade','minimizers_taxa','ncbi_taxa','sci_name','cov','species_level_taxa','level_1'], axis=1, inplace=True)
    filter_desired_krak_report['scientific name'] = filter_desired_krak_report['scientific name'].apply(lambda x: x.strip())
    
    # # Filter out rows where 'ncbi_taxa' matches any value from 'excluded_taxonomy_ids'
    # filter_desired_krak_report = filter_desired_krak_report[~filter_desired_krak_report['ncbi_taxa'].isin(args.exclude)]

    logger.info(f'Finishing filtering taxa with quality control indicators', status='complete')
    num_unique_species = len(filter_desired_krak_report['ncbi_taxa'].unique())
    logger.info(f'After filtering, found {num_unique_species} unique species and subspeceis level taxids', status='summary')
    # Save data
    logger.info(f'Saving the result', status='run')
    filter_desired_krak_report.to_csv(args.qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')


if __name__ == "__main__":
    main()