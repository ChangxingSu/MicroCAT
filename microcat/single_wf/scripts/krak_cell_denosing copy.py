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

# Create a logger object
logger = logging.getLogger('my_logger')

# Create a formatter object
log_format = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

# Create a handler and add the formatter to it
console_handler = logging.StreamHandler()  # Output logs to the console
console_handler.setFormatter(log_format)

# Add the handler to the logger object
logger.addHandler(console_handler)

# Customize logger.info function to include status
def custom_info(msg, status=None):
    if status:
        logger.info('%s (%s)', msg, status)  # Wrap status in parentheses
    else:
        logger.info(msg)

# Bind the custom_info function to the logger object
logger.info = custom_info

def calc_correlation(group):
    results = {}
    if len(group['unique_kmer_counts'].unique()) > 1 and len(group['microbiome_ub_counts'].unique()) > 1:
        corr, p_value = spearmanr(group['unique_kmer_counts'], group['microbiome_ub_counts'])
        results['corr_ub_counts'], results['p_value_ub_counts'] = corr, p_value

    if len(group['kmer_counts'].unique()) > 1 and len(group['unique_kmer_counts'].unique()) > 1:
        corr, p_value = spearmanr(group['kmer_counts'], group['unique_kmer_counts'])
        results['corr_kmer_counts'], results['p_value_kmer_counts'] = corr, p_value
        
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
    parser.add_argument('--bam', required=True,
        dest='bam_file', help='Input origin bam file for denosing')
    parser.add_argument('--output_file', required=True,
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
        default=0.5, help='minimum fraction of kmers directly assigned to taxid to use read [default=0.5]')
    parser.add_argument('--min_entropy', required=False,
        default=1.2, help='minimum entropy of sequences cutoff [default=1.2]')
    parser.add_argument('--min_dust', required=False,
        default=0.1, help='minimum dust score of sequences cutoff [default=1.2]')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')


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

    logger.info('Reading kraken2 classifier result infomation from ', status='run')
    krak_report = pd.read_csv(args.krak_report_file, sep="\t", names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
    krak2_inspect = pd.read_csv(args.inspect_file, sep="\t", names=['frac','minimizers_clade', 'minimizers_taxa', 'rank','ncbi_taxonomy','sci_name'])

    krak_report = krak_report.merge(krak2_inspect[['ncbi_taxonomy', 'minimizers_taxa', 'minimizers_clade']],
                                left_on='ncbi_taxa',
                                right_on='ncbi_taxonomy',
                                how='left')
    krak_report.drop(columns='ncbi_taxonomy', inplace=True)
    krak_report['cov'] = krak_report['uniqminimizers']/krak_report['minimizers_taxa']
    krak_report['dup'] = krak_report['minimizers']/krak_report['uniqminimizers']

    # filter kraken_file to species only
    desired_krak_report = krak_report.copy()[krak_report['classification_rank'].str.startswith(('S'), na=False)]
    desired_krak_report['species_level_taxa'] = desired_krak_report.apply(lambda x: taxid_to_desired_rank(str(x['ncbi_taxa']), 'species', child_parent, taxid_rank), axis=1)

    logger.info('Finished processing kraken2 classifier result', status='complete')

    # Reading kraken2 classifier output information
    logger.info('Reading kraken2 classifier output information', status='run')
    krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type','query_name', 'taxid_info', 'len','kmer_position'])
    krak2_output[['taxa', 'taxid']] = krak2_output['taxid_info'].str.extract(r'(.*) \(taxid (\d+)\)')

    # Remove the ')' and leading/trailing whitespace from the 'taxid' and 'name' columns
    krak2_output['taxid'] = krak2_output['taxid'].str.replace(r'\)', '').str.strip()

    logger.info('Pre-calculate and store taxid, kmer_count, and taxid index information', status='run')
    # Pre-calculate and store taxid, kmer_count, and taxid index information in a dictionary
    kmer_map_dict = {}
    for i, row in krak2_output.iterrows():
        kmer_position = row['kmer_position'].strip()
        kmer_info = kmer_position.split(" ")
        kmer_map_dict[i] = {
            'taxids': np.array([info.split(":")[0] for info in kmer_info]),
            'kmer_counts': np.array([int(info.split(":")[1]) for info in kmer_info]),
            "main_taxid": row['taxid'],
            'query_name': row['query_name'],
        }

    logger.info(f'Parsing taxid, kmer_count, and taxid index information with {args.num_processes} processes', status='run')
    kmer_df =[]
    desired_krak_report['ncbi_taxa'] = desired_krak_report['ncbi_taxa'].astype(str)
    desired_taxid_list = set(desired_krak_report['ncbi_taxa'].unique())
    # Specify the number of processes to use (e.g., 4)
    num_processes = args.num_processes
    manager = Manager()
    min_frac = manager.dict()
    min_frac['value'] = args.min_frac
    # Use multiprocessing to parallelize the execution for each taxa_name
    with Pool(processes=num_processes) as pool:
        # Pass additional arguments to process_taxa_kmer_dict
        kmer_df_batches = pool.starmap(process_taxa_kmer_dict, [(taxa_name, child_parent, taxid_rank, kmer_map_dict, args.kmer_len, min_frac) for taxa_name in desired_taxid_list])

    # Flatten the list of batches
    kmer_df = [item for sublist in kmer_df_batches for item in sublist]

    kraken_df = pd.DataFrame(kmer_df, columns=[
        "query_name", "taxid", "start_pos", "end_pos", "kmer_count", "species_taxid", "sequence_taxid"
    ])

    kraken_df = kraken_df.drop_duplicates(subset=['query_name', 'start_pos', 'end_pos'], keep='first')
    # Get species level taxid
    kraken_df['species_level_taxid'] = kraken_df.apply(lambda x: taxid_to_desired_rank(str(x['taxid']), 'species', child_parent, taxid_rank), axis=1)

    num_unique_species = len(kraken_df['species_level_taxid'].unique())

    logger.info(f'Finished parsing taxid, kmer_count, and taxid index information', status='complete')
    logger.info(f'Found {num_unique_species} unique species level taxids', status='summary')

    logger.info(f'Get the raw classified reads from bam file', status='run')
    # 将Kraken的DataFrame的query_name列转换为一个集合
    kraken_query_names = set(kraken_df["query_name"])

    # 将DataFrame转换为字典
    kraken_dict = kraken_df.groupby('query_name').apply(lambda x: x.to_dict('records')).to_dict()

    kmer_dict = defaultdict(list)
    skipped = 0
    # Move constant values outside the loop
    kmer_key = "kmer_count"
    start_pos_key = "start_pos"
    end_pos_key = "end_pos"
    species_level_taxid_key = "species_taxid"

    # 遍历bam文件的每一个read
    for sread in pysam.AlignmentFile(args.bam_file, "rb"):
        # 检查read的query_name是否在Kraken的DataFrame中
        if sread.query_name not in kraken_query_names:
            skipped += 1
            continue

        # 获取CB和UB
        try:
            sread_CB = sread.get_tag('CB')
            sread_UB = sread.get_tag('UB')
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

        # 获取对应的Kraken信息
        kraken_info = kraken_dict.get(sread.query_name, [])

        # 遍历每一个Kraken信息
        for row in kraken_info:
            # 获取序列的子串
            sequence = sread.seq[row[start_pos_key] - 1:row[end_pos_key]]

            # Check sequence length and skip if necessary
            if len(sequence) < args.kmer_len:
                print(f"Read: {sread.query_name}, start_pos: {row[start_pos_key]}, end_pos: {row[end_pos_key]}, sequence: {sequence}")
                print(len(sread.seq))

            # Store the result in the dictionary
            kmer_dict[sread.query_name].append([
                sread_CB, sread_UB, row[species_level_taxid_key], row[start_pos_key],
                row[end_pos_key], row[kmer_key], len(sequence), sequence
            ])

    # Create a defaultdict to map species_level_taxid to its set of all kmers
    taxid_to_kmers = defaultdict(set)
    # Create a dictionary to map CB and taxid to its set of all UB and kmers
    cb_taxid_to_ub_kmers = defaultdict(lambda: {"UB": [], "kmers": []})  # Using a nested defaultdict

    # Iterate over kmer_dict
    for entries in kmer_dict.values():
        for entry in entries:
            sread_CB, sread_UB, species_level_taxid, start_pos, end_pos, kmer_count, seq_len, sequence = entry
            # Use a sliding window of 1 to get all kmers
            kmers = (sequence[i:i+args.kmer_len] for i in range(seq_len - args.kmer_len + 1))
            # Limit the number of kmers to kmer_count
            kmers = list(itertools.islice(kmers, kmer_count))
            # Add kmers to the set of the corresponding species_level_taxid
            taxid_to_kmers[species_level_taxid].update(kmers)
            # Check if taxid is 0 or "A"
            if species_level_taxid not in {"0", "A"}:
                key = (sread_CB, species_level_taxid)
                cb_taxid_to_ub_kmers[key]["UB"].append(sread_UB)
                cb_taxid_to_ub_kmers[key]["kmers"].extend(kmers)

    # Create a list of dictionaries for DataFrame
    data = [{"CB": cb, "species_level_taxid": species_level_taxid, "UB": ub_kmers["UB"], "kmers": ub_kmers["kmers"]} 
            for (cb, species_level_taxid), ub_kmers in cb_taxid_to_ub_kmers.items()]

    # Create the DataFrame from the list of dictionaries
    cb_taxid_ub_kmer_count_df = pd.DataFrame(data)

    # Convert the DataFrame to long format, each row contains a kmer
    cb_taxid_ub_kmer_count_df = cb_taxid_ub_kmer_count_df.explode('kmers').explode('UB')

    logger.info(f'Finishing getting the raw classified reads from bam file', status='complete')

    logger.info(f'Calculating quality control indicators', status='run')
    # Calculate the total kmer counts first
    kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'species_level_taxid']).size().reset_index(name='kmer_counts')

    # Calculate the unique kmer counts
    unique_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'species_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'unique_kmer_counts'})

    # Merge the two dataframes
    result_df = pd.merge(kmer_counts, unique_kmer_counts, on=['CB', 'species_level_taxid'])

    # Calculate the number of unique UBs for each CB
    ub_counts = cb_taxid_ub_kmer_count_df.groupby('CB')['UB'].nunique().reset_index(name='microbiome_ub_counts')

    # Merge the ub_counts dataframe with the result_df dataframe
    result_df = pd.merge(result_df, ub_counts, on='CB')

    correlations = result_df.groupby('species_level_taxid').apply(calc_correlation).reset_index()

    # Filter out NaN p-values and adjust them for both sets
    for col in ['p_value_ub_counts', 'p_value_kmer_counts']:
        if col in correlations.columns:
            non_nan = correlations[col].notna()
            correlations.loc[non_nan, col] = multipletests(correlations.loc[non_nan, col], method='fdr_bh')[1]

    # 创建一个空的DataFrame来存储结果
    qc_df = pd.DataFrame(columns=['taxid', 'kmer_consistency', 'entropy', 'dust_score'])

    # 根据species_taxid对字典进行分组
    grouped_dict = {}
    for entries in kmer_dict.values():
        for entry in entries:
            if entry[-6] not in grouped_dict:
                grouped_dict[entry[-6]] = []
            grouped_dict[entry[-6]].append(entry)

    # Iterate through the grouped dictionary
    for taxid, entries in grouped_dict.items():
        if taxid == 0:
            continue
        # Get all sequences
        sequences = [entry[-1] for entry in entries]
        
        # Calculate the maximum length of sequences
        max_contig = max(len(sequence) for sequence in sequences)
        # Calculate the mean length of sequences
        mean_contig = sum(len(sequence) for sequence in sequences) / len(sequences)
        # Calculate k-mer consistency, entropy, and DUST score for each sequence
        kmer_consistency_scores = [kmer_consistency(sequence) for sequence in sequences]
        entropy_scores = [calculate_entropy(sequence) for sequence in sequences]
        dust_score_values = [dust_score(sequence) for sequence in sequences]
        # Add the results to the DataFrame
        qc_df = pd.concat([qc_df, pd.DataFrame({'taxid': [taxid],
                                                'kmer_consistency': [np.mean(kmer_consistency_scores)],
                                                'entropy': [np.mean(entropy_scores)],
                                                'dust_score': [np.mean(dust_score_values)],
                                                'max_contig': [max_contig],
                                                'mean_contig': [mean_contig]})], ignore_index=True)

    final_desired_krak_report = desired_krak_report.copy()
    # Convert 'ncbi_taxa' column to string data type
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    # final_desired_krak_report.drop('fraction', axis=1, inplace=True)
    final_desired_krak_report['cov'].replace([float('inf'), float('-inf')], float('nan'), inplace=True)
    final_desired_krak_report['max_cov'] = final_desired_krak_report.groupby('species_level_taxa')['cov'].transform('max')
    final_desired_krak_report['max_minimizers'] = final_desired_krak_report.groupby('species_level_taxa')['minimizers'].transform('max')
    final_desired_krak_report['max_uniqminimizers'] = final_desired_krak_report.groupby('species_level_taxa')['uniqminimizers'].transform('max')

    final_desired_krak_report = final_desired_krak_report.merge(qc_df,left_on='ncbi_taxa', right_on='taxid')
    final_desired_krak_report.drop('taxid', axis=1, inplace=True)

    final_desired_krak_report = final_desired_krak_report.merge(correlations,left_on='ncbi_taxa', right_on='species_level_taxid')
    final_desired_krak_report.drop('species_level_taxid', axis=1, inplace=True)

    final_desired_krak_report['superkingdom'] = final_desired_krak_report.apply(lambda x: taxid_to_desired_rank(str(x['ncbi_taxa']),'superkingdom', child_parent, taxid_rank), axis=1)


    logger.info(f'Finishging calculating quality control indicators', status='complete')

    num_unique_species = len(final_desired_krak_report['species_level_taxid'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids having qc indictor', status='summary')

    logger.info(f'Filtering taxa with quality control indicators', status='run')
    filter_desired_krak_report = final_desired_krak_report.copy()[
        (final_desired_krak_report['entropy'] > args.min_entropy) &
        (final_desired_krak_report['max_minimizers'] > 5) &
        (final_desired_krak_report['dust_score'] > args.min_dust) &
        (
            (
                ((final_desired_krak_report['superkingdom'] == '2') & (final_desired_krak_report['max_cov'] > 0.000001)) |
                ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov'] >  0.000001)) |
                ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >  0)) |
                ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0)) 
            )
        ) &
        (
            (
                (final_desired_krak_report['corr_kmer_counts'] > 0.5) &
            ((final_desired_krak_report['p_value_kmer_counts'] < 0.05) | (final_desired_krak_report['p_value_kmer_counts'].isna()))
            ) |
            (
                (final_desired_krak_report['corr_kmer_counts'].isna()) &
                (final_desired_krak_report['p_value_kmer_counts'].isna()) &
                (final_desired_krak_report['entropy'] > 1.5*args.min_entropy) &
                (final_desired_krak_report['max_contig'] > 38)
            )
        )
    ]
    # filter_desired_krak_report.drop(['frac','classification_rank','fraction','minimizers_clade','minimizers_taxa','ncbi_taxa','sci_name','cov','species_level_taxa','level_1'], axis=1, inplace=True)
    filter_desired_krak_report['scientific name'] = filter_desired_krak_report['scientific name'].apply(lambda x: x.strip())
    
    # # Filter out rows where 'ncbi_taxa' matches any value from 'excluded_taxonomy_ids'
    # filter_desired_krak_report = filter_desired_krak_report[~filter_desired_krak_report['ncbi_taxa'].isin(args.exclude)]

    logger.info(f'Finishing filtering taxa with quality control indicators', status='complete')
    num_unique_species = len(filter_desired_krak_report['species_level_taxid'].unique())
    logger.info(f'After filtering, found {num_unique_species} unique species level taxids', status='summary')
    # Save data
    logger.info(f'Saving the result', status='run')
    filter_desired_krak_report.to_csv(args.output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')


if __name__ == "__main__":
    main()

