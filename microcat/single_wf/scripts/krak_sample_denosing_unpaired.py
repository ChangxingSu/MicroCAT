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

#Tree Class 
#usage: tree node used in constructing taxonomy tree  
#   includes only taxonomy levels and genomes identified in the Kraken report
class Tree(object):
    'Tree node.'
    def __init__(self,  taxid, name, level_rank, level_num, p_taxid, parent=None,children=None):
        self.taxid = taxid
        self.name = name
        self.level_rank= level_rank
        self.level_num = int(level_num)
        self.p_taxid = p_taxid
        self.all_reads = 0
        self.lvl_reads = 0
        #Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)
        
    def taxid_to_desired_rank(self, desired_rank):
        # Check if the current node's level_id matches the desired_rank
        if self.level_rank == desired_rank:
            return self.taxid
        child, parent, parent_taxid = self, None, None
        while not parent_taxid == '1':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            if rank == desired_rank:
                return parent.taxid
            child = parent # needed for recursion
        # If no parent node is found, or the desired_rank is not reached, return error
        return 'error - taxid above desired rank, or not annotated at desired rank'
    def lineage_to_desired_rank(self, desired_parent_rank):
        lineage = [] 
        lineage.append(self.taxid)
        # Check if the current node's level_id matches the desired_rank
        if self.level_num == "1":
            return lineage
        if self.level_rank == "S":
            subspecies_nodes = self.children
            while len(subspecies_nodes) > 0:
                #For this node
                curr_n = subspecies_nodes.pop()
                lineage.append(curr_n.taxid)
        child, parent, parent_taxid = self, None, None
        
        while not parent_taxid == '1':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            lineage.append(parent_taxid)
            if rank == desired_parent_rank:
                return lineage
            child = parent # needed for recursion
        return lineage

    def is_microbiome(self):
        is_microbiome = False
        main_lvls = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        lineage_name = []
        #Create level name 
        level_rank = self.level_rank
        name = self.name
        name = name.replace(' ','_')
        lineage_name.append(name)
        if level_rank not in main_lvls:
            level_rank = "x"
        elif level_rank == "K":
            level_rank = "k"
        elif level_rank == "D":
            level_rank = "d"
        child, parent, parent_taxid = self, None, None
        
        while not parent_taxid == '1':
            parent = child.parent
            level_rank = parent.level_rank
            parent_taxid = parent.taxid
            name = parent.name
            name = name.replace(' ','_')
            lineage_name.append(name)
            child = parent # needed for recursion
        if 'Fungi' in lineage_name or 'Bacteria' in lineage_name or 'Viruses' in lineage_name:
            is_microbiome = True
        return is_microbiome

    def get_mpa_path(self):
        mpa_path = []
        main_lvls = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        #Create level name 
        level_rank = self.level_rank
        name = self.name
        name = name.replace(' ','_')
        if level_rank not in main_lvls:
            level_rank = "x"
        elif level_rank == "K":
            level_rank = "k"
        elif level_rank == "D":
            level_rank = "d"
        child, parent, parent_taxid = self, None, None
        level_str = level_rank.lower() + "__" + name
        mpa_path.append(level_str)

        while not parent_taxid == '1':
            parent = child.parent
            level_rank = parent.level_rank
            parent_taxid = parent.taxid
            name = parent.name
            name = name.replace(' ','_')
            try:
                if level_rank not in main_lvls:
                    level_rank = "x"
                elif level_rank == "K":
                    level_rank = "k"
                elif level_rank == "D":
                    level_rank = "d"
                level_str = level_rank.lower() + "__" + name
                mpa_path.append(level_str)
            except ValueError:
                raise
            child = parent # needed for recursion        

        mpa_path = "|".join(map(str, mpa_path[::-1]))
        return mpa_path

    def get_taxon_path(self):

        kept_levels = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        lineage_taxid = []
        lineage_name = []
        name = self.name
        rank = self.level_rank
        name = name.replace(' ','_')
        lineage_taxid.append(self.taxid)
        lineage_name.append(name)
        child, parent = self, None
        while not rank == 'D':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            name = parent.name
            name = name.replace(' ','_')
            if rank in kept_levels:
                lineage_taxid.append(parent_taxid)
                lineage_name.append(name)
            child = parent # needed for recursion
        taxid_path = "|".join(map(str, lineage_taxid[::-1]))
        taxsn_path = "|".join(map(str, lineage_name[::-1]))
        return [taxid_path, taxsn_path]

def make_dicts(ktaxonomy_file):
    #Parse taxonomy file 
    root_node = -1
    taxid2node = {}
    with open(ktaxonomy_file, 'r') as kfile:
        for line in kfile:
            [taxid, p_tid, rank, lvl_num, name] = line.strip().split('\t|\t')
            curr_node = Tree(taxid, name, rank, lvl_num, p_tid)
            taxid2node[taxid] = curr_node
            #set parent/kids
            if taxid == "1":
                root_node = curr_node
            else:
                curr_node.parent = taxid2node[p_tid]
                taxid2node[p_tid].add_child(curr_node)
            #set parent/kids
            if taxid == "1":
                root_node = curr_node
            else:
                curr_node.parent = taxid2node[p_tid]
                taxid2node[p_tid].add_child(curr_node)            
    return taxid2node

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
    parser.add_argument('--ktaxonomy', required=True,
        help='Kraken2 database ktaxonomy file path')
    parser.add_argument('--inspect', required=True,
        dest="inspect_file", help='Kraken2 database inspect file path')
    parser.add_argument('--kmer_len', required=False,
        default=35, help='Kraken classifer kmer length [default=35]')
    parser.add_argument('--exclude', required=False,
        default=9606, nargs='+',
        help='Taxonomy ID[s] of reads to exclude (space-delimited)')
    parser.add_argument('--cluster', required=True,
        help='barcode cluster file path')
    parser.add_argument('--nsample', required=False,
        default=2000,
        help='Max number of reads to sample per taxa')
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

    logger.info('Parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='run')
    try:
        taxid2node = make_dicts(args.ktaxonomy)
        logger.info('Successfully parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='complete')
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

    logger.info('Reading kraken2 bacteria, fungi, virus classifier rank infomation from mpa report', status='run')
    # krak_mpa_report = pd.read_csv(args.krak_mpa_report_file, sep='\t', names=['mpa_taxa','reads'])
    # krak_mpa_report['taxa'] = krak_mpa_report['mpa_taxa'].apply(lambda x: re.sub(r'[a-z]__', '', x.split('|')[-1]))

    # # we only focus on  k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea
    # keywords = ["k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea"]

    # krak_mpa_report_subset = krak_mpa_report[krak_mpa_report['mpa_taxa'].str.contains('|'.join(keywords))]
    # df = pd.merge(krak_mpa_report_subset, krak_report, left_on='taxa', right_on='scientific name')

    # filter kraken_file to species and subspecies only
    desired_krak_report = krak_report.copy()[krak_report['classification_rank'].str.startswith(('S'), na=False)]
    desired_krak_report['species_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("S"), axis=1)
    desired_krak_report['superkingdom'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("D"), axis=1)
    # desired_krak_report = desired_krak_report[desired_krak_report['scientific name'].isin(df['scientific name'])] 
    desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    desired_krak_report['ncbi_taxa'] = desired_krak_report['ncbi_taxa'].astype(str)
    ## select microbiome
    desired_krak_report['is_microbiome'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].is_microbiome(), axis=1)
    desired_krak_report = desired_krak_report[desired_krak_report["is_microbiome"]==True]
    desired_taxid_list = set(desired_krak_report['ncbi_taxa'].unique())
    desired_species_taxid_list = set(desired_krak_report['species_level_taxid'].unique())
    logger.info('Finished processing kraken2 classifier result', status='complete')
    # del df

    lineage_dict = {}
    for species_tax_id in desired_species_taxid_list:
        lineage_taxid_list = taxid2node[species_tax_id].lineage_to_desired_rank("D")
        lineage_dict[species_tax_id] = lineage_taxid_list
    child_dict = {}
    for species_tax_id in desired_species_taxid_list:
        lineage_taxid_list = taxid2node[species_tax_id].lineage_to_desired_rank("S")
        child_dict[species_tax_id] = lineage_taxid_list

    descendants_dict = {}
    for species_tax_id in desired_species_taxid_list:
        descendants_taxid_list = []
        descendants_taxid_list.append(species_tax_id)
        descendants_nodes_list = taxid2node[species_tax_id].children
        while len(descendants_nodes_list) > 0:
            #For this node
            curr_n = descendants_nodes_list.pop()
            descendants_taxid_list.append(curr_n.taxid)
        descendants_dict[species_tax_id] = descendants_taxid_list

    descendants_ascendants_dict = {}
    for species_tax_id in desired_species_taxid_list:
        descendants_ascendants_taxid_list = []
        descendants_ascendants_taxid_list.append(species_tax_id)
        descendants_ascendants_taxid_list.append(taxid2node[species_tax_id].parent.taxid)
        descendants_nodes_list = taxid2node[species_tax_id].children
        while len(descendants_nodes_list) > 0:
            #For this node
            curr_n = descendants_nodes_list.pop()
            descendants_ascendants_taxid_list.append(curr_n.taxid)
        descendants_ascendants_dict[species_tax_id] = descendants_ascendants_taxid_list
    
    # Reading kraken2 classifier output information
    logger.info('Reading kraken2 classifier output information', status='run')
    logger.info('Pre-calculate and store taxid, kmer_count, and taxid index information', status='run')
    taxid_counts = {}
    # save_readids = {}
    # save_readids2 = {} 
    kraken_data = {}
    with open(args.krak_output_file, 'r') as kfile:
        for kraken_line in kfile:
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                read_type,query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                tax_id = str(re.search(r'\(taxid (\d+)\)', taxid_info).group(1))
            except:
                # in this case, something is wrong!
                print("Here is an error. Queryname: {}".format(query_name))
                continue            

            if tax_id == "-1":
                continue
            #Skip if reads are human/artificial/synthetic
            if (tax_id in desired_taxid_list):
                if tax_id not in taxid_counts:
                    taxid_counts[tax_id] = 1
                else:
                    taxid_counts[tax_id] += 1
                if taxid_counts[tax_id] >= args.nsample:
                    continue 
                if tax_id in desired_species_taxid_list:
                    species_tax_id = tax_id
                else:
                    species_tax_id = taxid2node[tax_id].taxid_to_desired_rank("S")
                kraken_data[query_name] = [species_tax_id, read_len, kmer_position]
            else:
                continue
    logger.info('Finishing reading kraken2 classifier output information', status='complete')


    # Get species level taxid
    # kraken_df['species_level_taxid'] = kraken_df.apply(lambda x: taxid_to_desired_rank(str(x['taxid']), 'species', child_parent, taxid_rank), axis=1)
    # kraken_df['species_level_taxid'] = kraken_df.apply(lambda x: taxid_to_desired_rank(str(x['taxid']), 'species', child_parent, taxid_rank)['species_taxid'], axis=1)
    num_unique_species = len(desired_krak_report['species_level_taxid'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids', status='summary')

    logger.info(f'Get the raw classified reads from bam file', status='run')
    # Init bam count
    skipped = 0
    # Init bam count
    read_count = 0
    use_count = 0
    krak_count = 0
    # # Create a defaultdict to map species_level_taxid to its set of all kmers
    # taxid_to_kmers = defaultdict(set)
    # Create a dictionary to map CB and taxid to its set of all UB and kmers
    cb_taxid_to_ub_kmers = defaultdict(lambda: {"kmers": []})  # Using a nested defaultdict
    kmer_map = defaultdict()
    species_metrics_list =[]
    species_conf_list = []
    logger.info(f'Parsing the raw classified reads from bam file', status='run')

    with pysam.AlignmentFile(args.bam_file, "rb") as krak_bamfile:
        # 遍历BAM文件和krak2_output数据
        for sread in krak_bamfile:
            read_count += 1
            # Check if the read exists in the kraken file
            if sread.query_name not in kraken_data:
                skipped += 1
                # logging.warning("Read name {} not found in kraken file".format(sread.query_name))
                continue
                
            # Get cell barcode and UMI from bam file
            try:
                sread_CB = sread.get_tag(args.barcode_tag)
            except:
                # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
                skipped += 1
                continue
            # Use the kraken data for this read
            kread = kraken_data[sread.query_name]
            krak_count += 1
            if len(sread.seq) < int(kread[1])-1:
                continue
            if kread[2].strip() == "":
                r1_conf_score  = 0
                r1_rtl_score = 0
                r1_host_score = 0
                pass
            else:
                use_count += 1
                kmer_positions_tuple = np.array([list(map(str, info.split(":"))) for info in kread[2].strip().split()])
                total_kmer_count = np.sum(kmer_positions_tuple[:, 1].astype(int))
                # Calculate selected kmer counts for specific taxids
                selected_taxa = np.concatenate((["0"], lineage_dict[kread[0]]))
                selected_mask = np.isin(kmer_positions_tuple[:, 0], selected_taxa)
                selected_kmer_count = np.sum(kmer_positions_tuple[selected_mask, 1].astype(int))
                selected_rtl_taxa = descendants_ascendants_dict[kread[0]]
                selected_rtl_mask = np.isin(kmer_positions_tuple[:, 0], selected_rtl_taxa)
                selected_rtl_kmer_count = np.sum(kmer_positions_tuple[selected_rtl_mask, 1].astype(int))
                selected_host_taxa = ["9606","9605"]
                selected_host_mask = np.isin(kmer_positions_tuple[:, 0], selected_host_taxa)
                selected_host_kmer_count = np.sum(kmer_positions_tuple[selected_host_mask, 1].astype(int))
                selected_conf_taxa = descendants_dict[kread[0]]
                selected_conf_mask = np.isin(kmer_positions_tuple[:, 0], selected_conf_taxa)
                selected_conf_kmer_count = np.sum(kmer_positions_tuple[selected_conf_mask, 1].astype(int))
                # Calculate the percentage of selected k-mer counts out of total k-mer counts
                selected_percentage = selected_kmer_count / total_kmer_count
                r1_conf_score = selected_conf_kmer_count / total_kmer_count
                r1_rtl_score = selected_rtl_kmer_count / total_kmer_count
                r1_host_score = selected_host_kmer_count / total_kmer_count

                species_conf_list.append([
                                            kread[0], r1_conf_score, r1_rtl_score, r1_host_score
                                        ])
                # If the selected percentage is less than min_frac, skip
                if selected_percentage < args.min_frac:
                    pass
                else:
                    # Init position    
                    position = 0
                    for (tax, kmer_count) in kmer_positions_tuple:
                        kmer_count =int(kmer_count)
                        xmer = sread.seq[position:position + args.kmer_len + kmer_count -1]
                        if tax not in ("0", "28384", "1","A") :
                            kmers = [xmer[i:i + args.kmer_len] for i in range(0, len(xmer) - args.kmer_len + 1)]
                            for kmer in kmers:
                                if(kmer in kmer_map):
                                    kmer_map[kmer] = "D"
                                kmer_map[kmer] = tax
                            if tax in lineage_dict[kread[0]]:
                                seq_kmer_consistency = kmer_consistency(xmer)
                                seq_entropy = calculate_entropy(xmer)
                                seq_dust_score = dust_score(xmer)
                                seq_length = len(xmer)
                                species_metrics_list.append([
                                    kread[0], seq_kmer_consistency, seq_entropy, seq_dust_score, seq_length
                                ])
                                key = (sread_CB,kread[0])
                                # cb_taxid_to_ub_kmers[key]["kmers"].extend(kmers)
                                cb_taxid_to_ub_kmers[key]["kmers"].extend(kmers)

                        position = position + kmer_count

    taxMap = dict()
    count = 0
    for xmer, taxId in kmer_map.items():
        if taxId == "dup":
            continue
        if(taxId in taxMap):
            taxMap[taxId] = taxMap[taxId] + len(xmer)
        else:
            taxMap[taxId] = len(xmer)
    taxa_nucleotides_df = pd.DataFrame.from_dict(taxMap, orient='index', columns=['nucleotides'])
    taxa_nucleotides_df.reset_index(level=0, inplace=True)
    taxa_nucleotides_df.rename(columns={'index': 'ncbi_taxa'}, inplace=True)

    ## Get the final species_seq_metrics
    all_species_conf_metrics =[]
    # Convert the species_metrics_list to a DataFrame
    conf_columns = [
        'species_level_taxid', 'mean_seq_confidence_score',
        'mean_seq_rtl_score',"mean_seq_host_score"
    ]
    species_conf_df = pd.DataFrame(species_conf_list, columns=conf_columns)
    species_conf_df['mean_seq_confidence_score'] = species_conf_df['mean_seq_confidence_score'].astype(float)
    species_conf_df['mean_seq_rtl_score'] = species_conf_df['mean_seq_rtl_score'].astype(float)
    species_conf_df['mean_seq_host_score'] = species_conf_df['mean_seq_host_score'].astype(float)

    # 分组数据并计算分位数
    grouped = species_conf_df.groupby('species_level_taxid')
    q25 = grouped['mean_seq_confidence_score'].quantile(0.25)
    q75 = grouped['mean_seq_confidence_score'].quantile(0.75)

    # 创建包含分位数的新DataFrame
    species_conf_quantile_df = pd.DataFrame({'species_level_taxid': q25.index, 'mean_seq_confidence_score_q25': q25.values, 'mean_seq_confidence_score_q75': q75.values})

    # 同样的步骤计算另一个列的分位数
    q25_rtl = grouped['mean_seq_rtl_score'].quantile(0.25)
    q75_rtl = grouped['mean_seq_rtl_score'].quantile(0.75)

    # 将另一个列的分位数合并到新DataFrame中
    species_conf_quantile_df['mean_seq_rtl_score_q25'] = q25_rtl.values
    species_conf_quantile_df['mean_seq_rtl_score_q75'] = q75_rtl.values


    q25_host = grouped['mean_seq_host_score'].quantile(0.25)
    q75_host = grouped['mean_seq_host_score'].quantile(0.75)
    species_conf_quantile_df['mean_seq_host_score_q25'] = q25_host.values
    species_conf_quantile_df['mean_seq_host_score_q75'] = q75_host.values

    ## Get the final species_seq_metrics
    all_species_seq_metrics =[]
    # Convert the species_metrics_list to a DataFrame
    metrics_columns = [
        'species_level_taxid', 'seq_kmer_consistency',
        'seq_entropy', 'seq_dust_score', 'seq_length'
    ]
    species_metrics_df = pd.DataFrame(species_metrics_list, columns=metrics_columns)

    # Group by species and calculate statistics
    species_seq_metrics = species_metrics_df.groupby('species_level_taxid').agg({
        'seq_kmer_consistency': 'mean',
        'seq_entropy': 'mean',
        'seq_dust_score': 'mean',
        'seq_length': ['max', 'mean']
    }).reset_index()

    # Rename the columns
    species_seq_metrics.columns = [
        'species_level_taxid', 'average_kmer_consistency',
        'average_seq_entropy', 'average_seq_dust_score',
        'max_seq_length', 'mean_seq_length'
    ]

    # Append the metrics data for the current ID to the respective lists
    all_species_seq_metrics.append(species_seq_metrics)

    ## Get the final species_seq_metrics
    species_seq_metrics = pd.concat(all_species_seq_metrics, ignore_index=True)

    logger.info(f'Finished parsing the raw classified reads from bam file', status='run')
    logger.info(f'Total unmapped reads: {read_count}', status='summary')
    logger.info(f'Total classified Reads with CB and UB: {use_count}', status='summary')
    logger.info(f'Skipped reads: {skipped}', status='summary')

    data = [{"CB": cb, "species_level_taxid": species_level_taxid, "kmers": kmers["kmers"]} 
            for (cb, species_level_taxid), kmers in cb_taxid_to_ub_kmers.items()]

    # Create the DataFrame from the list of dictionaries
    cb_taxid_ub_kmer_count_df = pd.DataFrame(data)
    # Del data
    del data
    del cb_taxid_to_ub_kmers

    num_unique_CB = len(cb_taxid_ub_kmer_count_df['CB'].unique())
    if num_unique_CB > 300:
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
        del cb_taxid_ub_kmer_count_df

        cluster_df = pd.read_csv(args.cluster,sep="\t",)
        cluster_df.columns = ["barcode","leiden"]
        cb_cluster_taxid_kmer_count_df = pd.merge(cb_taxid_kmer_count_df, cluster_df, left_on='CB', right_on='barcode', how='left')
        cb_cluster_taxid_kmer_count_df['leiden'] = cb_cluster_taxid_kmer_count_df['leiden'].astype(str)
        total_unique_CB = cb_cluster_taxid_kmer_count_df['CB'].nunique()
        total_unique_cluster = cb_cluster_taxid_kmer_count_df['leiden'].nunique()
        # species_prevalence = cb_cluster_taxid_kmer_count_df.groupby('species_level_taxid')['CB'].nunique().reset_index()
        # # Rename the columns for clarity
        # species_prevalence.columns = ['species_level_taxid', 'unique_CB_count']
        # species_prevalence['CB_prevalence'] = species_prevalence['unique_CB_count'] / total_unique_CB
        # Calculate the prevalence for each species_level_taxid in CBs
        species_cb_prevalence = cb_cluster_taxid_kmer_count_df.groupby('species_level_taxid')['CB'].nunique().reset_index()
        species_cb_prevalence.columns = ['species_level_taxid', 'unique_CB_count']
        species_cb_prevalence['CB_prevalence'] = species_cb_prevalence['unique_CB_count'] / total_unique_CB

        # First, calculate the count of each species_level_taxid within each cluster
        species_cluster_counts = cb_cluster_taxid_kmer_count_df.groupby(['species_level_taxid', 'leiden']).size().reset_index(name='count')

        # Filter out the rows where a species_level_taxid appears in a cluster less than 3 times
        species_cluster_filtered = species_cluster_counts[species_cluster_counts['count'] >= 3]

        # Calculate the prevalence for each species_level_taxid in clusters, based on the filtered data
        species_cluster_prevalence = species_cluster_filtered.groupby('species_level_taxid')['leiden'].nunique().reset_index()
        species_cluster_prevalence.columns = ['species_level_taxid', 'unique_cluster_count']
        species_cluster_prevalence['cluster_prevalence'] = species_cluster_prevalence['unique_cluster_count'] / total_unique_cluster

        # Merge the two prevalence dataframes on species_level_taxid
        species_prevalence_combined = pd.merge(species_cb_prevalence, species_cluster_prevalence, on='species_level_taxid')
        species_prevalence_combined['species_level_taxid'] = species_prevalence_combined['species_level_taxid'].astype(str)

        del cb_cluster_taxid_kmer_count_df
        del cb_taxid_ub_global_unique_count_df
        
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
    # final_desired_krak_report = final_desired_krak_report.groupby('species_level_taxid').apply(calculate_g_score)
    ## Reset index
    final_desired_krak_report.reset_index(drop=True, inplace=True)
    final_desired_krak_report = final_desired_krak_report.merge(species_seq_metrics,left_on='species_level_taxid', right_on='species_level_taxid')
    # final_desired_krak_report.drop('taxid', axis=1, inplace=True)

    if num_unique_CB > 300:
        final_desired_krak_report = final_desired_krak_report.merge(cb_taxid_kmer_corr_df,left_on='species_level_taxid', right_on='species_level_taxid')
        final_desired_krak_report = final_desired_krak_report.merge(species_prevalence_combined,left_on='species_level_taxid', right_on='species_level_taxid')


    final_desired_krak_report = final_desired_krak_report.merge(species_conf_quantile_df,left_on='species_level_taxid', right_on='species_level_taxid')
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    taxa_nucleotides_df['ncbi_taxa'] = taxa_nucleotides_df['ncbi_taxa'].astype(str)
    final_desired_krak_report = final_desired_krak_report.merge(taxa_nucleotides_df, left_on='ncbi_taxa',right_on='ncbi_taxa', how='left')
    # final_desired_krak_report.drop('species_level_taxid', axis=1, inplace=True)


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
    if ntests > 5:
        # filter_desired_krak_report = final_desired_krak_report.copy()[
        #     (final_desired_krak_report['average_seq_entropy'] > args.min_entropy) &
        #     (final_desired_krak_report['max_minimizers'] > 5) &
        #     (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
        #     (
        #         (
        #             ((final_desired_krak_report['superkingdom'] == '2') & (final_desired_krak_report['max_cov']*50 >=  bacteria_mean_cov)) |
        #             ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov']*50 >=  archaea_mean_cov)) |
        #             ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >  0)) |
        #             ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0.01)) 
        #         )
        #     ) &
        #     (
        #         (
        #             (final_desired_krak_report['corr_kmer_uniq_counts'] > 0.5) &
        #         ((final_desired_krak_report['p_value_kmer_uniq_counts'] < float(0.05)) | (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()))
        #         ) |
        #         (
        #             (final_desired_krak_report['corr_kmer_uniq_counts'].isna()) &
        #             (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()) &
        #             (final_desired_krak_report['average_seq_entropy'] > 1.5*args.min_entropy) &
        #             (final_desired_krak_report['max_seq_length'] > float(38))
        #         )
        #     )
        # ]
        filter_desired_krak_report = final_desired_krak_report.copy()[
            (
                (
                (final_desired_krak_report['average_seq_entropy'] > 1.2) &
                (final_desired_krak_report['max_minimizers'] > 5) &
                (final_desired_krak_report['cluster_prevalence'] < 0.8) &
                (final_desired_krak_report['CB_prevalence'] < 0.15) &
                (final_desired_krak_report['unique_CB_count'] >= 5) &
                # (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
                # (final_desired_krak_report['mean_seq_rtl_score_q75'] > 0.3) &
                # (
                #     (
                #         ((final_desired_krak_report['superkingdom'] == '2') & (final_desired_krak_report['max_cov'] >= 0.00001)) |
                #         ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov'] >=  0.00001)) |
                #         ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >  0)) |
                #         ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0.001)) 
                #     )
                # ) &
                (
                    # (
                    #     (final_desired_krak_report['corr_kmer_uniq_counts'] > 0.5) &
                    #     ((final_desired_krak_report['p_value_kmer_uniq_counts'] < float(0.05)) | (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()))
                    # )
                    (
                        (final_desired_krak_report['corr_kmer_uniq_counts'] > 0) &
                        (final_desired_krak_report['p_value_kmer_uniq_counts'] < float(0.05)) &
                        (final_desired_krak_report['corr_kmer_global_uniq_counts'] > 0) 
                        )
                    )
                    # ) |
                    # (
                    #     (final_desired_krak_report['corr_kmer_uniq_counts'].isna()) &
                    #     (final_desired_krak_report['p_value_kmer_uniq_counts'].isna()) &
                    #     (final_desired_krak_report['nucleotides']> 100)
                    # )
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
                    ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] >  0.01)) 
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
    num_unique_species = len(filter_desired_krak_report['species_level_taxid'].unique())
    logger.info(f'After filtering, found {num_unique_species} unique species level taxids', status='summary')
    # Save data
    logger.info(f'Saving the result', status='run')
    filter_desired_krak_report.to_csv(args.qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')


if __name__ == "__main__":
    main()