import pysam
import logging
import os
import gzip
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


# Calculate k-mer consistency
def kmer_consistency(sequence, k=6):
    """
    Calculate the k-mer consistency of a DNA or RNA sequence.

    Parameters:
    - sequence (str): Input DNA or RNA sequence.
    - k (int): Length of k-mers, default is 6.

    Returns:
    - float: k-mer consistency score, indicating the proportion of different k-mers.
    """
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    kmer_counts = Counter(kmers)
    return len(kmer_counts) / len(kmers)

# Calculate DUST score
def dust_score(sequence):
    """
    Calculate the DUST (proposed by NCBI) score of a DNA or RNA sequence.

    Parameters:
    - sequence (str): Input DNA or RNA sequence.

    Returns:
    - float: DUST score, representing the proportion of different characters in the sequence.
    """
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

# Calculate entropy
def calculate_entropy(sequence):
    """
    Calculate the information entropy of a DNA or RNA sequence.

    Parameters:
    - sequence (str): Input DNA or RNA sequence.

    Returns:
    - float: Entropy of the sequence, indicating the uncertainty of base distribution.
    """
    nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    sequence_length = len(sequence)
    
    for nucleotide in sequence:
        if nucleotide in nucleotide_counts:
            nucleotide_counts[nucleotide] += 1
    
    nucleotide_probabilities = [count / sequence_length for count in nucleotide_counts.values()]
    nucleotide_probabilities = [p for p in nucleotide_probabilities if p > 0]  # Calculate only for non-zero probability nucleotides
    entropy = -np.sum(nucleotide_probabilities * np.log2(nucleotide_probabilities))
    
    return entropy

# Tree Class
# A class representing a node in a taxonomy tree used for constructing hierarchical taxonomic structures.
# This class includes taxonomy levels and genomes identified in the Kraken report.
class Tree(object):
    'Tree node.'
    
    def __init__(self, taxid, name, level_rank, level_num, p_taxid, parent=None, children=None):
        """
        Initializes a Tree node with taxonomic information.

        Parameters:
        - taxid (int): Taxonomic identifier.
        - name (str): Name of the taxonomic entity.
        - level_rank (str): Taxonomic level rank (e.g., 'D' for Domain, 'P' for Phylum).
        - level_num (int): Numeric representation of the taxonomic level.
        - p_taxid (int): Parent taxonomic identifier.
        - parent (Tree): Parent node in the tree.
        - children (List[Tree]): List of child nodes.
        """
        self.taxid = taxid
        self.name = name
        self.level_rank = level_rank
        self.level_num = int(level_num)
        self.p_taxid = p_taxid
        self.all_reads = 0
        self.lvl_reads = 0

        # Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        """
        Adds a child node to the current node.

        Parameters:
        - node (Tree): Child node to be added.
        """
        assert isinstance(node, Tree)
        self.children.append(node)

    def taxid_to_desired_rank(self, desired_rank):
        """
        Retrieves the taxonomic identifier at the desired rank.

        Parameters:
        - desired_rank (str): Desired taxonomic rank.

        Returns:
        - int or str: Taxonomic identifier at the desired rank or an error message.
        """
        # Check if the current node's level_rank matches the desired_rank
        if self.level_rank == desired_rank:
            return self.taxid

        child, parent, parent_taxid = self, None, None
        while not parent_taxid == '1':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            if rank == desired_rank:
                return parent.taxid
            child = parent  # needed for recursion

        # If no parent node is found or the desired_rank is not reached, return an error
        return 'error - taxid above desired rank, or not annotated at desired rank'

    def lineage_to_desired_rank(self, desired_parent_rank):
        """
        Retrieves the taxonomic lineage up to the desired parent rank.

        Parameters:
        - desired_parent_rank (str): Desired parent taxonomic rank.

        Returns:
        - List[int]: List of taxonomic identifiers in the lineage up to the desired parent rank.
        """
        lineage = []
        lineage.append(self.taxid)

        # Check if the current node's level_num is at the top level (1)
        if self.level_num == "1":
            return lineage

        if self.level_rank in {"S", "G"}:
            for child in self.children:
                lineage.extend(child.get_all_descendants())

        child, parent, parent_taxid = self, None, None
        while not parent_taxid == '1':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            lineage.append(parent_taxid)
            if rank == desired_parent_rank:
                return lineage
            child = parent  # needed for recursion
        return lineage

    def get_main_lvl_taxid(self):
        """
        Retrieves the taxonomic identifier at the main taxonomic level.

        Returns:
        - int: Taxonomic identifier at the main taxonomic level.
        """
        main_lvls = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        level_rank = self.level_rank
        child, parent, parent_taxid = self, None, None

        while level_rank not in main_lvls:
            parent = child.parent
            level_rank = parent.level_rank
            child = parent  # needed for recursion

        main_lvl_taxid = child.taxid
        return main_lvl_taxid
    def get_all_descendants(self):
        """
        Get the taxids of all descendants in the subtree rooted at the current node.

        Returns:
        - list: List of taxids for all descendants in the subtree.
        """
        descendants_taxids = []

        descendants_taxids.append(self.taxid)

        for child in self.children:
            descendants_taxids.extend(child.get_all_descendants())

        return descendants_taxids
    def get_mpa_path(self):
        """
        Retrieves the taxonomic path formatted for the Metagenomics Pathway Analysis (MPA) tool.

        Returns:
        - str: Formatted taxonomic path for MPA.
        """
        mpa_path = []
        main_lvls = ['D', 'P', 'C', 'O', 'F', 'G', 'S']

        # Create level name
        level_rank = self.level_rank
        name = self.name
        name = name.replace(' ', '_')

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
            name = name.replace(' ', '_')

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

            child = parent  # needed for recursion
        # Reverse the MPA path list and join its components with "|".
        mpa_path = "|".join(map(str, mpa_path[::-1])) 
        return mpa_path

    def is_microbiome(self):
        """
        Checks if the taxonomic node represents a microbiome entity.

        Returns:
        - bool: True if the node represents a microbiome, False otherwise.
        """
        is_microbiome = False
        main_lvls = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        lineage_name = []

        # Create level name
        level_rank = self.level_rank
        name = self.name
        name = name.replace(' ', '_')
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
            name = name.replace(' ', '_')
            lineage_name.append(name)
            child = parent  # needed for recursion

        if 'Fungi' in lineage_name or 'Bacteria' in lineage_name or 'Viruses' in lineage_name:
            is_microbiome = True
        return is_microbiome

    def get_taxon_path(self):
        """
        Retrieves the taxonomic path including taxonomic identifiers and names.

        Returns:
        - List[str]: List containing taxonomic path as taxonomic identifiers and names.
        """
        kept_levels = ['D', 'P', 'C', 'O', 'F', 'G', 'S']
        lineage_taxid = []
        lineage_name = []
        name = self.name
        rank = self.level_rank
        name = name.replace(' ', '_')
        lineage_taxid.append(self.taxid)
        lineage_name.append(name)

        child, parent = self, None
        while not rank == 'D':
            parent = child.parent
            rank = parent.level_rank
            parent_taxid = parent.taxid
            name = parent.name
            name = name.replace(' ', '_')
            if rank in kept_levels:
                lineage_taxid.append(parent_taxid)
                lineage_name.append(name)
            child = parent  # needed for recursion

        taxid_path = "|".join(map(str, lineage_taxid[::-1]))
        taxsn_path = "|".join(map(str, lineage_name[::-1]))
        return [taxid_path, taxsn_path]

def make_dicts(ktaxonomy_file):
    """
    Parse a Kraken taxonomy file and create a dictionary of Tree nodes.

    Parameters:
    - ktaxonomy_file (str): Path to the Kraken taxonomy file.

    Returns:
    - dict: Dictionary mapping taxonomic identifiers to Tree nodes.
    """
    root_node = -1  # Initialize the root node identifier.
    taxid2node = {}  # Dictionary to store Tree nodes mapped to their taxonomic identifiers.

    with open(ktaxonomy_file, 'r') as kfile:
        for line in kfile:
            # Parse the tab-separated values from each line of the Kraken taxonomy file.
            [taxid, p_tid, rank, lvl_num, name] = line.strip().split('\t|\t')
            
            # Create a Tree node for the current taxonomic entry.
            curr_node = Tree(taxid, name, rank, lvl_num, p_tid)
            
            # Add the current node to the taxid2node dictionary.
            taxid2node[taxid] = curr_node
            
            # Set parent and children relationships for the current node.
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
    parser.add_argument('--cluster', required=True,
        help='barcode cluster file path')
    parser.add_argument('--inspect', required=True,
        dest="inspect_file", help='Kraken2 database inspect file path')
    parser.add_argument('--kmer_len', required=False,
        default=35, help='Kraken classifer kmer length [default=35]')
    parser.add_argument('--exclude', required=False,
        default=9606, nargs='+',
        help='Taxonomy ID[s] of reads to exclude (space-delimited)')
    parser.add_argument('--min_frac', required=False,
        default=0.5, type=float, help='minimum fraction of kmers directly assigned to taxid to use read [default=0.5]')
    parser.add_argument('--nsample', required=False,
        default=2500,
        help='Max number of reads to sample per taxa')
    parser.add_argument('--min_entropy', required=False,
        default=1.2, type=float, help='minimum entropy of sequences cutoff [default=1.2]')
    parser.add_argument('--min_dust', required=False,
        default=0.05, type=float, help='minimum dust score of sequences cutoff [default=1.2]')
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
    except (ValueError, KeyError) as e:
        logger.error(f"An error occurred while processing the Kraken taxonomy file: {e}")
        sys.exit()

    logger.info('Reading kraken2 classifier result infomation from report', status='run')
    krak_report = pd.read_csv(args.krak_report_file, sep="\t", names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
    # remove space
    krak_report['scientific name'] = krak_report['scientific name'].str.strip() 
    # replace space
    krak_report['scientific name'] = krak_report['scientific name'].str.replace(r' ', '_')
    total_reads = krak_report['fragments'].iloc[0] + krak_report['fragments'].iloc[1]
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

    # filter kraken_file to species and genus only
    desired_krak_report = krak_report.copy()[krak_report['classification_rank'].str.startswith((('G', 'S')), na=False)]
    desired_krak_report['species_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("S"), axis=1)
    desired_krak_report['main_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].get_main_lvl_taxid(), axis=1)
    desired_krak_report['genus_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("G"), axis=1)
    desired_krak_report['superkingdom'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("D"), axis=1)
    desired_krak_report['is_microbiome'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].is_microbiome(), axis=1)

    ## select microbiome
    desired_krak_report = desired_krak_report[desired_krak_report["is_microbiome"]==True]
    # Transform data type
    desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    # desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    # desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    desired_krak_report['ncbi_taxa'] = desired_krak_report['ncbi_taxa'].astype(str)
    
    # desired_krak_report
    desired_taxid_list = set(desired_krak_report['ncbi_taxa'].unique())
    desired_main_taxid_list = set(desired_krak_report['main_level_taxid'].unique())
    logger.info('Finished processing kraken2 classifier result', status='complete')
    # del df

    lineage_dict = {}
    for main_tax_id in desired_main_taxid_list:
        try:
            lineage_taxid_list = taxid2node[main_tax_id].lineage_to_desired_rank("D")
            lineage_dict[main_tax_id] = lineage_taxid_list
        except (ValueError, KeyError) as e:
            print("Error occur:", e)
    descendants_dict = {}
    for main_tax_id in desired_main_taxid_list:
        try:
            descendants_taxid_list = taxid2node[main_tax_id].get_all_descendants()
            descendants_dict[main_tax_id] = descendants_taxid_list
        except (ValueError, KeyError) as e:
            print("Error occur:", e)


    conf_dict = {}
    for tax_id in desired_main_taxid_list:
        if tax_id == 'error - taxid above desired rank, or not annotated at desired rank':
            continue
        descendants_taxid_list = []
        descendants_taxid_list.append(tax_id)
        descendants_nodes_list = taxid2node[tax_id].children
        while len(descendants_nodes_list) > 0:
            #For this node
            curr_n = descendants_nodes_list.pop()
            descendants_taxid_list.append(curr_n.taxid)
        conf_dict[tax_id] = descendants_taxid_list

    try:
        taxid2node = make_dicts(args.ktaxonomy)
    except (ValueError, KeyError) as e:
        logger.error(f"An error occurred while processing the Kraken taxonomy file: {e}")
        sys.exit()

    rtl_dict = {}
    for species_tax_id in desired_main_taxid_list:
        descendants_ascendants_taxid_list = []
        descendants_ascendants_taxid_list.append(species_tax_id)
        descendants_ascendants_taxid_list.append(taxid2node[species_tax_id].parent.taxid)
        descendants_nodes_list = taxid2node[species_tax_id].children
        while len(descendants_nodes_list) > 0:
            #For this node
            curr_n = descendants_nodes_list.pop()
            descendants_ascendants_taxid_list.append(curr_n.taxid)
        rtl_dict[species_tax_id] = descendants_ascendants_taxid_list

    # Check file is pair end
    is_paired = False
    read_count = 0
    with pysam.AlignmentFile(args.bam_file, 'rb') as krak_bamfile:
            for kread in krak_bamfile:
                read_count += 1
                if kread.is_paired:
                    is_paired = True
                    break
                else:
                    continue
                if read_count >= 100:
                    break

    # Reading kraken2 classifier output information
    logger.info('Reading kraken2 classifier output information', status='run')
    taxid_counts = {}
    kraken_data = {}
    if is_paired:
        with open(args.krak_output_file, 'r') as kfile:
            for kraken_line in kfile:
                try:
                    # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                    # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                    # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                    read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                    tax_id = str(re.search(r'\(taxid (\d+)\)', taxid_info).group(1))
                except (ValueError, KeyError) as e:
                    # in this case, something is wrong!
                    logger.error(f"An error occurred while processing the Kraken output file: {e}")
                    logger.error(f"Here is an error. Queryname: {query_name}")
                    continue
                if tax_id == "-1":
                    continue
                #Skip if reads are human/artificial/synthetic
                if (tax_id in desired_taxid_list):
                    r1_len, r2_len = read_len.split('|')
                    r1_kmer_position, r2_kmer_position  = kmer_position.split(' |:| ')
                    if tax_id not in taxid_counts:
                        taxid_counts[tax_id] = 1
                    else:
                        taxid_counts[tax_id] += 1
                    if taxid_counts[tax_id] >= args.nsample:
                        continue 
                    if tax_id in desired_main_taxid_list:
                        main_lvl_tax_id = tax_id
                    else:
                        main_lvl_tax_id = taxid2node[tax_id].get_main_lvl_taxid()
                    kraken_data[query_name] = [main_lvl_tax_id, r1_len, r2_len,r1_kmer_position, r2_kmer_position]
                else:
                    continue
    else:
        with open(args.krak_output_file, 'r') as kfile:
            for kraken_line in kfile:
                try:
                    # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                    # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                    # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                    read_type,query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                    tax_id = str(re.search(r'\(taxid (\d+)\)', taxid_info).group(1))
                except (ValueError, KeyError) as e:
                    # in this case, something is wrong!
                    logger.error(f"An error occurred while processing the Kraken taxonomy file: {e}")
                    logger.error(f"Here is an error. Queryname: {query_name}")
                    continue

                if tax_id == "-1":
                    continue
                #Skip if reads are human/artificial/synthetic
                if (tax_id in desired_taxid_list):
                    if tax_id not in taxid_counts:
                        taxid_counts[tax_id] = 1
                    else:
                        taxid_counts[tax_id] += 1
                    if taxid_counts[tax_id] >= 1500:
                        continue 
                    if tax_id in desired_main_taxid_list:
                        main_lvl_tax_id = tax_id
                    else:
                        main_lvl_tax_id = taxid2node[tax_id].get_main_lvl_taxid()
                    kraken_data[query_name] = [main_lvl_tax_id, read_len, kmer_position]
                else:
                    continue
    logger.info('Finishing reading kraken2 classifier output information', status='complete')

    # Get species level taxid
    num_unique_species = len(desired_krak_report['species_level_taxid'].unique())
    num_unique_genus = len(desired_krak_report['genus_level_taxid'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids and {num_unique_genus} unique genus level taxids', status='summary')

    logger.info(f'Get the raw classified reads from bam file', status='run')

    # Init bam count
    skipped = 0
    read_count = 0
    use_count = 0
    cb_count = 0
    krak_count = 0
    # Create a dictionary to map CB and taxid to its set of all UB and kmers
    cb_taxid_to_ub_kmers = defaultdict(lambda: {"kmers": []})  # Using a nested defaultdict
    kmer_map = defaultdict()
    species_metrics_list =[]
    species_conf_list = []
    logger.info(f'Parsing the raw classified reads from bam file', status='run')
    if is_paired:
        with pysam.AlignmentFile(args.bam_file, "rb") as krak_bamfile:
            # iter bam and krak output
            for sread_r1,sread_r2 in zip(krak_bamfile, krak_bamfile):
                read_count += 1
                # if read_count % 100000 == 0:
                #     sys.stdout.write('\r\t%0.3f million reads processed' % float(read_count/1000000.))
                #     sys.stdout.flush()
                # Get cell barcode and UMI from bam file
                try:
                    sread_CB = sread_r1.get_tag(args.barcode_tag)
                except Exception as e:
                    # Some reads don't have a cell barcode or transcript barcode; they can be skipped.
                    skipped += 1
                    # Log the error
                    logger.error(f"Error extracting barcode from read {sread_r1.query_name}: {e}")
                    continue
                cb_count += 1
                if sread_r1.query_name not in kraken_data:
                    skipped += 1
                    continue
                # some read may be reverse
                if not sread_r1.is_read1:
                    sread_r1, sread_r2 = sread_r2, sread_r1

                # Use the kraken data for this read
                kread = kraken_data[sread_r1.query_name]

                krak_count += 1
                if len(sread_r2.seq) < int(kread[2])-1:
                    continue
                if kread[3] == "":
                    # Initialize scores
                    r1_conf_score  = 0
                    r1_rtl_score = 0
                    r1_host_score = 0
                    pass
                else:
                    r1_kmer_positions_tuple = np.array([list(map(str, info.split(":"))) for info in kread[3].strip().split()])
                    total_kmer_count = np.sum(r1_kmer_positions_tuple[:, 1].astype(int))

                    # Calculate selected kmer counts for specific taxids
                    selected_taxa = np.concatenate((["0"], lineage_dict[kread[0]]))
                    selected_mask = np.isin(r1_kmer_positions_tuple[:, 0], selected_taxa)
                    selected_kmer_count = np.sum(r1_kmer_positions_tuple[selected_mask, 1].astype(int))

                    selected_rtl_taxa = rtl_dict[kread[0]]
                    selected_rtl_mask = np.isin(r1_kmer_positions_tuple[:, 0], selected_rtl_taxa)
                    selected_rtl_kmer_count = np.sum(r1_kmer_positions_tuple[selected_rtl_mask, 1].astype(int))

                    selected_host_taxa = ["9606","9605"]
                    selected_host_mask = np.isin(r1_kmer_positions_tuple[:, 0], selected_host_taxa)
                    selected_host_kmer_count = np.sum(r1_kmer_positions_tuple[selected_host_mask, 1].astype(int))

                    selected_conf_taxa = conf_dict[kread[0]]
                    selected_conf_mask = np.isin(r1_kmer_positions_tuple[:, 0], selected_conf_taxa)
                    selected_conf_kmer_count = np.sum(r1_kmer_positions_tuple[selected_conf_mask, 1].astype(int))

                    # Calculate the percentage of selected k-mer counts out of total k-mer counts
                    selected_percentage = selected_kmer_count / total_kmer_count
                    r1_conf_score = selected_conf_kmer_count / total_kmer_count
                    r1_rtl_score = selected_rtl_kmer_count / total_kmer_count
                    r1_host_score = selected_host_kmer_count / total_kmer_count

                    # If the selected percentage is less than min_frac, skip
                    if selected_percentage < args.min_frac:
                        pass
                    else:
                        use_count += 1
                        # Init position    
                        position = 0
                        for (tax, kmer_count) in r1_kmer_positions_tuple:
                            kmer_count =int(kmer_count)
                            xmer = sread_r1.seq[position:position + args.kmer_len + kmer_count -1]
                            if tax not in ("0", "28384", "1","A") :
                                kmers = [xmer[i:i + args.kmer_len] for i in range(0, len(xmer) - args.kmer_len + 1)]

                                for kmer in kmers:
                                    if(kmer in kmer_map):                     
                                        kmer_map[kmer] = "D"
                                    kmer_map[kmer] = tax
                            if tax in descendants_dict[kread[0]]:
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
                if kread[4].strip() == "":
                    r2_conf_score  = 0
                    r2_rtl_score = 0
                    r2_host_score = 0 
                    species_conf_list.append([
                                                kread[0], r1_conf_score, r1_rtl_score, r1_host_score, r2_conf_score, r2_rtl_score, r2_host_score
                                            ])
                    pass
                else:
                    r2_kmer_positions_tuple = np.array([list(map(str, info.split(":"))) for info in kread[4].strip().split()])
                    total_kmer_count = np.sum(r2_kmer_positions_tuple[:, 1].astype(int))

                    # Calculate selected kmer counts for specific taxids
                    selected_taxa = np.concatenate((["0"], lineage_dict[kread[0]]))
                    selected_mask = np.isin(r2_kmer_positions_tuple[:, 0], selected_taxa)
                    selected_kmer_count = np.sum(r2_kmer_positions_tuple[selected_mask, 1].astype(int))
                    
                    selected_rtl_taxa = rtl_dict[kread[0]]
                    selected_rtl_mask = np.isin(r2_kmer_positions_tuple[:, 0], selected_rtl_taxa)
                    selected_rtl_kmer_count = np.sum(r2_kmer_positions_tuple[selected_rtl_mask, 1].astype(int))

                    selected_conf_taxa = conf_dict[kread[0]]
                    selected_conf_mask = np.isin(r2_kmer_positions_tuple[:, 0], selected_conf_taxa)
                    selected_conf_kmer_count = np.sum(r2_kmer_positions_tuple[selected_conf_mask, 1].astype(int))

                    selected_host_taxa = ["9606","9605"]
                    selected_host_mask = np.isin(r2_kmer_positions_tuple[:, 0], selected_host_taxa)
                    selected_host_kmer_count = np.sum(r2_kmer_positions_tuple[selected_host_mask, 1].astype(int))

                    # Calculate the percentage of selected k-mer counts out of total k-mer counts
                    selected_percentage = selected_kmer_count / total_kmer_count
                    r2_conf_score = selected_conf_kmer_count / total_kmer_count
                    r2_rtl_score = selected_rtl_kmer_count / total_kmer_count
                    r2_host_score = selected_host_kmer_count / total_kmer_count
                    species_conf_list.append([
                                                kread[0], r1_conf_score, r1_rtl_score, r1_host_score, r2_conf_score, r2_rtl_score, r2_host_score
                                            ])
                    # Calculate the percentage of selected k-mer counts out of total k-mer counts
                    selected_percentage = selected_kmer_count / total_kmer_count

                    # If the selected percentage is less than min_frac, skip
                    if selected_percentage < args.min_frac:
                        pass
                    else:
                        # Init position    
                        position = 0
                        for (tax, kmer_count) in r2_kmer_positions_tuple:
                            kmer_count =int(kmer_count)
                            xmer = sread_r2.seq[position:position + args.kmer_len + kmer_count -1]

                            if tax not in ("0", "28384", "1","A") :
                                kmers = [xmer[i:i + args.kmer_len] for i in range(0, len(xmer) - args.kmer_len + 1)]
                                for kmer in kmers:
                                    if(kmer in kmer_map):
                                        kmer_map[kmer] = "D"
                                    kmer_map[kmer] = tax
                                if tax in descendants_dict[kread[0]]:
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

    else:
        with pysam.AlignmentFile(args.bam_file, "rb") as krak_bamfile:
            # Iterate over reads in the BAM file and corresponding krak2_output data
            for sread in krak_bamfile:
                read_count += 1
                # Try to get cell barcode and UMI from BAM file
                try:
                    sread_CB = sread.get_tag(args.barcode_tag)
                except Exception as e:
                    # Some reads don't have a cell barcode or transcript barcode; they can be skipped.
                    skipped += 1
                    # Log the error
                    logger.error(f"Error extracting barcode from read {sread.query_name}: {e}")
                    continue
                cb_count += 1
                # Check if the read exists in the kraken file
                if sread.query_name not in kraken_data:
                    skipped += 1

                    continue

                
                # Use the kraken data for this read
                kread = kraken_data[sread.query_name]
                krak_count += 1

                # Continue processing only if the read sequence length is sufficient
                if len(sread.seq) < int(kread[1]) - 1:
                    continue

                # Initialize scores
                r1_conf_score = 0
                r1_rtl_score = 0
                r1_host_score = 0

                # Process kraken data for scoring
                if kread[2].strip() != "":
                    use_count += 1
                    kmer_positions_tuple = np.array([list(map(str, info.split(":"))) for info in kread[2].strip().split()])
                    total_kmer_count = np.sum(kmer_positions_tuple[:, 1].astype(int))

                    # Calculate selected kmer counts for specific taxids
                    selected_taxa = np.concatenate((["0"], lineage_dict[kread[0]]))
                    selected_mask = np.isin(kmer_positions_tuple[:, 0], selected_taxa)
                    selected_kmer_count = np.sum(kmer_positions_tuple[selected_mask, 1].astype(int))

                    # Calculate selected rtl, host, and confidence scores
                    selected_rtl_taxa = rtl_dict[kread[0]]
                    selected_rtl_mask = np.isin(kmer_positions_tuple[:, 0], selected_rtl_taxa)
                    selected_rtl_kmer_count = np.sum(kmer_positions_tuple[selected_rtl_mask, 1].astype(int))

                    selected_host_taxa = ["9606", "9605"]
                    selected_host_mask = np.isin(kmer_positions_tuple[:, 0], selected_host_taxa)
                    selected_host_kmer_count = np.sum(kmer_positions_tuple[selected_host_mask, 1].astype(int))

                    selected_conf_taxa = conf_dict[kread[0]]
                    selected_conf_mask = np.isin(kmer_positions_tuple[:, 0], selected_conf_taxa)
                    selected_conf_kmer_count = np.sum(kmer_positions_tuple[selected_conf_mask, 1].astype(int))

                    # Calculate scores
                    r1_conf_score = selected_conf_kmer_count / total_kmer_count
                    r1_rtl_score = selected_rtl_kmer_count / total_kmer_count
                    r1_host_score = selected_host_kmer_count / total_kmer_count

                    # Log species confidence score information
                    species_conf_list.append([
                        kread[0], r1_conf_score, r1_rtl_score, r1_host_score
                    ])

                    # If the selected percentage is less than min_frac, skip
                    if selected_kmer_count / total_kmer_count < args.min_frac:
                        pass
                    else:
                        # Initialize position
                        position = 0

                        # Process kmer positions in the sequence
                        for (tax, kmer_count) in kmer_positions_tuple:
                            kmer_count = int(kmer_count)
                            xmer = sread.seq[position:position + args.kmer_len + kmer_count - 1]

                            # Update kmer_map based on taxon
                            if tax not in ("0", "28384", "1", "A"):
                                kmers = [xmer[i:i + args.kmer_len] for i in range(0, len(xmer) - args.kmer_len + 1)]

                                for kmer in kmers:
                                    if kmer in kmer_map:
                                        kmer_map[kmer] = "D"
                                    kmer_map[kmer] = tax

                                # Log species metrics
                                if tax in descendants_dict[kread[0]]:
                                    seq_kmer_consistency = kmer_consistency(xmer)
                                    seq_entropy = calculate_entropy(xmer)
                                    seq_dust_score = dust_score(xmer)
                                    seq_length = len(xmer)

                                    species_metrics_list.append([
                                        kread[0], seq_kmer_consistency, seq_entropy, seq_dust_score, seq_length
                                    ])

                                    # Update cb_taxid_to_ub_kmers
                                    key = (sread_CB, kread[0])
                                    cb_taxid_to_ub_kmers[key]["kmers"].extend(kmers)

                            position = position + kmer_count

    if krak_count == 0:
        logger.error(f"No reads were classified by Kraken2. Please check the input files and try again.")
        sys.exit()
    
    ## Get the final species_seq_metrics
    all_species_seq_metrics =[]
    # Convert the species_metrics_list to a DataFrame
    metrics_columns = [
        'main_level_taxid', 'seq_kmer_consistency',
        'seq_entropy', 'seq_dust_score', 'seq_length'
    ]
    species_metrics_df = pd.DataFrame(species_metrics_list, columns=metrics_columns)

    # Group by species and calculate statistics
    species_seq_metrics = species_metrics_df.groupby('main_level_taxid').agg({
        'seq_kmer_consistency': 'mean',
        'seq_entropy': 'mean',
        'seq_dust_score': 'mean',
        'seq_length': ['max', 'mean']
    }).reset_index()

    # Rename the columns
    species_seq_metrics.columns = [
        'main_level_taxid', 'average_kmer_consistency',
        'average_seq_entropy', 'average_seq_dust_score',
        'max_seq_length', 'mean_seq_length'
    ]

    # Append the metrics data for the current ID to the respective lists
    all_species_seq_metrics.append(species_seq_metrics)

    ## Get the final species_seq_metrics
    species_seq_metrics = pd.concat(all_species_seq_metrics, ignore_index=True)
    species_seq_metrics['classification_rank'] = species_seq_metrics.apply(lambda x: taxid2node[str(x['main_level_taxid'])].level_rank, axis=1)
    species_seq_metrics['genus_level_taxid'] = species_seq_metrics.apply(lambda x: taxid2node[str(x['main_level_taxid'])].taxid_to_desired_rank("G"), axis=1)
    # Calculate average_kmer_consistency based on classification_rank
    species_seq_metrics['average_kmer_consistency'] = np.where(
        species_seq_metrics['classification_rank'].str.startswith('S'),
        species_seq_metrics.groupby('main_level_taxid')['average_kmer_consistency'].transform('mean'),
        species_seq_metrics.groupby('genus_level_taxid')['average_kmer_consistency'].transform('max')
    )

    # Calculate average_seq_entropy based on classification_rank
    species_seq_metrics['average_seq_entropy'] = np.where(
        species_seq_metrics['classification_rank'].str.startswith('S'),
        species_seq_metrics.groupby('main_level_taxid')['average_seq_entropy'].transform('mean'),
        species_seq_metrics.groupby('genus_level_taxid')['average_seq_entropy'].transform('max')
    )

    # Calculate average_seq_dust_score based on classification_rank
    species_seq_metrics['average_seq_dust_score'] = np.where(
        species_seq_metrics['classification_rank'].str.startswith('S'),
        species_seq_metrics.groupby('main_level_taxid')['average_seq_dust_score'].transform('mean'),
        species_seq_metrics.groupby('genus_level_taxid')['average_seq_dust_score'].transform('max')
    )

    # Calculate max_seq_length based on classification_rank
    species_seq_metrics['max_seq_length'] = np.where(
        species_seq_metrics['classification_rank'].str.startswith('S'),
        species_seq_metrics.groupby('main_level_taxid')['max_seq_length'].transform('max'),
        species_seq_metrics.groupby('genus_level_taxid')['max_seq_length'].transform('max')
    )

    # Calculate mean_seq_length based on classification_rank
    species_seq_metrics['mean_seq_length'] = np.where(
        species_seq_metrics['classification_rank'].str.startswith('S'),
        species_seq_metrics.groupby('main_level_taxid')['mean_seq_length'].transform('mean'),
        species_seq_metrics.groupby('genus_level_taxid')['mean_seq_length'].transform('max')
    )


    logger.info(f'Finished parsing the raw classified reads from bam file', status='run')
    logger.info(f'Total unmapped reads: {read_count}', status='summary')
    logger.info(f'Total unmapped reads with CB : {cb_count}', status='summary')
    logger.info(f'Total species and genus level classified Reads with CB : {use_count}', status='summary')
    logger.info(f'Skipped reads: {skipped}', status='summary')
    logger.info(f'Finishing getting the raw classified reads from bam file', status='complete')
    logger.info(f'Calculating quality control indicators', status='run')

    taxMap = dict()

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

    all_species_conf_metrics =[]
    if is_paired:
        # Convert the species_metrics_list to a DataFrame
        conf_columns = [
            'main_level_taxid', 'r1_confidence_score',
            'r1_rtl_score',"r1_host_score",'r2_confidence_score',
            'r2_rtl_score',"r2_host_score" 
        ]

        species_conf_df = pd.DataFrame(species_conf_list, columns=conf_columns)
        species_conf_df['mean_seq_confidence_score'] = species_conf_df[['r1_confidence_score', 'r2_confidence_score']].mean(axis=1)
        species_conf_df['mean_seq_rtl_score'] = species_conf_df[['r1_rtl_score', 'r2_rtl_score']].mean(axis=1)
        species_conf_df['mean_seq_host_score'] = species_conf_df[["r1_host_score", "r2_host_score"]].mean(axis=1)
        # 分组数据并计算分位数
        grouped = species_conf_df.groupby('main_level_taxid')
        q25 = grouped['mean_seq_confidence_score'].quantile(0.25)
        q75 = grouped['mean_seq_confidence_score'].quantile(0.75)

        # 创建包含分位数的新DataFrame
        species_conf_quantile_df = pd.DataFrame({'main_level_taxid': q25.index, 'mean_seq_confidence_score_q25': q25.values, 'mean_seq_confidence_score_q75': q75.values})

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

    else:
        # Convert the species_metrics_list to a DataFrame
        conf_columns = [
            'main_level_taxid', 'mean_seq_confidence_score',
            'mean_seq_rtl_score',"mean_seq_host_score"
        ]
        species_conf_df = pd.DataFrame(species_conf_list, columns=conf_columns)
        species_conf_df['mean_seq_confidence_score'] = species_conf_df['mean_seq_confidence_score'].astype(float)
        species_conf_df['mean_seq_rtl_score'] = species_conf_df['mean_seq_rtl_score'].astype(float)
        species_conf_df['mean_seq_host_score'] = species_conf_df['mean_seq_host_score'].astype(float)

        # Group data and calculate quantiles
        grouped = species_conf_df.groupby('main_level_taxid')
        q25 = grouped['mean_seq_confidence_score'].quantile(0.25)
        q75 = grouped['mean_seq_confidence_score'].quantile(0.75)

        # Create a new DataFrame containing quantiles
        species_conf_quantile_df = pd.DataFrame({'main_level_taxid': q25.index, 'mean_seq_confidence_score_q25': q25.values, 'mean_seq_confidence_score_q75': q75.values})

        # Repeat the steps to calculate quantiles for another column
        q25_rtl = grouped['mean_seq_rtl_score'].quantile(0.25)
        q75_rtl = grouped['mean_seq_rtl_score'].quantile(0.75)

        # Merge quantiles for another column into the new DataFrame
        species_conf_quantile_df['mean_seq_rtl_score_q25'] = q25_rtl.values
        species_conf_quantile_df['mean_seq_rtl_score_q75'] = q75_rtl.values

        q25_host = grouped['mean_seq_host_score'].quantile(0.25)
        q75_host = grouped['mean_seq_host_score'].quantile(0.75)
        species_conf_quantile_df['mean_seq_host_score_q25'] = q25_host.values
        species_conf_quantile_df['mean_seq_host_score_q75'] = q75_host.values

    # Add classification_rank and genus_level_taxid columns
    species_conf_quantile_df['classification_rank'] = species_conf_quantile_df.apply(lambda x: taxid2node[str(x['main_level_taxid'])].level_rank, axis=1)
    species_conf_quantile_df['genus_level_taxid'] = species_conf_quantile_df.apply(lambda x: taxid2node[str(x['main_level_taxid'])].taxid_to_desired_rank("G"), axis=1)

    # Calculate quantiles based on classification_rank
    species_conf_quantile_df['mean_seq_confidence_score_q25'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_confidence_score_q25'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_confidence_score_q25'].transform('max')
    )
    species_conf_quantile_df['mean_seq_rtl_score_q25'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_rtl_score_q25'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_rtl_score_q25'].transform('max')
    )
    species_conf_quantile_df['mean_seq_confidence_score_q75'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_confidence_score_q75'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_confidence_score_q75'].transform('max')
    )
    species_conf_quantile_df['mean_seq_rtl_score_q75'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_rtl_score_q75'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_rtl_score_q75'].transform('max')
    )
    species_conf_quantile_df['mean_seq_host_score_q25'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_host_score_q25'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_host_score_q25'].transform('mean')
    )
    species_conf_quantile_df['mean_seq_host_score_q75'] = np.where(
        species_conf_quantile_df['classification_rank'].str.startswith('S'),
        species_conf_quantile_df.groupby('main_level_taxid')['mean_seq_host_score_q75'].transform('mean'),
        species_conf_quantile_df.groupby('genus_level_taxid')['mean_seq_host_score_q75'].transform('mean')
    )
    # Create a list of dictionaries for
    data = [{"CB": cb, "main_level_taxid": main_level_taxid, "kmers": kmers["kmers"]} 
            for (cb, main_level_taxid), kmers in cb_taxid_to_ub_kmers.items()]

    # Create the DataFrame from the list of dictionaries
    cb_taxid_ub_kmer_count_df = pd.DataFrame(data)
    # Del data
    del data
    del cb_taxid_to_ub_kmers

    num_unique_CB = len(cb_taxid_ub_kmer_count_df['CB'].unique())
    if num_unique_CB > 300:

        cb_taxid_ub_kmer_count_df['genus_level_taxid'] = cb_taxid_ub_kmer_count_df.apply(lambda x: taxid2node[str(x['main_level_taxid'])].taxid_to_desired_rank("G"), axis=1)

        # Convert the DataFrame to long format, each row contains a kmer
        cb_taxid_ub_kmer_count_df = cb_taxid_ub_kmer_count_df.explode('kmers')

        # Calculate total kmer counts for each CB and species_level_taxid combination
        total_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'main_level_taxid']).size().reset_index(name='kmer_counts')

        # Calculate number of unique kmers for each CB and species_level_taxid combination 
        unique_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'main_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'unique_kmer_counts'})

        unique_genus_kmer_counts = cb_taxid_ub_kmer_count_df.groupby(['CB', 'genus_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'unique_genus_kmer_counts'})

        # 标识重复的 kmers，并获取不重复的行
        cb_taxid_ub_global_unique_count_df = cb_taxid_ub_kmer_count_df[~cb_taxid_ub_kmer_count_df.duplicated(subset=['kmers'], keep=False)]

        global_unique_kmer_counts =cb_taxid_ub_global_unique_count_df.groupby(['CB', 'main_level_taxid']).agg({'kmers': pd.Series.nunique}).reset_index().rename(columns={'kmers': 'global_unique_kmer_counts'})

        cb_taxid_kmer_count_df = pd.merge(total_kmer_counts, unique_kmer_counts, on=['CB', 'main_level_taxid'])
        cb_taxid_kmer_count_df = pd.merge(cb_taxid_kmer_count_df,global_unique_kmer_counts, on=['CB', 'main_level_taxid'])
        cb_taxid_kmer_count_df['genus_level_taxid'] = cb_taxid_kmer_count_df.apply(lambda x: taxid2node[str(x['main_level_taxid'])].taxid_to_desired_rank("G"), axis=1)
        cb_taxid_kmer_count_df = pd.merge(cb_taxid_kmer_count_df,unique_genus_kmer_counts, on=['CB', 'genus_level_taxid'])
        cb_taxid_kmer_count_df['classification_rank'] = cb_taxid_kmer_count_df.apply(lambda x: taxid2node[str(x['main_level_taxid'])].level_rank, axis=1)

        del cb_taxid_ub_kmer_count_df

        # Filter rows with classification_rank as S
        cb_species_kmer_count_df = cb_taxid_kmer_count_df[cb_taxid_kmer_count_df['classification_rank'] == 'S']

        # Initialize an empty list to store DataFrames for each group
        dfs = []

        grouped = cb_species_kmer_count_df.groupby('main_level_taxid') 
        for name, group in grouped:
            if len(group) > 1:
                corr_kmer_uniq, p_kmer_uniq = spearmanr(group['kmer_counts'], group['unique_kmer_counts'])
                corr_kmer_glob_uniq, p_kmer_glob_uniq = spearmanr(group['kmer_counts'], group['global_unique_kmer_counts'])
                # Create a DataFrame for the current group
                group_df = pd.DataFrame({
                    'main_level_taxid': [name],
                    'corr_kmer_uniq': [corr_kmer_uniq],
                    'p_value_kmer_uniq': [p_kmer_uniq],
                    'corr_kmer_glob_uniq': [corr_kmer_glob_uniq],
                    'p_value_kmer_glob_uniq': [p_kmer_glob_uniq]
                })
                
                # Append the DataFrame to the list
                dfs.append(group_df)
            else:
                corr_kmer_uniq, p_kmer_uniq, corr_kmer_glob_uniq, p_kmer_glob_uniq = np.nan, np.nan, np.nan, np.nan

                # Create a DataFrame for the current group
                group_df = pd.DataFrame({
                    'main_level_taxid': [name],
                    'corr_kmer_uniq': [corr_kmer_uniq],
                    'p_value_kmer_uniq': [p_kmer_uniq],
                    'corr_kmer_glob_uniq': [corr_kmer_glob_uniq],
                    'p_value_kmer_glob_uniq': [p_kmer_glob_uniq]
                })
                
                # Append the DataFrame to the list
                dfs.append(group_df)
                
        cb_taxid_kmer_count_df['genus_kmer_counts'] = cb_taxid_kmer_count_df.groupby(['CB','genus_level_taxid'])['kmer_counts'].transform('sum')
        cb_taxid_kmer_count_df['genus_global_unique_kmer_counts'] = cb_taxid_kmer_count_df.groupby(['CB','genus_level_taxid'])['global_unique_kmer_counts'].transform('sum')
        cb_taxid_kmer_count_df
        # Filter rows with classification_rank as S
        cb_genus_kmer_count_df = cb_taxid_kmer_count_df[cb_taxid_kmer_count_df['classification_rank'] == 'G']
        grouped = cb_genus_kmer_count_df.groupby('main_level_taxid') 
        for name, group in grouped:
            if len(group) > 1:
                corr_kmer_uniq, p_kmer_uniq = spearmanr(group['genus_kmer_counts'], group['unique_genus_kmer_counts'])
                corr_kmer_glob_uniq, p_kmer_glob_uniq = spearmanr(group['genus_kmer_counts'], group['genus_global_unique_kmer_counts'])
                # Create a DataFrame for the current group
                group_df = pd.DataFrame({
                    'main_level_taxid': [name],
                    'corr_kmer_uniq': [corr_kmer_uniq],
                    'p_value_kmer_uniq': [p_kmer_uniq],
                    'corr_kmer_glob_uniq': [corr_kmer_glob_uniq],
                    'p_value_kmer_glob_uniq': [p_kmer_glob_uniq]
                })
                
                # Append the DataFrame to the list
                dfs.append(group_df)
            else:
                corr_kmer_uniq, p_kmer_uniq, corr_kmer_glob_uniq, p_kmer_glob_uniq = np.nan, np.nan, np.nan, np.nan

                # Create a DataFrame for the current group
                group_df = pd.DataFrame({
                    'main_level_taxid': [name],
                    'corr_kmer_uniq': [corr_kmer_uniq],
                    'p_value_kmer_uniq': [p_kmer_uniq],
                    'corr_kmer_glob_uniq': [corr_kmer_glob_uniq],
                    'p_value_kmer_glob_uniq': [p_kmer_glob_uniq]
                })
                
                # Append the DataFrame to the list
                dfs.append(group_df)
                
        # Concatenate all DataFrames in the list
        cb_taxid_kmer_corr_df = pd.concat(dfs, ignore_index=True)

        cluster_df = pd.read_csv(args.cluster,sep="\t",)
        cluster_df.columns = ["barcode","leiden"]
        cb_cluster_taxid_kmer_count_df = pd.merge(cb_taxid_kmer_count_df, cluster_df, left_on='CB', right_on='barcode', how='inner')
        leiden_counts = cluster_df.groupby('leiden').size().reset_index(name='cluster_sum')
        cb_cluster_taxid_kmer_count_df['leiden'] = cb_cluster_taxid_kmer_count_df['leiden'].astype(str)
        total_unique_CB = cb_cluster_taxid_kmer_count_df['CB'].nunique()
        total_unique_cluster = cb_cluster_taxid_kmer_count_df['leiden'].nunique()
        # species_prevalence = cb_cluster_taxid_kmer_count_df.groupby('species_level_taxid')['CB'].nunique().reset_index()
        # # Rename the columns for clarity
        # species_prevalence.columns = ['species_level_taxid', 'unique_CB_count']
        # species_prevalence['CB_prevalence'] = species_prevalence['unique_CB_count'] / total_unique_CB
        # cb_cluster_taxid_kmer_count_df.to_csv("/data/scRNA_analysis/benchmark/Galeano2022_GEX/notebooks/cb_cluster_taxid_kmer_count_df.csv",index=False)
        # Calculate the prevalence for each species_level_taxid in CBs
        cb_cluster_species_taxid_kmer_count_df = cb_cluster_taxid_kmer_count_df[cb_cluster_taxid_kmer_count_df['classification_rank'] == 'S']

        # Calculate unique CB counts at species level
        species_cb_counts = cb_cluster_species_taxid_kmer_count_df.groupby('main_level_taxid')['CB'].nunique()

        # # Calculate unique valid CB counts at species level
        # valid_species_cb_counts = cb_cluster_species_taxid_kmer_count_df[cb_cluster_species_taxid_kmer_count_df['valid'] == True].groupby('main_level_taxid')['CB'].nunique()

        # Calculate unique cluster counts at species level
        species_cluster_counts = cb_cluster_species_taxid_kmer_count_df.groupby(['main_level_taxid', 'leiden']).size().reset_index(name='count')

        # Calculate unique CB counts at genus level
        genus_cb_counts = cb_cluster_taxid_kmer_count_df.groupby('genus_level_taxid')['CB'].nunique()

        # # Calculate unique valid CB counts at genus level
        # valid_genus_cb_counts = cb_cluster_taxid_kmer_count_df[cb_cluster_taxid_kmer_count_df['valid'] == True].groupby('genus_level_taxid')['CB'].nunique()

        # Calculate unique cluster counts at genus level
        genus_cluster_counts = cb_cluster_taxid_kmer_count_df.groupby(['genus_level_taxid', 'leiden']).size().reset_index(name='count')

        # Place the two Series data into DataFrames
        df_species = species_cb_counts.reset_index()
        df_genus = genus_cb_counts.reset_index()

        # Rename columns in DataFrames
        df_species.columns = ['main_level_taxid', 'unique_CB_count']
        df_genus.columns = ['main_level_taxid', 'unique_CB_count']
        genus_cluster_counts.columns = ['main_level_taxid', 'leiden', 'unique_cluster_count']
        species_cluster_counts.columns = ['main_level_taxid', 'leiden', 'unique_cluster_count']

        # Merge the two DataFrames
        taxa_cb_prevalence = pd.concat([df_species, df_genus], ignore_index=True)
        taxa_cb_prevalence['CB_prevalence'] = taxa_cb_prevalence['unique_CB_count'] / total_unique_CB
        taxa_cluster_counts = pd.concat([species_cluster_counts, genus_cluster_counts], ignore_index=True)
        # Add leiden clsuter sum to taxa_cluster_counts dataframe
        leiden_counts['leiden'] = leiden_counts['leiden'].astype(str)
        taxa_cluster_counts['leiden'] = taxa_cluster_counts['leiden'].astype(str)
        taxa_cluster_counts = pd.merge(taxa_cluster_counts, leiden_counts, on='leiden', how='left')

        # Calculate the min value of count * 0.01 and 3 as the threshold
        taxa_cluster_counts['threshold'] = np.minimum(taxa_cluster_counts['cluster_sum'] * 0.01, 3)

        # Filter the data based on unique_cluster_count and the calculated threshold
        taxa_cluster_filtered = taxa_cluster_counts[taxa_cluster_counts['unique_cluster_count'] >= taxa_cluster_counts['threshold']]

        taxa_cluster_prevalence = taxa_cluster_filtered.groupby('main_level_taxid')['leiden'].nunique().reset_index()
        taxa_cluster_prevalence.columns = ['main_level_taxid', 'unique_cluster_count']

        taxa_cluster_prevalence['cluster_prevalence'] = taxa_cluster_prevalence['unique_cluster_count'] / total_unique_cluster

        taxa_prevalence_combined = pd.merge(taxa_cb_prevalence, taxa_cluster_prevalence, on='main_level_taxid', how ='left')

        taxa_prevalence_combined['main_level_taxid'] = taxa_prevalence_combined['main_level_taxid'].astype(str)

        del cb_cluster_taxid_kmer_count_df
        del cb_taxid_ub_global_unique_count_df


        # non_na_rows = cb_taxid_kmer_corr_df[p_val_cols].notna().any(axis=1)

        # Calculate ntests using non-na rows
        # ntests = non_na_rows.sum()
        # Perform multiple testing correction only if ntests is non-zero
        # if ntests > 5:
        #     # Filter out NaN p-values and adjust them for both sets
        #     for col in p_val_cols:
        #         if col in cb_taxid_kmer_corr_df.columns:
        #             non_nan = cb_taxid_kmer_corr_df[col].notna()
        #             cb_taxid_kmer_corr_df.loc[non_nan, col] = multipletests(cb_taxid_kmer_corr_df.loc[non_nan, col], method='fdr_bh')[1]
        # else:
        #     pass
    else:
        # ntests = 0
        pass

    

    final_desired_krak_report = desired_krak_report.copy()
    # Convert 'ncbi_taxa' column to string data type
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    # final_desired_krak_report.drop('fraction', axis=1, inplace=True)
    final_desired_krak_report['cov'].replace([float('inf'), float('-inf')], float('nan'), inplace=True)
    final_desired_krak_report['max_cov'] = np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['cov'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['cov'].transform('max')
    )
    final_desired_krak_report['max_uniqminimizers'] = np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['uniqminimizers'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['uniqminimizers'].transform('max')
    )
    final_desired_krak_report['max_minimizers'] = np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['minimizers'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['minimizers'].transform('max')
    )
    final_desired_krak_report['max_minimizers'] = np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['minimizers'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['minimizers'].transform('max')
    )

    ## Reset index
    final_desired_krak_report.reset_index(drop=True, inplace=True)
    # Merge species_seq_metrics to final_desired_krak_report
    final_desired_krak_report = final_desired_krak_report.merge(
        species_seq_metrics[
            ["average_kmer_consistency", "average_seq_entropy", "average_seq_dust_score", "max_seq_length", "mean_seq_length", "main_level_taxid"]
        ],
        left_on='main_level_taxid', 
        right_on='main_level_taxid'
    )    
    # Merge species_conf_quantile_df to final_desired_krak_report
    final_desired_krak_report = final_desired_krak_report.merge(
        species_conf_quantile_df[
            ["mean_seq_confidence_score_q25", "mean_seq_confidence_score_q75", "mean_seq_rtl_score_q25", "mean_seq_rtl_score_q75", "mean_seq_host_score_q25", "mean_seq_host_score_q75", "main_level_taxid"]
        ],
        left_on='main_level_taxid', 
        right_on='main_level_taxid'
    )
    # Set ncbi taxa as str object
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    taxa_nucleotides_df['ncbi_taxa'] = taxa_nucleotides_df['ncbi_taxa'].astype(str)
    final_desired_krak_report = final_desired_krak_report.merge(taxa_nucleotides_df, left_on='ncbi_taxa',right_on='ncbi_taxa', how='left')
    # Get the max nucleotides
    final_desired_krak_report['max_nucleotides'] = np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['nucleotides'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['nucleotides'].transform('max')
    )

    # select the unique taxa
    taxa_nucleotides_unique = taxa_nucleotides_df[~taxa_nucleotides_df['ncbi_taxa'].isin(final_desired_krak_report['ncbi_taxa'])]
    taxa_nucleotides_unique.loc[:, "classification_rank"] = taxa_nucleotides_unique.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].level_rank, axis=1)
    taxa_nucleotides_unique.loc[:,"scientific name"] = taxa_nucleotides_unique.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].name, axis=1)
    taxa_nucleotides_unique = taxa_nucleotides_unique[taxa_nucleotides_unique['classification_rank'].str.startswith('S')]
    taxa_nucleotides_unique["genus_level_taxid"] = taxa_nucleotides_unique.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("G"), axis=1)
    taxa_nucleotides_unique['main_level_taxid'] = taxa_nucleotides_unique.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].get_main_lvl_taxid(), axis=1)
    taxa_nucleotides_unique['is_microbiome'] = taxa_nucleotides_unique.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].is_microbiome(), axis=1)

    # final_desired_krak_report.drop('taxid', axis=1, inplace=True)
    if num_unique_CB > 300:
        final_desired_krak_report = final_desired_krak_report.merge(cb_taxid_kmer_corr_df,left_on='main_level_taxid', right_on='main_level_taxid')
        final_desired_krak_report = final_desired_krak_report.merge(taxa_prevalence_combined,left_on='main_level_taxid', right_on='main_level_taxid')
        
    # final_desired_krak_report.drop('main_level_taxid', axis=1, inplace=True)
    else:
        pass

    logger.info(f'Finishging calculating quality control indicators', status='complete')

    num_unique_species = len(final_desired_krak_report['ncbi_taxa'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids having qc indictor', status='summary')

    # Save data
    logger.info(f'Saving the raw result', status='run')
    final_desired_krak_report.to_csv(args.raw_qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')

    logger.info(f'Filtering taxa with quality control indicators', status='run')
    final_desired_krak_report['superkingdom'] = final_desired_krak_report['superkingdom'].astype(str)

    bac_cov_cutoff = total_reads/10000000*0.0001
    viral_cov_cutoff = total_reads/10000000*0.0005
    
    ## For many corr
    if  num_unique_CB > 300:
        filter_desired_krak_report = final_desired_krak_report.copy()[
            (
                (
                (final_desired_krak_report['max_minimizers'] > 5) &
                (
                    ((final_desired_krak_report['cluster_prevalence'] < 0.8) &
                (final_desired_krak_report['CB_prevalence'] < 0.6)) |
                    ((final_desired_krak_report['max_cov'] >= 2) &(final_desired_krak_report['max_uniqminimizers'] >= 2000))
                )&
                (final_desired_krak_report['unique_CB_count'] >= 3) &
                (
                ( final_desired_krak_report['dup'] <= 80) | (final_desired_krak_report['max_uniqminimizers'] >= 100000 ) 
                ) &
                # (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
                # (final_desired_krak_report['mean_seq_confidence_score_q75'] > 0.3) &
                (
                    (
                        ((final_desired_krak_report['superkingdom'] == '2') &
                            (
                                ((final_desired_krak_report['max_cov'] >= bac_cov_cutoff)) 
                            | 
                                 ((final_desired_krak_report['average_seq_entropy'] >= 1.8 ))
                            )
                        )                         
                        |
                        ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov'] >=  bac_cov_cutoff)) 
                        
                        |
                        ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_cov'] >=  bac_cov_cutoff)) 
                        |
                        ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] > viral_cov_cutoff) & (final_desired_krak_report['average_seq_entropy'] >= 1.2 )) 
                    )
                ) &
                (
                    ~(
                        (final_desired_krak_report['corr_kmer_uniq'] < 0) &
                        (final_desired_krak_report['p_value_kmer_uniq'] > float(0.05)) 
                        # (final_desired_krak_report['corr_kmer_glob_uniq'] > 0 ) &
                        # (final_desired_krak_report['p_value_kmer_glob_uniq'] <float(0.05))
                    )
                )
                )
            )
            ]
    else:
        filter_desired_krak_report = final_desired_krak_report.copy()[
            (final_desired_krak_report['average_seq_entropy'] > 1.4) &
            (final_desired_krak_report['max_minimizers'] > 5) &
            (
               ( final_desired_krak_report['dup'] <= 80) | (final_desired_krak_report['max_uniqminimizers'] >= 100000 ) 
            ) &
            # (final_desired_krak_report['average_seq_dust_score'] > args.min_dust) &
            (
                (
                    ((final_desired_krak_report['superkingdom'] == '2') &
                            (
                                ((final_desired_krak_report['max_cov'] >= bac_cov_cutoff)) 
                            # | 
                            #     # ((final_desired_krak_report['average_seq_entropy'] >= 1.8 ) & (final_desired_krak_report['max_uniqminimizers'] >= 300 ) &(final_desired_krak_report['average_seq_dust_score'] >= 0.08 ) &(final_desired_krak_report['max_minimizers'] >= total_reads*0.01 ))
                            #      ((final_desired_krak_report['average_seq_entropy'] >= 1.85 ))
                            )) |
                    ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_cov'] >=  bac_cov_cutoff)) |
                    ((final_desired_krak_report['superkingdom'] == '2759') &  (final_desired_krak_report['max_cov'] >=  bac_cov_cutoff)) |
                    ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_cov'] > viral_cov_cutoff)) 
                )
            ) 
        ]
    # filter_desired_krak_report.drop(['frac','classification_rank','fraction','minimizers_clade','minimizers_taxa','ncbi_taxa','sci_name','cov','species_level_taxa','level_1'], axis=1, inplace=True)
    filter_desired_krak_report['scientific name'] = filter_desired_krak_report['scientific name'].apply(lambda x: x.strip())
    recall_cutoff = filter_desired_krak_report[filter_desired_krak_report['classification_rank'].str.startswith('S')]["max_nucleotides"].quantile(0.05)
    recall = taxa_nucleotides_unique[(taxa_nucleotides_unique["nucleotides"]>recall_cutoff) & (taxa_nucleotides_unique["is_microbiome"] ==True)].head(20)
    filter_desired_krak_report = pd.concat([filter_desired_krak_report, recall], axis=0)
    # # Filter out rows where 'ncbi_taxa' matches any value from 'excluded_taxonomy_ids'
    # filter_desired_krak_report = filter_desired_krak_report[~filter_desired_krak_report['ncbi_taxa'].isin(args.exclude)]

    logger.info(f'Finishing filtering taxa with quality control indicators', status='complete')
    num_unique_species = len(filter_desired_krak_report['ncbi_taxa'].unique())
    logger.info(f'After filtering, found {num_unique_species} unique species and subspeceis level taxids', status='summary')
    num_unique_species = len(filter_desired_krak_report['species_level_taxid'].unique())
    num_unique_genus = len(filter_desired_krak_report['genus_level_taxid'].unique())

    logger.info(f'After filtering, found {num_unique_species} unique species level taxids and {num_unique_genus} unique genus level taxids', status='summary')

    # Save data
    logger.info(f'Saving the result', status='run')
    filter_desired_krak_report.to_csv(args.qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')


if __name__ == "__main__":
    main()