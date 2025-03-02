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

#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--krak_report', required=True, 
        dest="krak_report_file", help='Input kraken report file for denosing')
    parser.add_argument('--krak_output', required=True,
        dest='krak_output_file', help='Input kraken output file for denosing')
    parser.add_argument('--raw_qc_output_file', required=True,
        help='Output denosed info at individual level')
    parser.add_argument('--qc_output_file', required=True,
        help='Output denosed info at individual level')
    parser.add_argument('--ktaxonomy', required=True,
        help='Kraken2 database ktaxonomy file path')
    parser.add_argument('--inspect', required=True,
        dest="inspect_file", help='Kraken2 database inspect file path')
    parser.add_argument('--exclude', required=False,
        default=9606, nargs='+',
        help='Taxonomy ID[s] of reads to exclude (space-delimited)')
    parser.add_argument('--nsample', required=False,
        default=2500,
        help='Max number of reads to sample per taxa')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
    parser.add_argument('--min_read_fraction', required=False,
        default=0.15, type=float, help='Minimum fraction of kmers directly assigned to taxid [default=0.15]')
    
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

    # filter kraken_file to species only
    desired_krak_report = krak_report.copy()[krak_report['classification_rank'].str.startswith(('S'), na=False)]
    desired_krak_report['species_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("S"), axis=1)
    desired_krak_report['main_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].get_main_lvl_taxid(), axis=1)
    desired_krak_report['genus_level_taxid'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("G"), axis=1)
    desired_krak_report['superkingdom'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].taxid_to_desired_rank("D"), axis=1)
    desired_krak_report['is_microbiome'] = desired_krak_report.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].is_microbiome(), axis=1)

    ## select microbiome
    desired_krak_report = desired_krak_report[desired_krak_report["is_microbiome"]==True]
    # Transform data type
    desired_krak_report['species_level_taxid'] = desired_krak_report['species_level_taxid'].astype(str)
    desired_krak_report['ncbi_taxa'] = desired_krak_report['ncbi_taxa'].astype(str)
    
    # desired_krak_report
    desired_taxid_list = set(desired_krak_report['ncbi_taxa'].unique())
    desired_main_taxid_list = set(desired_krak_report['main_level_taxid'].unique())
    logger.info('Finished processing kraken2 classifier result', status='complete')
    # del df
    del krak_report

    # lineage_dict = {}
    # for main_tax_id in desired_main_taxid_list:
    #     try:
    #         lineage_taxid_list = taxid2node[main_tax_id].lineage_to_desired_rank("D")
    #         lineage_dict[main_tax_id] = lineage_taxid_list
    #     except (ValueError, KeyError) as e:
    #         print("Error occur:", e)
    descendants_dict = {}
    for tax_id in desired_taxid_list:
        try:
            descendants_taxid_list = taxid2node[tax_id].get_all_descendants()
            descendants_dict[tax_id] = descendants_taxid_list
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

    # Auto check file read type
    is_paired = None
    with open(args.krak_output_file, 'r') as kfile:
        for i, kraken_line in enumerate(kfile):
            if i >= 5:  # only check first 5 lines
                break
            try:
                read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                current_paired = "|:|" in kmer_position
                if is_paired is None:
                    is_paired = current_paired
                elif is_paired != current_paired:
                    raise ValueError("Inconsistent paired-end status detected in first 5 lines, please check the kraken output file and fastq file")
            except Exception as e:
                logger.error(f"Error processing line {i+1}: {e}")
                raise
    if is_paired is None:
        raise ValueError("Could not determine paired-end status from first 5 lines, please check the kraken output file and fastq file")

    # Reading kraken2 classifier output information
    logger.info('Reading kraken2 classifier output information', status='run')
    taxid_counts = {}
    kraken_data = {}
    taxid_fraction_counts = defaultdict(lambda: defaultdict(int))
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
                # only select desired taxid
                if (tax_id in desired_taxid_list):
                    r1_len, r2_len = read_len.split('|')
                    r1_kmer_position, r2_kmer_position  = kmer_position.split(' |:| ')
                    if tax_id not in taxid_counts:
                        taxid_counts[tax_id] = 1
                    else:
                        taxid_counts[tax_id] += 1
                    for kmer in r1_kmer_position.split() + r2_kmer_position.split():
                        kmer_taxid, count = kmer.split(':')
                        taxid_fraction_counts[tax_id][kmer_taxid] += int(count)
                    # if taxid_counts[tax_id] >= args.nsample:
                    #     continue 
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
                # only select desired taxid
                if (tax_id in desired_taxid_list):
                    if tax_id not in taxid_counts:
                        taxid_counts[tax_id] = 1
                    else:
                        taxid_counts[tax_id] += 1
                    for kmer in kmer_position.split():
                        kmer_taxid, count = kmer.split(':')
                        taxid_fraction_counts[tax_id][kmer_taxid] += int(count)
                    # if taxid_counts[tax_id] >= args.nsample:
                    #     continue 
                else:
                    continue
        
    logger.info("Calculating taxid counts ratio...", status='run')
    taxid_counts_ratio = {}
    for target_taxid in desired_taxid_list:
        sorted_taxids = sorted(taxid_fraction_counts[target_taxid].items(), key=lambda x: x[1], reverse=True)[:5]
        # calculate total kmer count
        total_kmers = sum(kmer_count for _, kmer_count in sorted_taxids)

        target_taxid_descendants = descendants_dict[target_taxid]
        # get kmer count of taxid, solve the case that taxid is not in descendants_dict
        target_kmer = next((kmer_count for taxid, kmer_count in sorted_taxids if taxid in target_taxid_descendants), 0)

        # calculate ratio
        if target_kmer > 0:
            ratio = target_kmer / total_kmers
        else:
            ratio = 0

        taxid_counts_ratio[target_taxid] = ratio

    logger.info('Finishing reading kraken2 classifier output information', status='complete')

    # Get species level taxid
    num_unique_species = len(desired_krak_report['species_level_taxid'].unique())
    num_unique_genus = len(desired_krak_report['genus_level_taxid'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids and {num_unique_genus} unique genus level taxids', status='summary')    

    # copy desired_krak_report as final_desired_krak_report
    final_desired_krak_report = desired_krak_report.copy()
    # Convert 'ncbi_taxa' column to string data type
    final_desired_krak_report['ncbi_taxa'] = final_desired_krak_report['ncbi_taxa'].astype(str)
    final_desired_krak_report['read_fraction'] = final_desired_krak_report['ncbi_taxa'].map(taxid_counts_ratio)
    final_desired_krak_report['cov'].replace([float('inf'), float('-inf')], float('nan'), inplace=True)

    # calculate max cov, max uniqminimizers, max minimizers, max read fraction
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
    final_desired_krak_report['max_read_fraction'] =  np.where(
        final_desired_krak_report['classification_rank'].str.startswith('S'),
        final_desired_krak_report.groupby('main_level_taxid')['read_fraction'].transform('max'),
        final_desired_krak_report.groupby('genus_level_taxid')['read_fraction'].transform('max')
    )

    logger.info(f'Finishging calculating quality control indicators', status='complete')

    num_unique_species = len(final_desired_krak_report['ncbi_taxa'].unique())
    logger.info(f'Found {num_unique_species} unique species level taxids having qc indictor', status='summary')

    # save raw result
    logger.info(f'Saving the raw result', status='run')
    final_desired_krak_report.to_csv(args.raw_qc_output_file, sep="\t", index=False)
    logger.info(f'Finishing saving the result', status='complete')

    logger.info(f'Filtering taxa with quality control indicators', status='run')
    final_desired_krak_report['superkingdom'] = final_desired_krak_report['superkingdom'].astype(str)

    # get min_read_fraction from args
    min_read_fraction = args.min_read_fraction
    
    # filter desired_krak_report
    # TODO: each superkingdom have different min_read_fraction
    # TODO: MAX COVERAGE
    filter_desired_krak_report = final_desired_krak_report.copy()[
            (
            (final_desired_krak_report['max_minimizers'] > 5) &
            (final_desired_krak_report['dup'] < 100) &
            (
                (
                    ((final_desired_krak_report['superkingdom'] == '2') &
                        (
                            ((final_desired_krak_report['max_read_fraction'] >= min_read_fraction)) 
                        )
                    )                         
                    |
                    ((final_desired_krak_report['superkingdom'] == '2157')& (final_desired_krak_report['max_read_fraction'] >= min_read_fraction)) 

                    |
                    ((final_desired_krak_report['superkingdom'] == '2759') & (final_desired_krak_report['max_read_fraction'] >= min_read_fraction)) 
                    |
                    ((final_desired_krak_report['superkingdom'] == '10239') & (final_desired_krak_report['max_read_fraction'] >= min_read_fraction)) 
                )
            )
        )
        ]

    # remove space
    filter_desired_krak_report['scientific name'] = filter_desired_krak_report['scientific name'].apply(lambda x: x.strip())

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