import pysam
import sys
import gzip
import argparse
import numpy as np
import logging
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import pandas as pd
import collections.abc
import os
import csv
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
    
def most_frequent(List):
    """Finds the most frequent element in a list"""
    return max(set(List), key = List.count)

def map_nested_dicts(ob, func):
    """ Applys a map to the inner item of nested dictionaries """
    for k, v in ob.items():
        if isinstance(v, collections.abc.Mapping):
            map_nested_dicts(v, func)
        else:
            ob[k] = func(v)

def twist_dict_UMI(nested, taxid2node):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = {}
    # Process cell-transcript relationships
    for ckey, tdict in nested.items():
        # for UMI situation
        for tkey, kvalue in tdict.items():
            if ckey in newdict:
                if kvalue in newdict[ckey]:
                    newdict[ckey][kvalue] += 1
                else:
                    newdict[ckey][kvalue] = 1
            else:
                newdict[ckey] = {kvalue: 1}

    # Process genus-species relationships
    for ckey, tdict in newdict.items():
        for kvalue, count in list(tdict.items()):
            # Get the lineage using taxid2node
            if kvalue in taxid2node:
                lineage = taxid2node[kvalue].get_taxon_path()[0].split('|')
                # Add count to parent if it exists
                for parent in lineage:
                    if parent != kvalue and parent in taxid2node:
                        if parent in newdict[ckey]:
                            newdict[ckey][parent] += count
                        else:
                            newdict[ckey][parent] = count
            
    return(newdict)

def twist_dict(nested, taxid2node):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = nested.copy()  # Create a copy to avoid modifying during iteration
    
    # Process genus-species relationships
    for ckey, tdict in newdict.items():
        for kvalue, count in list(tdict.items()):
            # Get the lineage using taxid2node
            if kvalue in taxid2node:
                lineage = taxid2node[kvalue].get_taxon_path()[0].split('|')
                # Add count to parent if it exists
                for parent in lineage:
                    if parent != kvalue and parent in taxid2node:
                        if parent in newdict[ckey]:
                            newdict[ckey][parent] += count
                        else:
                            newdict[ckey][parent] = count
            
    return(newdict)

def dict2lists(nested):
    """ Returns lists for sparse matrix """
    rows = [] # cell coordinate
    columns = [] # taxonomy id coordinate
    values = [] # count

    cell_list = [] # same order as rows
    taxid_list = [] # same order as columns

    j = 0

    for ckey, taxdict in nested.items():
        for taxkey, count in taxdict.items():
            try:
                k = taxid_list.index(taxkey)
            except:
                taxid_list.append(taxkey)
                k = taxid_list.index(taxkey)
                
            rows.append(k)
            columns.append(j)
            values.append(count) 
            
        # increase cell coordinate by 1
        cell_list.append(ckey)
        j += 1
    
    return rows, columns, values, cell_list, taxid_list

#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser(description='This script is used to output bam classified microbial data in cellranger format as feature.tsv,barcodes.tsv,matrix.mtx \n This requires the additional packages pysam(If your python version is up to 3.9)\n')
    parser.add_argument('--cb_bam', help='Input align SAM or BAM file with CB', required=True)
    parser.add_argument('--align_result', help='Input align tsv containing read name, taxonomy id', required=True)
    parser.add_argument('--ktaxonomy_file', help='Input kraken2 ktaxonomy file', required=True)
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
    parser.add_argument('--profile_tsv', help='Output microbiome read tsv', required=True)
    parser.add_argument('--matrixfile', help='Output microbiome matrix', required=True)
    parser.add_argument('--cellfile', help='Output cell barcodes', required=True)
    parser.add_argument('--taxfile', help='Output taxonomy IDs', required=True)
    parser.add_argument('--log_file', dest='log_file', required=True, help="File to write the log to")
    parser.add_argument('--output_read_tsv', help='Output read assignment tsv with barcode information', required=False)
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

    # Loading taxonomy database
    logger.info('Loading taxonomy database from Kraken taxonomy file', status='run')
    try:
        taxid2node = make_dicts(args.ktaxonomy_file)
        logger.info('Successfully loaded taxonomy database', status='complete')
    except Exception as e:
        logger.error(f"An error occurred while loading taxonomy database: {e}")
        sys.exit()

    logger.info('Prasing bam file', status='run')
    total_count = 0
    use_count = 0
    failed_count = 0
    read_taxid_info_dict = dict()
    taxid_info_dict = dict()
    with open(args.align_result, "r") as taxa_file:
        # Skip header line if present
        header = next(taxa_file, None)
        if header is None:
            logger.warning(f"Alignment result file {args.align_result} is empty.", status='warning')
            # Optionally exit if the file is completely empty
            # sys.exit(f"Error: Alignment result file {args.align_result} is empty.")

        for line in taxa_file:
            total_count += 1
            try:
                read_name, tax_id, prob, tax_name, tax_genus = line.strip().split("\t")
            except:
                failed_count += 1
                continue
            if tax_genus != "Unknown":
                read_taxid_info_dict[read_name] = {"taxid": tax_id, "taxname": tax_name, "taxgenus": tax_genus}
                taxid_info_dict[tax_id] = tax_name

                use_count +=1
    # Check if any usable reads were found
    if use_count == 0:
        logger.error(f"No usable alignment results found in {args.align_result}. Cannot proceed.", status='error')
        # Create empty output files to satisfy downstream dependencies if necessary
        logger.info("Creating empty output files.", status='run')
        open(args.matrixfile, 'w').close()
        open(args.cellfile, 'w').close()
        open(args.taxfile, 'w').close()
        open(args.profile_tsv, 'w').close()
        if args.output_read_tsv:
             open(args.output_read_tsv, 'w').close()
        logger.info("Empty output files created. Exiting.", status='complete')
        sys.exit(0) # Exit gracefully
    logger.info(f'Prasing bam file complete, total reads: {total_count}, use reads: {use_count}, failed reads: {failed_count}, failed rate: {failed_count/total_count}', status='complete')

    logger.info('Checking barcode bam file type', status='run')
    is_cb = False
    is_ub = False
    is_rg = False
    read_count = 0
    with pysam.AlignmentFile(args.cb_bam, "rb", check_sq = False) as barcode_file:

        for bread in barcode_file:
            # 检查是否存在目标标签
            if bread.has_tag("CB"):
                is_cb = True
            if bread.has_tag("UB"):
                is_ub = True
            if bread.has_tag("RG"):
                is_rg = True
            # if didnt find type in first 250 reads, exit
            read_count += 1
            if read_count >= 250:
                break
    
    if is_cb and is_ub :
        logger.info('Detect cellbarcode and UMI tag, use UMI identity', status='complete')
        mode = "CB_UMI"
    elif is_cb and not is_ub:
        logger.info('Only Detect cellbarcode tag, use cellbarcode identity', status='complete')
        mode = "CB"
    elif is_rg and not is_cb and not is_ub:
        logger.info('Only Detect read group tag, use read group identity', status='complete')
        mode = "RG"
    else:
        logger.error('Didnt detect any read identity, exit', status='error')
        sys.exit()
    
    logger.info("Parsing bam file", status='run')
    total_count = 0
    use_count = 0
    skipped = 0
    if mode == "CB_UMI":
        # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
        nested_dict = {}
        with pysam.AlignmentFile(args.cb_bam, "rb", check_sq = False) as cb_file:

            for bread in cb_file:
                total_count += 1
                # Check if the read exists in the kraken file
                if bread.query_name not in read_taxid_info_dict:
                    skipped += 1
                    # logging.warning("Read name {} not found in kraken file".format(sread.query_name))
                    continue    
                try:
                    bread_CB = bread.get_tag("CB")
                    bread_UB = bread.get_tag("UB")
                except:
                    skipped += 1
                    continue
                
                read_taxid = read_taxid_info_dict[bread.query_name]["taxid"]
                # add barcode information to read_taxid_info_dict
                read_taxid_info_dict[bread.query_name]["Barcode"] = bread_CB
                
                # Make nested dictionary with cells and transcripts
                # {cellbarcode: {transcriptbarcode: krakentaxonomyID}
                if bread_CB in nested_dict:
                    # If cell and transcript exist, add taxonomy ID to list
                    if bread_UB in nested_dict[bread_CB]:
                        nested_dict[bread_CB][bread_UB].append(read_taxid)
                    # Otherwise create transcript dictionary for cell
                    else:
                        nested_dict[bread_CB][bread_UB] = [read_taxid]
                else:
                    # if cell doesn't exist, create cell and transcript dictionary with kraken id
                    nested_dict[bread_CB] = {bread_UB: [read_taxid]}

                use_count += 1
        # Find most frequent taxonomy for each transcript
        map_nested_dicts(nested_dict, most_frequent)
        # Make sparse matrix
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict_UMI(nested_dict, taxid2node))
    if mode == "CB":
        # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
        nested_dict = {}
        with pysam.AlignmentFile(args.cb_bam, "rb",check_sq = False) as cb_file:

            for bread in cb_file:
                total_count += 1
                # Check if the read exists in the kraken file
                if bread.query_name not in read_taxid_info_dict:
                    skipped += 1
                    # logging.warning("Read name {} not found in kraken file".format(sread.query_name))
                    continue    
                try:
                    bread_CB = bread.get_tag("CB")
                except:
                    skipped += 1
                    continue
                
                read_taxid = read_taxid_info_dict[bread.query_name]["taxid"]
                # add barcode information to read_taxid_info_dict
                read_taxid_info_dict[bread.query_name]["Barcode"] = bread_CB
                
                # Make nested dictionary with RG and taxonomy IDs
                # {cellbarcode: {taxonomyID}
                # If CB exists, add taxonomy ID to list 
                if bread_CB in nested_dict:
                    if read_taxid in nested_dict[bread_CB]:
                        nested_dict[bread_CB][read_taxid] += 1
                    else:
                        nested_dict[bread_CB][read_taxid] = 1                    
                # If CB doesn't exist, create list with taxonomy ID
                else:
                    nested_dict[bread_CB] = {read_taxid: 1}

                use_count += 1
                
        # Make sparse matrix
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(nested_dict, taxid2node))        
    if mode == "RG":
        # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
        nested_dict = {}
        with pysam.AlignmentFile(args.cb_bam, "rb",check_sq = False) as cb_file:

            for bread in cb_file:
                total_count += 1
                # Check if the read exists in the kraken file
                if bread.query_name not in read_taxid_info_dict:
                    skipped += 1
                    # logging.warning("Read name {} not found in kraken file".format(sread.query_name))
                    continue    
                try:
                    bread_RG = bread.get_tag("RG")
                except:
                    skipped += 1
                    continue
                
                read_taxid = read_taxid_info_dict[bread.query_name]["taxid"]
                # add barcode information to read_taxid_info_dict
                read_taxid_info_dict[bread.query_name]["Barcode"] = bread_RG
                
                # Make nested dictionary with RG and taxonomy IDs
                # {cellbarcode: {taxonomyID}
                # If RG exists, add taxonomy ID to list 
                if bread_RG in nested_dict:
                    if read_taxid in nested_dict[bread_RG]:
                        nested_dict[bread_RG][read_taxid] += 1
                    else:
                        nested_dict[bread_RG][read_taxid] = 1                    
                # If RG doesn't exist, create list with taxonomy ID
                else:
                    nested_dict[bread_RG] = {read_taxid: 1}

                use_count += 1

        # Make sparse matrix
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(nested_dict, taxid2node))
    logger.info(f'Parsing bam file complete, total reads: {total_count}, use reads: {use_count}, skipped reads: {skipped}', status='complete')


    sparsematrix =  csr_matrix((vals, (rows, cols)))
    # Get mpa name for taxonomy ID
    taxname_list = [taxid2node[k].name for k in taxid_list]
    # store sparse matrix
    mmwrite(args.matrixfile, sparsematrix)
    taxa_df = pd.DataFrame(data=csr_matrix.todense(sparsematrix))
    taxa_df.index = taxname_list
    taxa_df.columns = cell_list

    taxa_df.to_csv(args.profile_tsv,sep=",")

    # Store list of cell barcodes
    with open(args.cellfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\n')
        tsv_output.writerow(cell_list)
    
    # Store list of taxonomy IDs
    data = zip(taxid_list, taxname_list)
    with open(args.taxfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for idx, tax in data:
            tsv_output.writerow([idx, tax])

    # Output detailed results with barcode information
    if args.output_read_tsv:
        logger.info(f'Building taxonomic lineage information for detailed output', status='run')
        
        # Build a dictionary mapping each taxonomy ID to its full taxonomic path
        taxid_to_lineage = {}
        taxonomic_levels = {'D': 'Domain', 'P': 'Phylum', 'C': 'Class', 'O': 'Order', 
                           'F': 'Family', 'G': 'Genus', 'S': 'Species'}
        
        # Process all taxonomy IDs found in read assignments
        unique_taxids = set(taxid_info_dict.keys())
        for taxid in unique_taxids:
            if taxid in taxid2node:
                node = taxid2node[taxid]
                # Get taxonomic lineage path
                lineage = {}
                
                # Start with the current node
                lineage[node.level_rank] = node.name
                
                # Traverse up to get all ancestors
                current = node
                while current.parent is not None:
                    current = current.parent
                    if current.level_rank in taxonomic_levels:
                        lineage[current.level_rank] = current.name
                
                taxid_to_lineage[taxid] = lineage
        
        # Convert dictionary to list for DataFrame creation
        detailed_results = []
        for read_name, info in read_taxid_info_dict.items():
            if "Barcode" in info:  # Only output reads with barcode information
                taxid = info["taxid"]
                lineage = taxid_to_lineage[taxid]
                result = {
                    "Read_Name": read_name,
                    "Taxonomy_ID": info["taxid"],
                    "Barcode": info["Barcode"],
                    "Species": lineage.get("S", np.nan),
                    "Genus": lineage.get("G", np.nan),
                    "Family": lineage.get("F", np.nan),
                    "Order": lineage.get("O", np.nan),
                    "Class": lineage.get("C", np.nan),
                    "Phylum": lineage.get("P", np.nan),
                    "Domain": lineage.get("D", np.nan)
                }
                
                if mode == "CB_UMI" and "UMI" in info:
                    result["UMI"] = info["UMI"]
                
                detailed_results.append(result)
        
        # Create DataFrame and output as TSV
        detailed_df = pd.DataFrame(detailed_results)
        detailed_df.to_csv(args.output_read_tsv, sep="\t", index=False)
        logger.info(f"Detailed TSV with taxonomic lineage information saved to {args.output_read_tsv}", status="complete")

    logger.info(f'Finish Saving the result', status='Complete')


if __name__ == "__main__":
    main()