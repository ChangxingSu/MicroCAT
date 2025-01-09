import pysam
import re
import collections.abc
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import csv
import logging
import os
import argparse
import pandas as pd
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

        while not level_rank == 'D':
            parent = child.parent
            level_rank = parent.level_rank
            parent_taxid = parent.taxid
            name = parent.name
            name = name.replace(' ', '_')
            if level_rank in main_lvls:
                try:
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

def twist_dict(nested):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = {}
    for ckey, tdict in nested.items():
        for tkey, kvalue in tdict.items():
            if ckey in newdict:
                if kvalue in newdict[ckey]:
                    newdict[ckey][kvalue] += 1
                else:
                    newdict[ckey][kvalue] = 1
            else:
                newdict[ckey] = {kvalue: 1}
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


def twist_dict(nested_dict):
    """Make count dictionary with {cellbarcode: {taxonomyID: count}}"""
    
    newdict = {}
    
    for ckey, tdict in nested_dict.items():
        for tkey in tdict:  # Assuming tdict is a list of taxonomy IDs
            if ckey in newdict:
                if tkey in newdict[ckey]:
                    newdict[ckey][tkey] += 1
                else:
                    newdict[ckey][tkey] = 1
            else:
                newdict[ckey] = {tkey: 1}

    # Process genus-species relationships
    for ckey, tdict in newdict.items():
        for kvalue, count in tdict.items():
            # Check if kvalue has a genus in the same cellbarcode
            rank = taxid2node[str(kvalue)].level_rank
            genus_taxid = taxid2node[str(kvalue)].taxid_to_desired_rank("G")
            if rank == 'S' and genus_taxid in newdict[ckey] and genus_taxid != "error - taxid above desired rank, or not annotated at desired rank":
                newdict[ckey][genus_taxid] += count

    return newdict


def dict2lists(nested):
    """ Returns lists for sparse matrix """
    rows = []  # cell coordinate
    columns = []  # taxonomy id coordinate
    values = []  # count

    cell_list = list(nested.keys())  # cell barcodes in the same order as rows
    taxid_list = []  # list to store unique tax IDs

    for ckey, taxdict in nested.items():
        for taxkey, count in taxdict.items():
            if taxkey not in taxid_list:
                taxid_list.append(taxkey)
            
            k = taxid_list.index(taxkey)
            rows.append(k)
            columns.append(cell_list.index(ckey))
            values.append(count)

    return rows, columns, values, cell_list, taxid_list



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='This script is used to output kraken2 classified microbial data in cellranger format as feature.tsv,barcodes.tsv,matrix.mtx \n This requires the additional packages pysam(If your python version is up to 3.9)\n')
    parser.add_argument('--kraken_output',help='Kraken output file.')
    parser.add_argument('--bam', dest='bam', 
                        help="The bam file after human host comparison, as input to kraken")
    parser.add_argument('--outdir', dest='outdir', default='krak2sc', 
                        help="name of the folder to download the genomes to. If this already exists, the result will be added to it. By default this is krak2sc")
    parser.add_argument('--log_file', dest='log_file', default='logfile_krak2sc.log',
                        help="File to write the log to")
    parser.add_argument("--krak_study_denosing_file", action="store", help="path to krak_study_denosing file")
    parser.add_argument('--ktaxonomy', required=True,
        help='Kraken2 database ktaxonomy file path')
    parser.add_argument('--tsv_output', required=True,
        help='Microbiome profile tsv output')
    parser.add_argument('--verbose', action='store_true', help='Detailed print')

    args = parser.parse_args()
    bamfile = args.bam
    outdir = args.outdir
    log_file = args.log_file
    ktaxonomy = args.ktaxonomy
    kraken_output = args.kraken_output
    tsv_output = args.tsv_output
    krak_study_denosing_file = args.krak_study_denosing_file
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
        taxid2node = make_dicts(ktaxonomy)
        logger.info('Successfully parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='complete')
    except (ValueError, KeyError) as e:
        logger.error(f"An error occurred while processing the Kraken taxonomy file: {e}")
        sys.exit()

    logger.info('Reading kraken2 classifier result infomation from denosed report', status='run')
    # Read taxa file (krak_study_denosing)
    krak_study_denosing = pd.read_csv(krak_study_denosing_file, sep="\t")

    # krak2_output_copy['taxid'] =krak2_output_copy['taxid'].astype(str)
    krak_study_denosing['ncbi_taxa'] = krak_study_denosing['ncbi_taxa'].astype(str)
    desired_taxid_list = set(krak_study_denosing['ncbi_taxa'].unique())
    logger.info('Finished processing kraken2 classifier result', status='complete')

    if not os.path.exists(outdir):
        os.system('mkdir -p '+outdir)

    all_desired_taxid_list = []  # List to store descendants_taxid_list for each taxid
    taxinfo = {}  # Dictionary to store all descendants and their corresponding taxids

    for taxid in desired_taxid_list:
        descendants_taxid_list = []
        rank = taxid2node[str(taxid)].level_rank

        # Retrieve descendants based on rank
        if rank == "S":
            descendants_taxid_list = taxid2node[str(taxid)].get_all_descendants()
        if rank == "G":
            descendants_taxid_list.append(taxid)
            descendants_nodes_list = taxid2node[str(taxid)].children
            while len(descendants_nodes_list) > 0:
                # For this node
                curr_n = descendants_nodes_list.pop()
                curr_rank = curr_n.level_rank
                if curr_rank != "S":
                    descendants_taxid_list.append(curr_n.taxid)
        
        # Extend the all_desired_taxid_list with elements of descendants_taxid_list
        all_desired_taxid_list.extend(descendants_taxid_list)

        # If descendants_taxid_list is not empty, create a dictionary
        if descendants_taxid_list:
            new_descendants_dict = {descendant: str(taxid) for descendant in descendants_taxid_list}
            # Update taxinfo dictionary with new key-value pairs
            taxinfo.update(new_descendants_dict)

    # Generate variables based on input
    matrixfile = os.path.join(outdir,'matrix.mtx')
    cellfile = os.path.join(outdir, 'barcodes.tsv')
    taxfile = os.path.join(outdir,'features.tsv')

    # Extract taxonomy IDs for each transcript
    logger.info('Extracting taxonomy IDs for each transcript', status='run')

    line = 0
    skipped = 0
    use_count = 0
    # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
    nested_dict = {}

    # Load the kraken file into memory
    kraken_data = {}
    with open(kraken_output, "r") as krakenfile:
        for kread in krakenfile:
            read_type, query_name, taxid_info, read_len, kmer_position = kread.strip().split('\t')
            if read_type == "U":
                continue
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                kread_taxid = re.search(r'\(taxid (\d+)\)', taxid_info).group(1)                
            except (ValueError, KeyError) as e:
                # in this case, something is wrong!
                logger.error(f"An error occurred while processing the Kraken output file: {e}")
                logger.error(f"Here is an error. Queryname: {query_name}")
                continue
            if kread_taxid in all_desired_taxid_list:
                # Store as main level id  
                main_level_id = taxinfo[kread_taxid]
                kraken_data[query_name] = main_level_id

    # Iterate through the bam file
    for sread in pysam.AlignmentFile(bamfile, "rb"):
        # count the total number of reads analysed
        line += 1

        # Check if the read exists in the kraken file
        if sread.query_name not in kraken_data:
            skipped += 1
            # logging.warning("Read name {} not found in kraken file".format(sread.query_name))
            continue

        # Use the kraken data for this read
        kread_taxid = kraken_data[sread.query_name]

        # Get cell barcode and UMI from bam file
        try:
            sread_CB = sread.get_tag('RG')
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

        # Make nested dictionary with RG and taxonomy IDs
        if sread_CB in nested_dict:
            # If RG exists, add taxonomy ID to list 
            nested_dict[sread_CB].append(kread_taxid)
        else:
            # If RG doesn't exist, create list with taxonomy ID
            nested_dict[sread_CB] = [kread_taxid]
        
        use_count += 1

    logger.info(f'Finished extracting taxonomy IDs for each transcript', status='complete')
    logger.info(f'Total extracted reads: {line}', status='summary')
    logger.info(f'Total classified Reads : {use_count}', status='summary')
    logger.info(f'Skipped reads: {skipped}', status='summary')

    logger.info(f'Saving the result', status='run')

    # Make sparse matrix
    rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(nested_dict))
    sparsematrix = csr_matrix((vals, (rows, cols)))

    # # Get mpa name for taxonomy ID
    taxname_list = [taxid2node[k].name for k in taxid_list]
    # store sparse matrix
    mmwrite(matrixfile, sparsematrix)
    taxa_df = pd.DataFrame(data=csr_matrix.todense(sparsematrix))
    taxa_df.index = taxname_list
    taxa_df.columns = cell_list

    taxa_df.to_csv(tsv_output,sep="\t")
    
    # Store list of cell barcodes
    with open(cellfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\n')
        tsv_output.writerow(cell_list)
    
    # Store list of taxonomy IDs
    data = zip(taxid_list, taxname_list)
    with open(taxfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for idx, tax  in data:
            tsv_output.writerow([idx, tax])

    logger.info(f'Finish Saving the result', status='Complete')