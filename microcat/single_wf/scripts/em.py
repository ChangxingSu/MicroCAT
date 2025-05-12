import os
from flatten_dict import unflatten
import pathlib
import csv
from operator import add, mul
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from operator import add
import math
import pysam
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import logging
import re

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
    """
    Custom logging function that includes an optional status message.

    Args:
        level (int): Logging level (e.g., logging.INFO).
        msg (str): The main log message.
        *args: Arguments for the log message formatting.
        status (str, optional): An optional status indicator (e.g., 'run', 'done', 'error'). Defaults to None.
    """
    if status:
        msg = f'({status}) {msg}'  # Concatenate the message and status
    logger.log(level, msg, *args)

# Bind the custom_log function to the logger object for different log levels
# This allows calling logger.info(..., status='run') etc.
logger.info = lambda msg, *args, status=None: custom_log(logging.INFO, msg, *args, status=status)
logger.warning = lambda msg, *args, status=None: custom_log(logging.WARNING, msg, *args, status=status)
logger.error = lambda msg, *args, status=None: custom_log(logging.ERROR, msg, *args, status=status)
logger.debug = lambda msg, *args, status=None: custom_log(logging.DEBUG, msg, *args, status=status)

def get_path_to_root(node):
    """
    Gets the full path from the current node to the root node in the taxonomy tree.
    获取从当前节点到根节点的完整路径

    Parameters:
    - node (Tree): The starting node.
                   起始节点

    Returns:
    - List[Tree]: A list of all nodes on the path from the current node to the root.

    """
    path = []
    current = node
    while current is not None:
        path.append(current)
        current = current.parent
    return path
    
def calculate_taxonomic_penalty(taxid2node, taxid1, taxid2, lambda_param=1.0, max_distance=7):
    """
    Calculates the taxonomic distance between two taxa and converts it into a penalty factor.
    The penalty is based on the distance to their lowest common ancestor (LCA) in the taxonomy tree.
    计算两个分类单元之间的分类学距离，并转换为惩罚因子

    Parameters:
    - taxid2node (dict): A dictionary mapping taxid to Tree node objects.
                         taxid到Tree节点的映射字典
    - taxid1 (str): The taxid of the first taxon.
                    第一个分类单元的taxid
    - taxid2 (str): The taxid of the second taxon.
                    第二个分类单元的taxid
    - lambda_param (float): Penalty strength parameter. Higher values mean stronger penalties for larger distances.
                            惩罚强度参数
    - max_distance (int): The maximum distance used for normalization. Distances beyond this are capped.
                          最大距离，用于归一化

    Returns:
    - float: The penalty factor (a value between 0 and 1). 1.0 means no penalty (identical taxids or close relatives),
             values closer to 0 mean higher penalty (distant relatives or missing taxids).
             惩罚因子 (0-1之间的值)
    """
    # If taxids are the same, no penalty
    # 如果taxid相同，无需惩罚
    if taxid1 == taxid2:
        return 1.0
        
    # Handle cases where taxid is "0" (unclassified) - apply maximum penalty
    # 处理taxid为0的情况
    if taxid1 == "0" or taxid2 == "0":
        return math.exp(-lambda_param)  # Maximum penalty
    
    # Handle cases where taxid is not found in the tree - apply maximum penalty
    # 处理taxid不在树中的情况
    if taxid1 not in taxid2node or taxid2 not in taxid2node:
        return math.exp(-lambda_param)  # Maximum penalty
    
    # Get the corresponding Tree nodes
    # 获取对应的Tree节点
    node1 = taxid2node[taxid1]
    node2 = taxid2node[taxid2]
    
    # Get paths to the root node
    # 获取到根节点的路径
    path1 = get_path_to_root(node1)
    path2 = get_path_to_root(node2)
    
    # Find the lowest common ancestor (LCA)
    # 找到最近共同祖先
    lca = None
    for n1 in path1:
        for n2 in path2:
            if n1.taxid == n2.taxid:  # Compare by taxid, not object identity
                lca = n1
                break
        if lca:
            break
            
    if lca is None:
        return math.exp(-lambda_param)  # No common ancestor
        
    # Calculate distances to the LCA
    # 计算距离
    dist1 = path1.index(lca)  # Distance from node1 to LCA
    dist2 = path2.index(lca)  # Distance from node2 to LCA
    distance = dist1 + dist2
    
    # Normalize the distance and calculate the penalty factor using an exponential decay function
    # 归一化距离并计算惩罚因子
    normalized_distance = min(distance / max_distance, 1.0)
    penalty = math.exp(-lambda_param * normalized_distance)
    
    return penalty

def precompute_taxonomic_penalties(read_kraken_taxids, log_p_rgs, taxid2node, lambda_param=1.0):
    """
    Precomputes taxonomic penalty factors for all relevant read-candidate species pairs.
    This avoids redundant calculations within the EM loop.
    预计算所有读段-物种对的分类学惩罚因子

    Args:
        read_kraken_taxids (dict): Mapping from read ID to its Kraken-assigned taxid {read_id: taxid}.
                                   读段到Kraken分类ID的映射 {read_id: taxid}
        log_p_rgs (dict): Dictionary of log likelihoods for reads mapped to candidate species {read_id: ([taxids], [scores])}.
                          读段对应的物种对数似然字典 {read_id: ([taxids], [scores])}
        taxid2node (dict): Mapping from taxid to Tree node objects.
                           分类ID到分类树节点的映射字典
        lambda_param (float): Penalty strength parameter for calculate_taxonomic_penalty.
                              惩罚强度参数

    Returns:
        dict: A dictionary storing precomputed penalties {(read_id, candidate_taxid): penalty}.
              预计算的惩罚因子字典 {(read_id, candidate_taxid): penalty}
    """
    penalty_dict = {}
    logger.info(f"Starting precomputation of taxonomic penalties", status="run")
    
    # To reduce redundant calculations, first compute penalties for all unique taxid pairs
    # 为了减少重复计算，先计算所有taxid对之间的惩罚
    taxid_pair_penalties = {}
    
    # Collect all unique (kraken_taxid, candidate_taxid) pairs
    # 收集所有唯一的taxid对
    unique_pairs = set()
    for r, (taxids, _) in log_p_rgs.items():
        if r in read_kraken_taxids:
            kraken_taxid = str(read_kraken_taxids[r]) # Ensure string comparison
            for candidate_taxid in taxids:
                # Ensure candidate_taxid is also a string for consistent pairing
                unique_pairs.add((kraken_taxid, str(candidate_taxid)))
    
    # Precompute penalties for all unique pairs
    # 预计算所有唯一taxid对的惩罚
    logger.info(f"Calculating penalties for {len(unique_pairs)} unique taxid pairs", status="run")
    for kraken_taxid, candidate_taxid in unique_pairs:
        taxid_pair_penalties[(kraken_taxid, candidate_taxid)] = calculate_taxonomic_penalty(
            taxid2node, kraken_taxid, candidate_taxid, lambda_param)
    logger.info(f"Finished calculating unique pair penalties", status="done")

    # Assign the precomputed penalties to each read and its candidate species
    # 为每个读段和其候选物种计算惩罚
    for r, (taxids, _) in log_p_rgs.items():
        if r in read_kraken_taxids:
            kraken_taxid = str(read_kraken_taxids[r])
            for candidate_taxid in taxids:
                candidate_taxid_str = str(candidate_taxid)
                pair_key = (kraken_taxid, candidate_taxid_str)
                # Use the precomputed penalty; default to 1.0 if pair somehow wasn't precomputed (shouldn't happen)
                penalty = taxid_pair_penalties.get(pair_key, 1.0)
                # Store penalty using the original candidate_taxid type from log_p_rgs (likely int)
                penalty_dict[(r, candidate_taxid)] = penalty
    
    logger.info(f"Precomputation complete. Stored {len(penalty_dict)} read-species penalties.", status="done")
    return penalty_dict


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
        return lineage # Return the lineage collected so far if root is reached

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

def get_align_stats(alignment):
    """Retrieve list of inquired cigar stats (I,D,S,X) for alignment

        alignment (pysam.AlignmentFile): align of interest
        return (list(int)): list of counts for each cigar operation defined in (I,D,S,X)
    """
    cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    return [cigar_stats[1], cigar_stats[2], cigar_stats[4], n_mismatch]

def get_cigar_op_log_probabilities(sam_path):
    """P(align_type) for each type in CIGAR_OPS by counting how often the corresponding
            operations occur in the primary alignments and by normalizing over the total
            sum of operations.

        sam_path(str): path to sam file of interest
        return: log probabilities (list(float)) for each cigar operation defined in CIGAR_OPS,
                where p > 0
            zero_locs (list(int)): list of indices (int) where probability == 0
            dict_longest_align (dict[str]:(int)): dict of max alignment length
                for each query read
    """
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    dict_longest_align = {}
    # pylint: disable=maybe-no-member
    sam_pysam = pysam.AlignmentFile(sam_path)
    # add alignment lengths and adjust existing alignment lengths in dict if necessary
    for alignment in sam_pysam.fetch(until_eof=True):
        align_len = get_align_len(alignment)
        if align_len not in dict_longest_align:
            dict_longest_align[alignment.query_name] = align_len
        if not alignment.is_secondary and not alignment.is_supplementary \
                and alignment.reference_name:
            cigar_stats_primary = list(map(add, cigar_stats_primary, get_align_stats(alignment)))
            # calculate cigar stats for alignment
            if dict_longest_align[alignment.query_name] < align_len:
                dict_longest_align[alignment.query_name] = align_len
    # check if any probabilities are 0, if so, remove
    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]
    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]
    n_char = sum(cigar_stats_primary)
    return [math.log(x) for x in np.array(cigar_stats_primary)/n_char], zero_locs, \
           dict_longest_align

def get_align_len(alignment):
    """Retrieve number of columns in alignment

        alignment (pysam.AlignmentFile): align of interest
        return (int): number of columns in alignment
    """
    return sum(alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS_ALL)

def compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len):
    """ 
    log(L(r|s)) = log(P(cigar_op)) × n_cigar_op for CIGAR_OPS

    Args:
        alignment(pysam.AlignmentFile): pysam alignment to score
        cigar_stats(list(int)): list of cigar stats to compute
        log_p_cigar_op(list(float)): list of cigar_op probabilities corresponding to cigar_stats;
                                    computed from primary alignments
        dict_longest_align (dict[str]:(int)): dict of max alignment length for each query read
        align_len (int): number of columns in the alignment
    Returns:
        log_score (float): log(L(r|s))
            query_name (str): query name in alignment
            species_tid (int): species-level taxonomy id corresponding to ref_name
    """

    ref_name, query_name = alignment.reference_name, alignment.query_name
    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats))) * \
                (dict_longest_align[query_name]/align_len)
    species_tid = int(seqid2tax_map[ref_name])
    return log_score, query_name, species_tid

def log_prob_rgs_dict(sam_path, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    """dict containing log(L(read|seq)) for all pairwise alignments in sam file

        sam_path(str): path to sam file
        log_p_cigar_op(list(float)): probability for each cigar operation defined in CIGAR_OPS,
                                         where p > 0
        dict_longest_align (dict[str]:(int)): dict of max alignment length for each query read
        zero_locs(list(int)): list of indices (int) where probability == 0
        return ({[str,int]:float}): dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
            int: unassigned read count
            int: assigned read count
    """
    # calculate log(L(read|seq)) for all alignments
    log_p_rgs, unassigned_set = {}, set()
    # pylint: disable=maybe-no-member
    sam_filename = pysam.AlignmentFile(sam_path, 'rb')

    if not p_cigar_op_zero_locs:
        for alignment in sam_filename.fetch(until_eof=True):
            align_len = get_align_len(alignment)
            identity = (alignment.query_alignment_length - alignment.get_tag("NM")) / dict_longest_align[alignment.query_name]
            if alignment.reference_name and align_len and identity >= 0.5:
                cigar_stats = get_align_stats(alignment)
                log_score, query_name, species_tid = \
                    compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op,
                                        dict_longest_align, align_len)

                if query_name not in log_p_rgs:
                    log_p_rgs[query_name] = ([species_tid], [log_score])
                elif query_name in log_p_rgs:
                    if species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] + [species_tid],
                                                 log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score

            else:
                unassigned_set.add(alignment.query_name)
    else:
        for alignment in sam_filename.fetch(until_eof=True):
            align_len = get_align_len(alignment)
            identity = (alignment.query_alignment_length - alignment.get_tag("NM")) / dict_longest_align[alignment.query_name]
            if alignment.reference_name and align_len and identity >= 0.5:
                cigar_stats = get_align_stats(alignment)
                if sum(cigar_stats[x] for x in p_cigar_op_zero_locs) == 0:
                    for i in sorted(p_cigar_op_zero_locs, reverse=True):
                        del cigar_stats[i]
                    log_score, query_name, species_tid = \
                        compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op,
                                            dict_longest_align, align_len)

                    if query_name not in log_p_rgs:
                        log_p_rgs[query_name] = ([species_tid], [log_score])
                    elif query_name in log_p_rgs and species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] +[species_tid],
                                                 log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        # keep the highest log_score
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score
            else:
                unassigned_set.add(alignment.query_name)

    assigned_reads = set(log_p_rgs.keys())
    unassigned_reads = unassigned_set - assigned_reads
    unassigned_count = len(unassigned_reads)
    print(f"Unassigned read count: {unassigned_count}\n")

    ## remove low likelihood alignments?
    ## remove if p(r|s) < 0.01
    #min_p_thresh = math.log(0.01)
    #log_p_rgs = {r_map: val for r_map, val in log_p_rgs.items() if val > min_p_thresh}
    return log_p_rgs, unassigned_count, len(assigned_reads)

def parse_kraken_qc(kraken_qc_file):
    """
    解析 Kraken QC 文件，提取物种 ID 和相对丰度
    
    Args:
        kraken_qc_file (str): Kraken QC 文件路径
    
    Returns:
        dict: {物种ID (int): 相对丰度 (float)}
    """
    
    try:
        krak2_qc = pd.read_csv(kraken_qc_file, sep='\t')
        krak2_qc_species = krak2_qc[krak2_qc['classification_rank'] == 'S']
        total_reads = krak2_qc_species['fragments'].sum()
        # divide by total reads to get relative abundance
        krak2_qc_species['relative_abundance'] = krak2_qc_species['fragments'] / total_reads
        krak2_qc_species['ncbi_taxa'] = krak2_qc_species['ncbi_taxa'].astype(int)
        species_abundance = dict(zip(krak2_qc_species['ncbi_taxa'], krak2_qc_species['relative_abundance']))
    except FileNotFoundError:
        logger.error(f"找不到 Kraken QC 文件: {kraken_qc_file}", status="error")
        sys.exit()
    except Exception as e:
        logger.error(f"解析 Kraken QC 文件时出错: {e}", status="error")
        sys.exit()

    return species_abundance

def expectation_maximization_iterations(log_p_rgs, freq, lli_thresh, input_threshold, read_kraken_taxids=None, taxid2node=None, lambda_param=1.0, save_iterations=None):
    """
    执行EM算法迭代，包含迭代质量保持阈值步骤以解决长尾效应，并加入Kraken2分类学惩罚
    
    Args:
        log_p_rgs: 读段对应的物种对数似然字典 {read_id: ([taxids], [scores])}
        freq: 初始物种丰度估计
        lli_thresh: 对数似然增量阈值
        input_threshold: 物种丰度阈值
        read_kraken_taxids: 读段到Kraken分类ID的映射字典 {read_id: taxid}
        taxid2node: 分类ID到分类树节点的映射字典
        lambda_param: 惩罚强度参数
        save_iterations: 需要保存结果的迭代次数列表，例如[1,5,10,15,20]，默认为None
    
    Returns:
        freq_full: 所有物种的丰度
        freq_set_thresh: 经过阈值处理后的物种丰度
        p_sgr: 读段-物种分配概率字典
        iteration_results: 保存的迭代结果字典 {迭代次数: {物种ID: 丰度}}
    """
    n_reads = len(log_p_rgs)
    if n_reads == 0:
        return freq, freq, {}, {}
    
    # 初始化用于保存迭代结果的字典
    iteration_results = {}
    
    # 如果有Kraken分类信息和分类树信息，预计算惩罚因子
    penalty_dict = {}
    if read_kraken_taxids and taxid2node:
        logger.info("预计算分类学惩罚因子", status="run")
        penalty_dict = precompute_taxonomic_penalties(
            read_kraken_taxids, log_p_rgs, taxid2node, lambda_param)
        logger.info(f"预计算完成，共{len(penalty_dict)}个惩罚因子", status="done")
    
    # 初始化物种有效性标记
    strain_valid = {strain: True for strain in freq}
    
    freq_thresh = 1/n_reads
    if n_reads > 1000:
        freq_thresh = 10/n_reads
    
    # 初始化
    freq_full = freq.copy()
    prev_log_likelihood = -float('inf')
    
    # 每10次EM迭代执行一次阈值处理
    thresholding_iter_step = 10
    can_help = True
    counter = 0
    
    # 确保初始频率值都为正数
    for g in freq:
        if freq[g] <= 0:
            freq[g] = 1e-10  # 设置一个非常小但非零的值
    
    while True:
        counter += 1
        
        # 判断是否执行阈值处理
        if counter % thresholding_iter_step == 0 and can_help:
            strain_valid, potentially_removable, can_help = apply_set_cover(
                log_p_rgs, freq, strain_valid, max(freq_thresh, input_threshold))
        
        # 执行常规EM步骤，但仅考虑有效物种
        strain_read_count = {strain: 0 for strain in freq}
        p_sgr = {}  # p(s|r)
        
        # Expectation步骤
        log_likelihood = 0
        for r, (taxids, log_scores) in log_p_rgs.items():
            # 过滤掉无效物种和频率为零的物种
            valid_indices = []
            for i, g in enumerate(taxids):
                if strain_valid[g] and freq[g] > 0:
                    valid_indices.append(i)
            
            # 如果该读段没有有效的映射物种，则跳过
            if not valid_indices:
                continue
                
            # 计算读段r归属于各个物种的概率
            p_r = {}  # 归一化概率p(s|r)
            max_log_p = -float('inf')
            
            # 应用分类学惩罚并找出最大对数概率
            for i in valid_indices:
                g = taxids[i]
                log_p = log_scores[i]
                
                # 如果有预计算的惩罚因子，应用它
                penalty = 1.0
                if (r, g) in penalty_dict:
                    penalty = penalty_dict[(r, g)]
                
                try:
                    # 加上log(penalty)到对数似然
                    log_penalty = math.log(max(penalty, 1e-300))
                    log_value = log_p + math.log(max(freq[g], 1e-300)) + log_penalty
                    # print(f"log_value: {log_value}, log_p: {log_p}, freq[g]: {math.log(max(freq[g], 1e-300))}, penalty: {log_penalty}")
                    if log_value > max_log_p:
                        max_log_p = log_value
                except ValueError:
                    # 如果出现数值错误，记录并跳过
                    print(f"Warning: math domain error for g={g}, freq[g]={freq[g]}")
                    continue
            
            # 如果无法计算最大对数概率，则跳过
            if max_log_p == -float('inf'):
                continue
                
            # 计算分母(归一化因子)
            log_denom = 0
            first = True
            for i in valid_indices:
                g = taxids[i]
                log_p = log_scores[i]
                
                # 应用分类学惩罚
                penalty = 1.0
                if (r, g) in penalty_dict:
                    penalty = penalty_dict[(r, g)]
                
                try:
                    log_penalty = math.log(max(penalty, 1e-300))
                    log_num = log_p + math.log(max(freq[g], 1e-300)) + log_penalty - max_log_p
                    if first:
                        log_denom = log_num
                        first = False
                    else:
                        log_denom = np.logaddexp(log_denom, log_num)
                except ValueError:
                    continue
            
            # 计算每个有效物种的后验概率p(s|r)
            p_r = {}
            for i in valid_indices:
                g = taxids[i]
                log_p = log_scores[i]
                
                # 应用分类学惩罚
                penalty = 1.0
                if (r, g) in penalty_dict:
                    penalty = penalty_dict[(r, g)]
                
                try:
                    log_penalty = math.log(max(penalty, 1e-300))
                    log_num = log_p + math.log(max(freq[g], 1e-300)) + log_penalty - max_log_p
                    p_r[g] = math.exp(log_num - log_denom)
                    strain_read_count[g] += p_r[g]
                except (ValueError, OverflowError):
                    continue
            
            if p_r:  # 只有当p_r非空时才添加
                p_sgr[r] = p_r
                # 更新对数似然
                log_likelihood += max_log_p + log_denom
        
        # Maximization步骤 - 更新丰度
        total_reads_assigned = sum(strain_read_count.values())
        if total_reads_assigned > 0:  # 确保分母非零
            for g in freq:
                if strain_valid[g]:
                    freq[g] = max(strain_read_count[g] / total_reads_assigned, 1e-10)
                else:
                    freq[g] = 0
        
        # 保存当前迭代的结果（如果需要）
        if save_iterations is not None and counter in save_iterations:
            # 创建一个副本以避免引用问题
            current_freq = freq.copy()
            # 对频率进行归一化处理
            sum_freq = sum(v for v in current_freq.values())
            if sum_freq > 0:
                current_freq = {k: v/sum_freq for k, v in current_freq.items()}
            # 保存当前迭代结果
            iteration_results[counter] = current_freq
            logger.info(f"保存第{counter}次迭代结果", status="done")
        
        # 检查收敛性
        if abs(log_likelihood - prev_log_likelihood) < lli_thresh:
            # 如果最后一次迭代不在save_iterations中，保存它
            if save_iterations is not None and counter not in save_iterations:
                current_freq = freq.copy()
                sum_freq = sum(v for v in current_freq.values())
                if sum_freq > 0:
                    current_freq = {k: v/sum_freq for k, v in current_freq.items()}
                iteration_results["final"] = current_freq
                logger.info(f"保存最终迭代结果（第{counter}次）", status="done")
            break
        
        prev_log_likelihood = log_likelihood
    
    # 最终过滤结果
    freq_full = freq.copy()
    freq_set_thresh = {k: v for k, v in freq.items() if strain_valid[k]}
    
    # 归一化最终结果
    total = sum(freq_set_thresh.values())
    if total > 0:
        freq_set_thresh = {k: v/total for k, v in freq_set_thresh.items()}
    
    return freq_full, freq_set_thresh, p_sgr, iteration_results


def apply_set_cover(log_p_rgs, freq, strain_valid, min_count):
    """
    执行集合覆盖算法，寻找最小必要的物种集合以覆盖所有的读段
    
    Args:
        log_p_rgs: 读段对应的物种对数似然字典 {read_id: ([taxids], [scores])}
        freq: 当前物种丰度估计
        strain_valid: 物种有效性标记字典
        min_count: 最小丰度阈值
    
    Returns:
        strain_valid: 更新后的物种有效性标记
        strain_potentially_removable: 潜在可移除物种标记
        can_help: 是否可以继续执行阈值处理
    """
    previously_valid = sum(1 for s, valid in strain_valid.items() if valid)
    
    # 1. 标记潜在可移除物种 (PR)
    strain_potentially_removable = {s: (strain_valid[s] and freq[s] <= min_count) 
                                   for s in strain_valid}
    
    # 2. 找出拥有唯一"有效"读段的物种
    # 这些物种的读段无法被其他物种解释，因此不能被移除
    unique_strain_reads = {}
    for r, (taxids, _) in log_p_rgs.items():
        valid_strains = [s for s in taxids if strain_valid[s]]
        if len(valid_strains) == 1:
            strain = valid_strains[0]
            if strain not in unique_strain_reads:
                unique_strain_reads[strain] = set()
            unique_strain_reads[strain].add(r)
    
    # 具有唯一读段的物种不能被移除
    for strain in unique_strain_reads:
        strain_potentially_removable[strain] = False
    
    # 3. 识别关键读段（所有可比对物种都是PR的读段）
    critical_reads = {}
    for r, (taxids, _) in log_p_rgs.items():
        valid_strains = [s for s in taxids if strain_valid[s]]
        if all(strain_potentially_removable[s] for s in valid_strains):
            critical_reads[r] = set(valid_strains)
    
    # 如果没有关键读段，可以直接移除所有PR物种
    if len(critical_reads) == 0:
        for strain in strain_valid:
            if strain_potentially_removable[strain]:
                strain_valid[strain] = False
    else:
        # 4. 创建物种到读段的映射（用于集合覆盖问题）
        strain_to_reads = {}
        for strain in freq:
            if strain_valid[strain] and strain_potentially_removable[strain]:
                strain_to_reads[strain] = set()
                for r in critical_reads:
                    if strain in log_p_rgs[r][0]:  # 检查物种是否在taxids列表中
                        strain_to_reads[strain].add(r)
        
        # 5. 求解集合覆盖问题（使用贪心算法）
        must_keep = set()
        remaining_reads = set(critical_reads.keys())
        
        while remaining_reads and strain_to_reads:
            # 找到覆盖最多剩余读段的物种
            best_strain = None
            best_coverage = 0
            best_score = 0
            
            for strain, reads in strain_to_reads.items():
                covered_reads = reads.intersection(remaining_reads)
                
                if len(covered_reads) == 0:
                    continue
                
                # 计算该物种的得分（可以考虑丰度及其覆盖度）
                coverage = len(covered_reads)
                score = freq[strain] * coverage
                
                if score > best_score:
                    best_score = score
                    best_coverage = coverage
                    best_strain = strain
            
            if best_strain is None:
                break
            
            # 将最佳物种添加到必须保留的集合
            must_keep.add(best_strain)
            
            # 更新剩余读段
            covered_reads = strain_to_reads[best_strain].intersection(remaining_reads)
            remaining_reads -= covered_reads
            
            # 移除已处理的物种
            del strain_to_reads[best_strain]
        
        # 6. 更新物种有效性
        for strain in strain_valid:
            if strain_potentially_removable[strain] and strain not in must_keep:
                strain_valid[strain] = False
    
    # 检查是否有变化
    current_valid = sum(1 for s, valid in strain_valid.items() if valid)
    can_help = previously_valid != current_valid
    
    return strain_valid, strain_potentially_removable, can_help


# static global variables
CIGAR_OPS = [1, 2, 4, 10]
CIGAR_OPS_ALL = [0, 1, 2, 4]


# Move all execution code inside this condition
if __name__ == "__main__":
    # 添加参数解析
    parser = argparse.ArgumentParser(description='执行EM算法进行物种丰度估计')
    parser.add_argument('--bam_file', required=True, help='BAM文件路径')
    parser.add_argument('--taxonomy_file', required=True, help='分类文件路径')
    parser.add_argument('--output', required=True, help='输出文件路径')
    parser.add_argument('--kraken_output', required=True, help='kraken输出文件路径')
    parser.add_argument('--ktaxonomy_file', required=True, help='Kraken2 ktaxonomy file')
    parser.add_argument('--log_file', default='emu_em.log', help='日志文件路径')
    parser.add_argument('--verbose', action='store_true', help='是否输出详细日志信息')
    parser.add_argument('--save_iterations', type=str, help='要保存的迭代结果，逗号分隔的数字，例如 "1,5,10,15,20"')
    parser.add_argument('--iteration_output', help='迭代结果输出目录')
    parser.add_argument('--kraken_qc', help='Kraken QC 文件路径，用于初始化物种丰度')
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


    bam_file = args.bam_file
    taxonomy_file = args.taxonomy_file
    output = args.output
    read_classifications_file = output
    # log_level = args.log_level

    # Load the BAM file
    logger.info('Loading BAM file', status='run')
    log_prob_cigar_op, locs_p_cigar_zero, longest_align_dict = get_cigar_op_log_probabilities(bam_file)

    logger.info('Loading log probabilities of CIGAR operations', status='done')
    logger.info('Loading zero locations', status='done')
    logger.info('Loading longest alignment dictionary', status='done')


    # Loading taxonomy database
    logger.info('Loading taxonomy database', status='run')
    # use int for taxid, to accelerate the lookup
    seqid2tax_df = pd.read_csv(taxonomy_file, sep="\t",dtype={'seqid': str, 'taxid': int})
    seqid2tax_map = dict(zip(seqid2tax_df['seqid'], seqid2tax_df['taxid']))
    logger.info(f'Loading taxonomy database with {len(seqid2tax_map)} reference sequences, {seqid2tax_df["taxid"].nunique()} taxa', status='done')

    # Calculate log likelihood of reads given sequences
    logger.info('Calculating log likelihood of reads given sequences', status='run')
    log_prob_rgs, counts_unassigned, counts_assigned = log_prob_rgs_dict(
        bam_file, log_prob_cigar_op, longest_align_dict, locs_p_cigar_zero)


    logger.info('Calculating log likelihood of reads given sequences', status='done')

    # Run the EM algorithm
    db_ids = list(seqid2tax_map.values())
    n_db = len(db_ids)
    n_reads = len(log_prob_rgs)
    print("Assigned read count: {}\n".format(n_reads))
    # check if there are enough reads


    counter = 1

    # 处理迭代保存参数
    if args.save_iterations:
        save_iterations = [int(x) for x in args.save_iterations.split(',')]
        logger.info(f"将保存以下迭代次数的结果: {save_iterations}")
    else:
        save_iterations = None

    # 在解析完BAM文件和分类数据库后
    # 加载Kraken分类结果
    logger.info('加载Kraken分类结果', status='run')
    read_kraken_taxids = {}

    # 假设我们从参数中获取kraken输出文件路径
    if args.kraken_output:
        with open(args.kraken_output, 'r') as kfile:
            for kraken_line in kfile:
                try:
                    read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                    if read_type == "C":  # 只处理已分类的读段
                        kread_taxid = re.search(r'\(taxid (\d+)\)', taxid_info).group(1)
                        read_kraken_taxids[query_name] = int(kread_taxid)
                except (ValueError, AttributeError) as e:
                    logger.error(f"处理Kraken输出时出错: {e}")
                    continue
        logger.info(f'成功加载{len(read_kraken_taxids)}个Kraken分类结果', status='done')
    else:
        logger.info('未提供Kraken输出文件，将不使用分类学惩罚', status='warning')

    # 构建分类树
    taxid2node = {}
    # 这里需要添加构建分类树的代码，类似于extract_specific_kraken_reads.py中的make_dicts函数
    logger.info('Parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='run')
    try:
        taxid2node = make_dicts(args.ktaxonomy_file)
        logger.info('Successfully parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='complete')
    except:
        logger.error("Couldn't get the taxonmy full lineage infomation from NCBI nodes.dump")
        sys.exit()

    # 使用 Kraken QC 文件初始化物种丰度，如果文件不存在则使用均匀分布
    if args.kraken_qc:
        logger.info('从 Kraken QC 文件加载初始物种丰度', status='run')
        kraken_abundances = parse_kraken_qc(args.kraken_qc)
        
        # 初始化所有物种的丰度为很小的值
        freq = dict.fromkeys(db_ids, 1e-5)
        # TODO 现在版本中，由于是全部样本的微生物进行建库比对，所以有的样本qc里会没有reference的微生物
        # 更新从 Kraken QC 文件中提取的物种丰度
        for taxid, abundance in kraken_abundances.items():
            if taxid in freq:
                freq[taxid] = max(abundance, 1e-5)  # 确保丰度不为0

        # 归一化丰度，确保总和为1
        total = sum(freq.values())
        freq = {k: v/total for k, v in freq.items()}
        
        # # 检测kraken_abundances是否包含所有db_ids
        # missing_taxa = set(db_ids) - set(kraken_abundances.keys())
        # if missing_taxa:
        #     logger.error(f"Kraken QC 文件中缺少以下物种: {missing_taxa}")
        #     sys.exit()

        
        logger.info(f'成功从 Kraken QC 文件加载初始物种丰度 ({len(kraken_abundances)} 个物种)', status='done')
    else:
        # 使用均匀分布
        freq = dict.fromkeys(db_ids, 1 / n_db)

    # 在调用EM算法时传入Kraken分类信息和分类树
    f_full, f_set_thresh, read_dist, iteration_results = expectation_maximization_iterations(
        log_prob_rgs, freq, 0.01, 0.00001, read_kraken_taxids, taxid2node, 2.5, save_iterations)

    print(f"Number of EM iterations: {counter}\n")
    print(f"Number of species in f_full: {len(f_full)}")
                                                                      
    # results_df = pd.DataFrame(zip(list(f_full.keys()) + ['unassigned'],
    #                                 list(f_full.values()) + [0]),
    #                             columns=["tax_id", "abundance"])

    # print(results_df)

    # results_df.to_csv(output, sep="\t", index=False)

    logger.info('EM算法完成,开始处理read分类', status='run')

    # 创建输出目录（如果不存在）
    # output_dir = os.path.dirname(args.output)
    # read_classifications_file = os.path.join(output_dir, "read_classifications.tsv")

    # 为每个read分配最可能的分类并直接写入文件
    logger.info(f"开始将read分类结果写入文件: {read_classifications_file}")

    # print(read_dist.keys())
    # with open(read_classifications_file, 'w', newline='') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     writer.writerow(["Read_Name", "Taxonomy_ID", "Probability"])
    
    #     for (read_name, tax_id), prob in read_dist.items():
    #         if read_name not in read_classifications or prob > read_classifications[read_name][1]:
    #             read_classifications[read_name] = (tax_id, prob)
    #             writer.writerow([read_name, tax_id, prob])
    # 为每个read选择最高概率的分类
    best_classifications = {}
    # 修复这段代码，正确解析read_dist结构
    for read_name, tax_probs in read_dist.items():
        best_tax_id = None
        best_prob = -1
        for tax_id, prob in tax_probs.items():
            if prob > best_prob:
                best_tax_id = tax_id
                best_prob = prob
        if best_tax_id is not None:
            best_classifications[read_name] = (best_tax_id, best_prob)


    # a very simple taxonomy name parser
    tax_id_to_name = {}
    with open(args.ktaxonomy_file, 'r') as kfile:
        for line in kfile:
            [taxid, p_tid, rank, lvl_num, name] = line.strip().split('\t|\t')
            tax_id_to_name[int(taxid)] = name

    # 将结果写入文件
    logger.info(f"开始将read分类结果写入文件: {read_classifications_file}")

    # TODO: may be dont need to write species and genus information, which will handle in the bam2mtx.py
    with open(output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Read_Name", "Taxonomy_ID", "Probability", "Species", "Genus"])
        
        for read_name, (tax_id, prob) in best_classifications.items():
            tax_name = tax_id_to_name.get(tax_id, "Unknown")
            tax_genus = tax_name.split(" ")[0] if tax_name != "Unknown" else "Unknown"
            writer.writerow([read_name, tax_id, prob, tax_name, tax_genus])


    # 如果指定了迭代输出目录，则保存每个迭代的结果
    if args.iteration_output and save_iterations:
        os.makedirs(args.iteration_output, exist_ok=True)
        for iter_num, iter_freq in iteration_results.items():
            iter_output_file = os.path.join(args.iteration_output, f"em_iteration_{iter_num}.tsv")
            # 将迭代结果转换为DataFrame并保存
            iter_df = pd.DataFrame([(tax_id, abundance) for tax_id, abundance in iter_freq.items() if abundance > 0],
                                    columns=["tax_id", "abundance"])
            
            # 添加物种名称
            iter_df["species"] = iter_df["tax_id"].apply(lambda x: tax_id_to_name.get(x, "Unknown"))
            iter_df.sort_values("abundance", ascending=False, inplace=True)
            
            iter_df.to_csv(iter_output_file, sep="\t", index=False)
            logger.info(f"保存迭代{iter_num}的结果到: {iter_output_file}", status="done")
    logger.info(f"Read分类结果已保存到文件: {read_classifications_file}")
