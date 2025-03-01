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
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import logging

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
    # species_tid = int(ref_name.split(":")[0])
    species_tid = int(acc2tax[ref_name])
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
            identity = (alignment.query_alignment_length - alignment.get_tag("NM")) / alignment.infer_query_length()
            if alignment.reference_name and align_len and identity > 0:
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
            if alignment.reference_name and align_len:
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


def expectation_maximization_iterations(log_p_rgs, freq, lli_thresh, input_threshold):
    """
    执行EM算法迭代，包含迭代质量保持阈值步骤以解决长尾效应
    
    Args:
        log_p_rgs: 读段对应的物种对数似然字典 {read_id: ([taxids], [scores])}
        freq: 初始物种丰度估计
        lli_thresh: 对数似然增量阈值
        input_threshold: 物种丰度阈值
    
    Returns:
        freq_full: 所有物种的丰度
        freq_set_thresh: 经过阈值处理后的物种丰度
        p_sgr: 读段-物种分配概率字典
    """
    n_reads = len(log_p_rgs)
    if n_reads == 0:
        return freq, freq, {}
    
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
    
    while True:
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
            # 计算读段r归属于各个物种的概率
            p_r = {}  # 归一化概率p(s|r)
            max_log_p = -float('inf')
            
            # 找出最大对数概率以防数值溢出
            for i, g in enumerate(taxids):
                if strain_valid[g]:
                    log_p = log_scores[i]
                    if log_p + math.log(freq[g]) > max_log_p:
                        max_log_p = log_p + math.log(freq[g])
            
            # 计算分母(归一化因子)
            log_denom = 0
            first = True
            for i, g in enumerate(taxids):
                if strain_valid[g]:
                    log_p = log_scores[i]
                    log_num = log_p + math.log(freq[g]) - max_log_p
                    if first:
                        log_denom = log_num
                        first = False
                    else:
                        log_denom = np.logaddexp(log_denom, log_num)
            
            # 计算每个有效物种的后验概率p(s|r)
            p_r = {}
            for i, g in enumerate(taxids):
                if strain_valid[g]:
                    log_p = log_scores[i]
                    log_num = log_p + math.log(freq[g]) - max_log_p
                    p_r[g] = math.exp(log_num - log_denom)
                    strain_read_count[g] += p_r[g]
            
            p_sgr[r] = p_r
            
            # 更新对数似然
            log_likelihood += max_log_p + log_denom
        
        # Maximization步骤 - 更新丰度
        total_reads_assigned = sum(strain_read_count.values())
        for g in freq:
            if strain_valid[g]:
                freq[g] = strain_read_count[g] / total_reads_assigned
            else:
                freq[g] = 0
        
        # 检查收敛性
        if abs(log_likelihood - prev_log_likelihood) < lli_thresh:
            break
        
        prev_log_likelihood = log_likelihood
        counter += 1
    
    # 最终过滤结果
    freq_full = freq.copy()
    freq_set_thresh = {k: v for k, v in freq.items() if strain_valid[k]}
    
    # 归一化最终结果
    total = sum(freq_set_thresh.values())
    if total > 0:
        freq_set_thresh = {k: v/total for k, v in freq_set_thresh.items()}
    
    return freq_full, freq_set_thresh, p_sgr


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


# 添加参数解析
parser = argparse.ArgumentParser(description='执行EM算法进行物种丰度估计')
parser.add_argument('--bam_file', required=True, help='BAM文件路径')
parser.add_argument('--taxonomy_file', required=True, help='分类文件路径')
parser.add_argument('--output', required=True, help='输出文件路径')
# parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='日志级别')
parser.add_argument('--log_file', default='emu_em.log', help='日志文件路径')
parser.add_argument('--verbose', action='store_true', help='是否输出详细日志信息')
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
acc2tax = {}
with open(taxonomy_file, "r") as f:
    for line in f:
        # skip header
        if line.startswith("accession"):
            continue
            
        value_list = line.strip().split("\t")

        if len(value_list) == 3:
            # accession, taxid, _ = value_list[0], value_list[1], value_list[2]
            acc, accession, taxid = value_list[0], value_list[1], value_list[2]
        elif len(value_list) == 2:
            accession, taxid = value_list[0], value_list[1]
        else:
            print(f"Invalid taxonomy file format at line: {line}")
            continue        
        
        acc2tax[accession] = int(taxid)


logger.info('Loading taxonomy database', status='done')

# Calculate log likelihood of reads given sequences
logger.info('Calculating log likelihood of reads given sequences', status='run')
log_prob_rgs, counts_unassigned, counts_assigned = log_prob_rgs_dict(
    bam_file, log_prob_cigar_op, longest_align_dict, locs_p_cigar_zero)


logger.info('Calculating log likelihood of reads given sequences', status='done')

# Run the EM algorithm
db_ids = list(acc2tax.values())
n_db = len(db_ids)
n_reads = len(log_prob_rgs)
print("Assigned read count: {}\n".format(n_reads))
# check if there are enough reads


freq, counter = dict.fromkeys(db_ids, 1 / n_db), 1

logger.info('Running EM algorithm', status='run')
f_full, f_set_thresh, read_dist = expectation_maximization_iterations(log_prob_rgs,
                                                                      freq, 0.01, 0.00001)

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

# get the taxonomy name from the taxonomy id
def parse_names_dmp(file_path):
    tax_id_to_name = {}
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('|')
            tax_id = int(fields[0].strip())
            name = fields[1].strip()
            name_class = fields[3].strip()
            if name_class == 'scientific name':
                tax_id_to_name[tax_id] = name
    return tax_id_to_name


names_dmp_file = "/data/comics-sucx/database/microcat/taxonomy/names.dmp"

# 解析names.dmp文件
tax_id_to_name = parse_names_dmp(names_dmp_file)


# 将结果写入文件
logger.info(f"开始将read分类结果写入文件: {read_classifications_file}")

with open(output, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["Read_Name", "Taxonomy_ID", "Probability", "Species", "Genus"])
    
    for read_name, (tax_id, prob) in best_classifications.items():
        tax_name = tax_id_to_name.get(tax_id, "Unknown")
        tax_genus = tax_name.split(" ")[0] if tax_name != "Unknown" else "Unknown"
        writer.writerow([read_name, tax_id, prob, tax_name, tax_genus])

logger.info(f"Read分类结果已保存到文件: {read_classifications_file}")