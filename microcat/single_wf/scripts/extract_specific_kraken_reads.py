import argparse
import pandas as pd
import re
import os
import sys
import logging
import pysam
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--krak_output_file", action="store", help="path to kraken output file")
    parser.add_argument("--fastq_r1", action="store", help="extract filename path")
    parser.add_argument("--fastq_r2", action="store", help="extract filename path")
    parser.add_argument("--fastq", action="store", help="extract filename path")
    parser.add_argument('--ktaxonomy', required=True,
        help='Kraken2 database ktaxonomy file path')
    parser.add_argument("--keep_original", action="store", default=True, help="delete original bam file? T/F")
    parser.add_argument('--candidate', dest='candidate', 
                        help="candidate species")
    parser.add_argument('--processes', type=int, default=1,help='Number of processes to use')
    parser.add_argument('--input_bam_file', required=True,
        dest='input_bam_file', help='Input origin bam file for denosing')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
    parser.add_argument('--sample_name',
                        type=str,
                        help='One sample name corresponding to the input files')    
    args = parser.parse_args()
    processes = args.processes
    # Set log level based on command line arguments
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # Create a file handler and add the formatter to it
    file_handler = logging.FileHandler(args.log_file)  # Output logs to the specified file
    file_handler.setFormatter(log_format)
    logger.addHandler(file_handler)

    
    logger.info('Reading Kraken report file', status='run')
    candidate_species = pd.read_csv(args.candidate,sep="\t")
    # extract sample name from kraken report file
    candidate_species['sample'] = candidate_species['sample'].str.replace('_krak_sample_denosing', '')
    # select the sample
    candidate_species = candidate_species.loc[candidate_species['sample'] == args.sample_name]
    # print(candidate_species)
    # sys.exit()
    desired_taxid_list = set(candidate_species["main_level_taxid"].unique())

    logger.info('Finishing reading Kraken report file', status='complete')
    logger.info('Parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='run')
    try:
        taxid2node = make_dicts(args.ktaxonomy)
        logger.info('Successfully parsing taxonmy full lineage infomation from Kraken ktaxonomy', status='complete')
    except:
        logger.error("Couldn't get the taxonmy full lineage infomation from NCBI nodes.dump")
        sys.exit()

    # lineage_dict = {}
    unique_taxid_set = set()
    for main_tax_id in desired_taxid_list:
        try:
            lineage_taxid_list = taxid2node[str(main_tax_id)].lineage_to_desired_rank("G")
            unique_taxid_set.update(lineage_taxid_list)
        except (ValueError, KeyError) as e:
            print("Error occur:", e)
    # unique_taxid_set.update(["78331","3073562","120793","1768","1866885","2642494","1936029"])
    # rtl_dict = {}
    # for species_tax_id in desired_taxid_list:
    #     descendants_ascendants_taxid_list = []
    #     descendants_ascendants_taxid_list.append(species_tax_id)
    #     descendants_ascendants_taxid_list.append(taxid2node[species_tax_id].parent.taxid)
    #     descendants_nodes_list = taxid2node[species_tax_id].children
    #     while len(descendants_nodes_list) > 0:
    #         #For this node
    #         curr_n = descendants_nodes_list.pop()
    #         descendants_ascendants_taxid_list.append(curr_n.taxid)
    #     rtl_dict[species_tax_id] = descendants_ascendants_taxid_list
    kraken_output_query_names = set()
    logger.info('Extracting Bacteria, Fungi, Viruses ncbi taxID', status='run')
    logger.info('Reading kraken output file', status='run')
    # Read krak2 output file and create a copy
    # krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type', 'query_name', 'taxid_info', 'len', 'kmer_position'])

    with open(args.krak_output_file, 'r') as kfile:
        for kraken_line in kfile:
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
            except (ValueError, KeyError) as e:
                # in this case, something is wrong!
                logger.error(f"An error occurred while processing the Kraken output file: {e}")
                logger.error(f"Here is an error. Queryname: {query_name}")
                continue
            if read_type == "C":
                try:
                    kread_taxid = re.search(r'\(taxid (\d+)\)', taxid_info).group(1)
                except (ValueError, KeyError) as e:
                    # in this case, something is wrong!
                    logger.error(f"An error occurred while processing the Kraken output file: {e}")
                    logger.error(f"Here is an error. Queryname: {query_name}")
                    continue
                if kread_taxid in unique_taxid_set:
                    kraken_output_query_names.add(query_name)
            # else:                        
            #     kmer_taxid = set(info.split(':')[0] for info in kmer_position.strip().split())
            #     #Skip if reads are human/artificial/synthetic
            #     if len(kmer_taxid) == 1 and "0" in kmer_taxid:
            #         continue
            #     if "9606" in kmer_taxid or "9605" in kmer_taxid:
            #         continue
            #     if kmer_taxid & unique_taxid_set:
            #         kraken_output_query_names.add(query_name)

    # # 将Kraken的DataFrame的query_name列转换为一个集合
    # kraken_output_query_names = set(krak2_output_filtered["query_name"])
    
    logger.info(f'Extract classified reads from bam file', status='run')
    read_count = 0
    krak_count = 0
    # 打开源BAM文件和目标BAM文件
    with pysam.AlignmentFile(args.input_bam_file, "rb") as source_bam:
        # 初始化计数器
        read_count = 0

        # 检查文件是单端还是双端
        is_paired_end = False
        for sread in source_bam:
            if sread.is_paired:
                is_paired_end = True
                break

    if is_paired_end:
        with pysam.AlignmentFile(args.input_bam_file, "rb",threads=processes) as source_bam:
            # 初始化计数器
            read_count = 0
            krak_count = 0

            # 创建输出的 Fastq 文件（read1 和 read2）
            with open(args.fastq_r1, "a") as output_fastq_r1, open(args.fastq_r2, "a") as output_fastq_r2:
                # 遍历源BAM文件的每一个read
                for sread in source_bam:
                    read_count += 1
                    if sread.query_name in kraken_output_query_names:
                        if sread.is_read1:
                            # 将 read1 写入到第一个 Fastq 文件
                            output_fastq_r1.write(f"@{sread.query_name}\n{sread.query_sequence}\n+\n{sread.qual}\n")
                        elif sread.is_read2:
                            # 将 read2 写入到第二个 Fastq 文件
                            output_fastq_r2.write(f"@{sread.query_name}\n{sread.query_sequence}\n+\n{sread.qual}\n")
                        krak_count += 1
        with open(args.fastq, "w") as output_fastq:
            pass
    else:
        # 打开源BAM文件和目标BAM文件
        with pysam.AlignmentFile(args.input_bam_file, "rb",threads=processes) as source_bam:
            # 初始化计数器
            read_count = 0
            krak_count = 0
            with open(args.fastq, "a") as output_fastq:
                for sread in source_bam:
                    read_count += 1
                    if sread.query_name in kraken_output_query_names:
                        output_fastq.write(f"@{sread.query_name}\n{sread.query_sequence}\n+\n{sread.qual}\n")    
                        krak_count += 1
        with open(args.fastq_r1, "w") as output_fastq_r1, open(args.fastq_r2, "w") as output_fastq_r2:
            pass

    # 日志记录分类的reads提取完成
    logger.info(f'Extract classified reads from bam file', status='complete')
    logger.info(f'Total unmapped reads: {read_count}', status='summary')
    logger.info(f'Total unmapped reads classified as bactreia, virus ,archaea and fungi by Kraken: {krak_count}', status='summary')
    print('Done')

if __name__ == "__main__":
    main()
