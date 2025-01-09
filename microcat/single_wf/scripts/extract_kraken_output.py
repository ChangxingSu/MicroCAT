import argparse
import pandas as pd
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os
import sys
def grep_batch_taxids(batch_taxids, krak_output_file, extract_file_path):
    flattened_taxids = [taxid for sublist in batch_taxids for taxid in sublist]
    taxid_commands = "\|".join([f"(taxid {t})" for t in flattened_taxids])
    command = f"grep -w \"{taxid_commands}\" {krak_output_file}"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, check=True)

    with open(extract_file_path, 'w') as extract_file:
        extract_file.write(result.stdout.decode())

def batch_search_taxids(taxids, krak_output_file, extract_krak_file, ntaxid, num_cores=None):
    num_cores = num_cores or os.cpu_count()
    
    batch_starts = range(0, len(taxids), ntaxid) 
    batch_ends = [min(start + ntaxid, len(taxids)) for start in batch_starts]
    batch_taxids = [taxids[start:end] for start, end in zip(batch_starts, batch_ends)]

    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        executor.map(grep_batch_taxids, 
                     batch_taxids, 
                     [krak_output_file]*len(batch_taxids),
                     [extract_krak_file]*len(batch_taxids)
                     )
                     
        executor.shutdown()
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
    parser.add_argument("--kraken_report",dest="krak_report_file",action="store", help="path to kraken report")
    parser.add_argument('--ktaxonomy', required=True,
        help='Kraken2 database ktaxonomy file path')
    parser.add_argument("--extract_krak_file", action="store", help="extract filename path")
    parser.add_argument("--keep_original", action="store", default=True, help="delete original fastq file? T/F")
    parser.add_argument("--ntaxid", action="store", type=int, default=8000, help="number of taxids to extract at a time")
    parser.add_argument("--cores", action="store", type=int, default=8, help="number of cores at a time")

    args = parser.parse_args()
    args.krak_output_file = os.path.abspath(args.krak_output_file)
    args.extract_krak_file = os.path.abspath(args.extract_krak_file)

    kr = pd.read_csv(args.krak_report_file, sep='\t',names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
    # removing root and unclassified taxa
    kr = kr.iloc[2:]
    # 去除两端空格
    kr['scientific name'] = kr['scientific name'].str.strip() 
    try:
        taxid2node = make_dicts(args.ktaxonomy)
    except:
        print("Couldn't get the taxonmy full lineage infomation from NCBI nodes.dump")
        sys.exit()

    kr['is_microbiome'] = kr.apply(lambda x: taxid2node[str(x['ncbi_taxa'])].is_microbiome(), axis=1)
    krak_filtered = kr[kr["is_microbiome"]==True]
    taxid = krak_filtered['ncbi_taxa'].tolist()


    if os.path.exists(args.extract_krak_file):
        os.remove(args.extract_krak_file)

    if not os.path.exists(args.extract_krak_file):
        open(args.extract_krak_file, 'w').close()
    batch_search_taxids(taxid_list,krak_output_file = args.krak_output_file,extract_krak_file = args.extract_krak_file,ntaxid = args.ntaxid, num_cores=args.cores)  # Specify the number of cores here

    with open(args.krak_output_file, 'r') as kfile_in:
        with open(args.extract_krak_file, 'r') as kfile_out:

            for kraken_line in kfile:
                read_count += 1
                try:
                    # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                    # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                    # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                    read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                    tax_id = str(re.search(r'\(taxid (\d+)\)', taxid_info).group(1))
                except:
                    # in this case, something is wrong!
                    print("Here is an error. Queryname: {}".format(query_name))
                    # sys.exit()
                    continue
                if tax_id in :
                    krak_count += 1
                    kfile_out.write(kraken_line)

    if not args.keep_original:
        os.remove(args.krak_output_file)

    print('Done')

if __name__ == "__main__":
    main()
