import pysam
import re
import collections.abc
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import csv
import logging
import os
import argparse
import multiprocessing as mp
from multiprocessing import freeze_support
import sys



parser = argparse.ArgumentParser(description='This script is used to output kraken2 classified microbial data in cellranger format as feature.tsv,barcodes.tsv,matrix.mtx \n This requires the additional packages pysam(If your python version is up to 3.9)\n')
parser.add_argument('--kraken_output', dest='krak_output',
                    help='Kraken output file.')
parser.add_argument('--bam', dest='bam', 
                    help="The bam file after human host comparison, as input to kraken")
parser.add_argument('--outdir', dest='outdir', default='krak2sc', 
                    help="name of the folder to download the genomes to. If this already exists, the result will be added to it. By default this is krak2sc")
parser.add_argument('--log_file', dest='log_file', default='logfile_krak2sc.log',
                    help="File to write the log to")
parser.add_argument('--processors', dest='proc', default=1,
                    help="Number of processors to use to rename genome files")
parser.add_argument('--nodes_dump', required=True,
    help='Kraken2 database node tree file path')
parser.add_argument('--inspect', required=True,
    dest="inspect_file", help='Kraken2 database inspect file path')

args = parser.parse_args()
bamfile = args.bam
mgfile = args.krak_output
outdir = args.outdir
log_file = args.log_file
n_processors =args.proc


if not os.path.exists(outdir):
    os.system('mkdir -p '+outdir)

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

def mg2sc(bamfile, mgfile, dbfile, outdir,child_parent,taxid_rank):
    """ Main Function. 
    Creates a sparse matrix with transcript count per organism for each cell."""

    # Generate variables based on input
    matrixfile = os.path.join(outdir,'matrix.mtx')
    cellfile = os.path.join(outdir, 'barcodes.tsv')
    taxfile = os.path.join(outdir,'features.tsv')
    # dbfile_out =os.path.join(outdir,'hierarchy.txt')

    # Extract taxonomy IDs for each transcript
    mg_dict = extract_ids(bamfile, mgfile, child_parent, taxid_rank)

    # # Find most frequent taxonomy for each transcript
    # map_nested_dicts(mg_dict, most_frequent)

    # Make sparse matrix
    rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(mg_dict))
    sparsematrix = csr_matrix((vals, (rows, cols)))

    # Get ncbi name for taxonomy ID
    taxdict = krakenID2dict(dbfile, taxid_list)
    # get tax names
    taxname_list = pool.map(taxdict.get, taxid_list)
    # taxname_list = [taxdict[k] for k in taxid_list]

    # store sparse matrix
    mmwrite(matrixfile, sparsematrix)
    
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
    
    # # Store reference database hierarchy
    # with open(dbfile) as f:
    #     with open(dbfile_out, "w") as f1:
    #         for line in f:
    #             f1.write(line) 

## Fix filter kraken
# def extract_ids(bamfile, krakenfile): 
#     """
#     Builds a nested dictionary with KRAKEN2 taxonomy code for each transcript and the cell it belongs to.
#     Input:  Output from KRAKEN2, .bam file with unmapped reads
#     Output: {cellbarcode: {transcriptbarcode: krakentaxonomyID}}
#     """
#     line = 0
#     skipped = 0
#     # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
#     nested_dict = {}
    
#     # Iterate simultanously through bam and kraken file
#     for sread,kread in zip(pysam.AlignmentFile(bamfile, "rb"),open(krakenfile,"r")):
        
#         # count the total number of reads analysed
#         line += 1
        
#         # Check that read names in kraken and bam file match
#         if sread.query_name != kread.split('\t')[1]:
#             skipped += 1
#             logging.warning("sam file read name and metagenomicsfile read name don't match and are therefore excluded: sam: {}, kraken: {}".format(sread.query_name, kread.split('\t')[1]))
#             continue

#         # Get cell barcode and UMI from bam file
#         try:
#             sread_CB = sread.get_tag('CB')
#             sread_UB = sread.get_tag('UB')
#         except:
#             # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
#             skipped += 1
#             continue
            
#         # Get taxonomy ID from kraken file
#         kread_taxid = kread.split('\t')[2]
#         if (type(kread_taxid) != int) and (kread_taxid.isdigit() == False):
#             try:
#                 # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
#                 kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
#             except:
#                 # in this case, something is wrong!
#                 logging.debug("Here is an error. TaxID: {}".format(kread_taxid))
#                 sys.exit()

#         # Make nested dictionary with cells and transcripts
#         if sread_CB in nested_dict:
#             # If cell and transcript exist, add taxonomy ID to list
#             if sread_UB in nested_dict[sread_CB]:
#                 nested_dict[sread_CB][sread_UB].append(kread_taxid)
#             # Otherwise create transcript dictionary for cell
#             else:
#                 nested_dict[sread_CB][sread_UB] = [kread_taxid]
#         else:
#             # if cell doesn't exist, create cell and transcript dictionary with kraken id
#             nested_dict[sread_CB] = {sread_UB: [kread_taxid]}
    
#     # Output control values
#     logging.info("total reads: {}, skipped reads: {}".format(line,skipped))
    
#     return nested_dict


def extract_ids(bamfile, krakenfile,child_parent, taxid_rank): 
    """
    Builds a nested dictionary with KRAKEN2 taxonomy code for each transcript and the cell it belongs to.
    Input:  Output from KRAKEN2, .bam file with unmapped reads
    Output: {cellbarcode: {transcriptbarcode: krakentaxonomyID}}
    """
    line = 0
    skipped = 0
    # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
    nested_dict = {}

    # Load the kraken file into memory
    kraken_data = {}
    with open(krakenfile, "r") as krakenfile:
        for kread in krakenfile:
            parts = kread.split('\t')
            if len(parts) > 1:
                kraken_data[parts[1]] = kread.strip()

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
        kread = kraken_data[sread.query_name]

        # Get cell barcode and UMI from bam file
        try:
            sread_CB = sread.get_tag('RG')
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue
            
        # Get taxonomy ID from kraken file
        kread_taxid = kread.split('\t')[2]
        if (type(kread_taxid) != int) and (kread_taxid.isdigit() == False):
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                kread_taxid = re.search(r'\(taxid (\d+)\)', kread_taxid).group(1)
                # Store as species level id
                kread_taxid = taxid_to_desired_rank(str(kread_taxid), 'species', child_parent, taxid_rank)
            except:
                # in this case, something is wrong!
                logging.debug("Here is an error. Queryname: {}".format(sread.query_name))
                # sys.exit()
                continue

        # Make nested dictionary with RG and taxonomy IDs
        if sread_CB in nested_dict:
            # If RG exists, add taxonomy ID to list 
            nested_dict[sread_CB].append(kread_taxid)
        else:
            # If RG doesn't exist, create list with taxonomy ID
            nested_dict[sread_CB] = [kread_taxid]
        # Output control values
    logging.info("total reads: {}, skipped reads: {}".format(line,skipped))
    
    return nested_dict


def process_line(line, taxids):
    taxid_db, taxname = line[4], line[5].lstrip()
    if taxid_db in taxids:
        return taxid_db, taxname
    else:
        return None

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

def krakenID2dict(dbfile, taxid_list):
    """
    Get name for each taxonomy ID from kraken database
    """
    # create a set of taxid_list
    taxids = set(taxid_list)
    
    # process lines in parallel using multiple processes
    with open(dbfile, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        with mp.Pool() as pool:
            results = pool.starmap(process_line, ((line, taxids) for line in reader))
        
        # filter out None values and create a dictionary
        taxdict = dict(filter(lambda x: x is not None, results))
    
    # add unclassified taxon
    taxdict['0'] = 'unclassified'
    
    return taxdict

# def krakenID2dict(dbfile, taxid_list):
#     """
#     Get name for each taxonomy ID from kraken database
#     """
#     # iterate through inspect file and lookup taxonomy ids
#     k=0
#     taxdict = {'0': 'unclassified'}

#     with open(dbfile) as f:
#         for line in f:
#             if line.startswith("#"):
#                 continue
            
#             # string to list
#             line = line[:-1].split('\t')
#             taxid_db = line[4]
#             taxname = line[5].lstrip()
            
#             if taxid_db in taxid_list:
#                 taxdict[taxid_db] = taxname
    
#     return taxdict

# def extract_taxref(file):
#     """ 
#     Extract taxonomy reference for each read.
#     Input:  viral track output .bam file
#     Output: dictionary with {readname: taxonomy ID}, list of unique taxonomy IDs
#     """
#     # extract taxref for each read
#     tdict = {}
#     line = 0
#     skipped = 0
#     taxref_list = set('0')
    
#     for read in pysam.AlignmentFile(file, "rb"):
#         # count the total number of reads analysed
#         line += 1
#         try:
#             # Extract readname and taxonomy reference
#             taxref = read.to_dict().get('ref_name').split('|')[1]
#             taxref_list.add(taxref)
#             tdict[read.query_name] = taxref
#         except:
#             # in case some reads are unmapped or don't work
#             skipped += 1
#     logging.info("Reads in ViralTrack output: {}, reads without taxonomy reference or that failed: {}".format(line, skipped))
#     return(tdict, taxref_list)


def extract_bc(file):
    """ 
    Extracts cellbarcode and UMI for each readname
    Input:  unmapped .bam file
    Output: dictionary with {readname: [cellbarcode, UMI]}
    """
    # extract UB and CB for each read
    bcdict = {}
    line = 0
    skipped = 0

    for read in pysam.AlignmentFile(file, "rb"):
        # count the total number of reads analysed
        line += 1
        # Get cell barcode and UMI from bam file
        try:
            # Extract readname, cell barcode and UMI
            bcdict[read.query_name] = [read.get_tag('CB'),read.get_tag('UB')]
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

    logging.info("Reads in original bam file: {}, reads without cellbarcode or UMI: {}".format(line, skipped))
    return(bcdict)

if __name__ == "__main__":
    freeze_support()   # required to use multiprocessing
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    # 创建日志处理器
    handler = logging.FileHandler(log_file)

    # 设置日志处理器的级别和格式
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    # 添加日志处理器到logger对象
    logger = logging.getLogger()
    logger.addHandler(handler)
    logger.info('Parsing taxonmy full lineage infomation from NCBI nodes.dump')
    try:
        child_parent, taxid_rank = make_dicts(args.nodes_dump)
        logger.info('Successfully parsing taxonmy full lineage infomation from NCBI nodes.dump')
    except:
        logger.error("Couldn't get the taxonmy full lineage infomation from NCBI nodes.dump")
        sys.exit()
    # multi process
    pool = mp.Pool(processes=int(n_processors))
    mg2sc(bamfile=bamfile, mgfile=mgfile, dbfile=args.inspect_file, outdir=outdir,child_parent=child_parent,taxid_rank=taxid_rank)
    # 关闭进程池
    pool.close()
    pool.join()

    
