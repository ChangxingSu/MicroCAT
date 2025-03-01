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
import taxopy
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

def twist_dict_UMI(nested,taxdb):
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
        for kvalue, count in tdict.items():
            # Check if kvalue has a genus in the same cellbarcode
            taxon = taxopy.Taxon(int(kvalue), taxdb)
            lineage = taxon.taxid_lineage            
            # Add count to parent if it exists
            for parent in lineage:
                if parent in newdict[ckey] and parent != kvalue:
                    newdict[ckey][parent] += count
            
    return(newdict)
def twist_dict(nested,taxdb):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = nested
    # Process genus-species relationships
    for ckey, tdict in newdict.items():
        for kvalue, count in tdict.items():
            # Check if kvalue has a genus in the same cellbarcode
            taxon = taxopy.Taxon(int(kvalue), taxdb)
            lineage = taxon.taxid_lineage            
            # Add count to parent if it exists
            for parent in lineage:
                if parent in newdict[ckey] and parent != kvalue:
                    newdict[ckey][parent] += count
            
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
    parser.add_argument('--nodes', help='Path for nodes.dmp file', required=True)
    parser.add_argument('--names', help='Path for names.dmp file', required=True)
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
    logger.info('Loading taxonomy database', status='run')
    try:
        taxdb = taxopy.TaxDb(nodes_dmp=args.nodes,names_dmp=args.names)
        logger.info('Successfully loading taxonomy database', status='complete')
    except (ValueError, KeyError) as e:
        logger.error(f"An error occurred while loading taxonomy database: {e}")
        sys.exit()

    logger.info('Prasing bam file', status='run')
    total_count = 0
    use_count = 0
    failed_count = 0
    read_taxid_info_dict = dict()
    taxid_info_dict = dict()
    with open(args.align_result, "r") as taxa_file:

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

    logger.info(f'Prasing bam file complete, total reads: {total_count}, use reads: {use_count}, failed reads: {failed_count}, failed rate: {failed_count/total_count}', status='complete')

    logger.info('Checking barcode bam file type', status='run')
    is_cb = False
    is_ub = False
    is_rg = False
    read_count = 0
    with pysam.AlignmentFile(args.cb_bam, "rb") as barcode_file:

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
        with pysam.AlignmentFile(args.cb_bam, "rb") as cb_file:

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
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict_UMI(nested_dict,taxdb))
    if mode == "CB":
        # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
        nested_dict = {}
        with pysam.AlignmentFile(args.cb_bam, "rb") as cb_file:

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
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(nested_dict,taxdb))        
    if mode == "RG":
        # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
        nested_dict = {}
        with pysam.AlignmentFile(args.cb_bam, "rb") as cb_file:

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
        rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(nested_dict,taxdb))
    logger.info(f'Parsing bam file complete, total reads: {total_count}, use reads: {use_count}, skipped reads: {skipped}', status='complete')


    sparsematrix =  csr_matrix((vals, (rows, cols)))
    # Get mpa name for taxonomy ID
    taxname_list = [taxid_info_dict[k] for k in taxid_list]
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
        # Convert dictionary to list for DataFrame creation
        detailed_results = []
        for read_name, info in read_taxid_info_dict.items():
            if "Barcode" in info:  # Only output reads with barcode information
                result = {
                    "Read_Name": read_name,
                    "Taxonomy_ID": info["taxid"],
                    "Species": info["taxname"],
                    "Genus": info["taxgenus"],
                    "Barcode": info["Barcode"]
                }
                if mode == "CB_UMI" and "UMI" in info:
                    result["UMI"] = info["UMI"]
                detailed_results.append(result)
        
        # Create DataFrame and output as TSV
        detailed_df = pd.DataFrame(detailed_results)
        detailed_df.to_csv(args.output_read_tsv, sep="\t", index=False)
        logger.info(f"Detailed TSV with barcode information saved to {args.output_read_tsv}", status="complete")

    logger.info(f'Finish Saving the result', status='Complete')


if __name__ == "__main__":
    main()