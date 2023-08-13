import argparse
import pandas as pd
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os
import sys
# def grep_batch_taxids(batch_taxids, krak_output_file, extract_file):
#     with open(extract_file_path, 'w') as extract_file:
#         try:
#             taxid_str = "|".join([f"(taxid {t})" for t in batch_taxids])
#             subprocess.run(["grep", "-E", "-w", taxid_str, krak_output_file], stdout=extract_file, check=True)
#         except subprocess.CalledProcessError as e:
#             print(f"Error occurred while searching for batch of taxids: {e}")


# def batch_search_taxids(taxids, krak_output_file, extract_krak_file, ntaxid, num_cores=None):
#     num_cores = num_cores or os.cpu_count()  # Default to the number of available CPU cores
#     with open(extract_krak_file, "w") as extract_file:
#         with ThreadPoolExecutor(max_workers=num_cores) as executor:
#             batch_starts = range(0, len(taxids), ntaxid)
#             batch_ends = [min(start + ntaxid, len(taxids)) for start in batch_starts]
#             batch_taxids = [taxids[start:end] for start, end in zip(batch_starts, batch_ends)]

#             executor.map(grep_batch_taxids, batch_taxids, [krak_output_file]*len(batch_taxids))
#             # wait
#             executor.shutdown()
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--krak_output_file", action="store", help="path to kraken output file")
    parser.add_argument("--kraken_report",dest="krak_report_file",action="store", help="path to kraken report")
    parser.add_argument("--mpa_report",dest="mpa_report_file", action="store", help="path to standard kraken mpa report")
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

    # 用下划线替换中间空格
    kr['scientific_name'] = kr['scientific name'].str.replace(r' ', '_')
    mpa = pd.read_csv(args.mpa_report_file, sep='\t', names=['mpa_taxa','reads'])
    mpa['taxa'] = mpa['mpa_taxa'].apply(lambda x: re.sub(r'[a-z]__', '', x.split('|')[-1]))

    # we only focus on  k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea
    keywords = ["k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea"]

    mpa_subset = mpa[mpa['mpa_taxa'].str.contains('|'.join(keywords))]
    df = pd.merge(mpa_subset, kr, left_on='taxa', right_on='scientific_name')
    krak_filtered = kr[kr['scientific_name'].isin(df['scientific_name'])] 
    taxid = krak_filtered['ncbi_taxa'].tolist()

    taxid_list = [taxid[i:i+args.ntaxid] for i in range(0, len(taxid), args.ntaxid)]

    if os.path.exists(args.extract_krak_file):
        os.remove(args.extract_krak_file)

    if not os.path.exists(args.extract_krak_file):
        open(args.extract_krak_file, 'w').close()
    batch_search_taxids(taxid_list,krak_output_file = args.krak_output_file,extract_krak_file = args.extract_krak_file,ntaxid = args.ntaxid, num_cores=args.cores)  # Specify the number of cores here

    if not args.keep_original:
        os.remove(args.krak_output_file)

    print('Done')

if __name__ == "__main__":
    main()
