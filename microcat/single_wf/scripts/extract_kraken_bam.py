import argparse
import pandas as pd
import re
import os
import sys
import pysam 
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



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--krak_output_file", action="store", help="path to kraken output file")
    parser.add_argument("--kraken_report",dest="krak_report_file",action="store", help="path to kraken report")
    parser.add_argument("--mpa_report",dest="mpa_report_file", action="store", help="path to standard kraken mpa report")
    parser.add_argument("--extract_krak_file", action="store", help="extract filename path")
    parser.add_argument("--keep_original", action="store", default=True, help="delete original bam file? T/F")
    parser.add_argument('--input_bam_file', required=True,
        dest='input_bam_file', help='Input origin bam file for denosing')
    parser.add_argument('--extracted_bam_file', required=True,
        dest='extracted_bam_file', help='Input origin bam file for denosing')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
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

    
    logger.info('Reading Kraken report file', status='run')
    kr = pd.read_csv(args.krak_report_file, sep='\t',names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
    # removing root and unclassified taxa
    kr = kr.iloc[2:]
    # 去除两端空格
    kr['scientific name'] = kr['scientific name'].str.strip() 
    logger.info('Finishing reading Kraken report file', status='complete')
    # 用下划线替换中间空格
    kr['scientific_name'] = kr['scientific name'].str.replace(r' ', '_')
    logger.info('Reading Kraken mpa style report file', status='run')
    mpa = pd.read_csv(args.mpa_report_file, sep='\t', names=['mpa_taxa','reads'])
    mpa['taxa'] = mpa['mpa_taxa'].apply(lambda x: re.sub(r'[a-z]__', '', x.split('|')[-1]))
    logger.info('Extracting Bacteria, Fungi, Viruses and Archaea ncbi taxID', status='run')
    # we only focus on  k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea
    keywords = ["k__Bacteria", "k__Fungi", "k__Viruses","k__Archaea"]

    mpa_subset = mpa[mpa['mpa_taxa'].str.contains('|'.join(keywords))]
    df = pd.merge(mpa_subset, kr, left_on='taxa', right_on='scientific_name')
    krak_filtered = kr[kr['scientific_name'].isin(df['scientific_name'])] 
    logger.info('Reading kraken output file', status='run')
    # Read krak2 output file and create a copy
    krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type', 'query_name', 'taxid_info', 'len', 'kmer_position'])

    # Extract 'taxa' and 'taxid' from 'taxid_info' column
    krak2_output[['taxa', 'taxid']] = krak2_output['taxid_info'].str.extract(r'(.*) \(taxid (\d+)\)')
    krak2_output['taxid'] = krak2_output['taxid'].str.replace(r'\)', '').str.strip()
    krak2_output['taxid'] =krak2_output['taxid'].astype(str)

    krak_filtered['ncbi_taxa'] = krak_filtered['ncbi_taxa'].astype(str)
    # Filter krak2_output to keep only rows with taxid appearing in krak_study_denosing ncbi_taxa
    krak2_output_filtered = krak2_output[krak2_output['taxid'].isin(krak_filtered['ncbi_taxa'])]

    # 将Kraken的DataFrame的query_name列转换为一个集合
    kraken_output_query_names = set(krak2_output_filtered["query_name"])

    logger.info(f'Extract classified reads from bam file', status='run')
    read_count = 0
    krak_count = 0
    # 打开源BAM文件和目标BAM文件
    with pysam.AlignmentFile(args.input_bam_file, "rb") as source_bam, \
        pysam.AlignmentFile(args.extracted_bam_file, "wb", header=source_bam.header) as output_bam:
        
        # 遍历源BAM文件的每一个read
        for sread in source_bam:
            read_count += 1
            if sread.query_name in kraken_output_query_names:
                # 如果query name匹配Kraken的输出，将read写入目标BAM文件
                output_bam.write(sread)
                krak_count += 1
    # 日志记录分类的reads提取完成
    logger.info(f'Extract classified reads from bam file', status='complete')
    logger.info(f'Total unmapped reads: {read_count}', status='summary')
    logger.info(f'Total unmapped reads classified as bactreia, virus ,archaea and fungi by Kraken: {krak_count}', status='summary')
    print('Done')

if __name__ == "__main__":
    main()
