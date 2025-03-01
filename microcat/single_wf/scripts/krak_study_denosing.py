import argparse
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
import numpy as np
from statsmodels.stats.multitest import multipletests
import logging
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


def read_kraken_reports(files, sample_names=None, study_name=None, min_reads=2, min_uniq=2):
    """
    Read Kraken reports from files and return a DataFrame with the data.

    Parameters:
        files (list): List of file paths containing Kraken reports.
        sample_names (list, optional): List of sample names corresponding to the input files. Default is None.
        study_name (str, optional): Name of the study. Default is None.
        min_reads (int, optional): Minimum number of reads per taxon. Default is 2.
        min_uniq (int, optional): Minimum number of unique sequences per taxon. Default is 2.
        path (str, optional): Path to the files. Default is '.'.

    Returns:
        pd.DataFrame: DataFrame containing the combined data from Kraken reports.
    """
    if sample_names is None:
        sample_names = [f.stem for f in files]  # Use file names without extension as sample names
    if study_name is None:
        study_name = [None] * len(files)

    df = []
    for i, f in enumerate(files):
        try:
            tmp = pd.read_csv(f, sep="\t")
        except pd.errors.EmptyDataError:
            logger.warning(f"Empty file: {f}. Skipping this file.")
            continue

        tmp['scientific name'] = tmp['scientific name'].str.strip()
        tmp_df = pd.DataFrame({
            'study': study_name[i],
            'sample': sample_names[i],
            'rank': tmp['classification_rank'],
            'ncbi_taxa': tmp['ncbi_taxa'],
            'main_level_taxid': tmp['main_level_taxid'],
            'sci_name': tmp['scientific name'],
            'reads': tmp['fragments'],
            'minimizers': tmp['max_minimizers'],
            'uniqminimizers': tmp['max_uniqminimizers'],
            'classification_rank': tmp['classification_rank'],
            'genus_level_taxid': tmp['genus_level_taxid'],
            'superkingdom': tmp['superkingdom']
        })

        df.append(tmp_df)

    df = pd.concat(df, ignore_index=True)
    logger.info(f"Successfully read {len(df)} records from {len(files)} files.",status="summary")
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path',
                        type=str,
                        help='One or more file path containing custom style Kraken reports')
    parser.add_argument('--out_path',
                        type=str,
                        help='Result output path')
    # parser.add_argument('--sample_name',
    #                     type=str,
    #                     help='One sample name corresponding to the input files')
    parser.add_argument('--study_name',
                        type=str,
                        help='Name of the study')
    parser.add_argument('--min_reads',
                        type=int,
                        default=2,
                        help='Minimum number of reads per taxon')
    parser.add_argument('--min_uniq',
                        type=int,
                        default=2,
                        help='Minimum number of unique sequences per taxon')
    # parser.add_argument('--cell_line',
    #                     type=str,
    #                     help='Cell line path')
    parser.add_argument('--raw_file_list', nargs='+',help='sample raw file list path')
    parser.add_argument('--file_list', nargs='+',help='sample denosing file list path')
    parser.add_argument('--log_file', dest='log_file', 
        required=True, default='logfile_download_genomes.txt',
        help="File to write the log to")
    parser.add_argument('--verbose', action='store_true', help='Detailed print')
    
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


    if args.path:
        path = args.path
        path = Path(path)
        files = list(path.glob('**/*_krak_sample_denosing.txt'))  
    elif args.file_list:
        # Convert the file paths to Path objects
        files = [Path(file_path) for file_path in args.file_list]
    else:
        raise ValueError("Either --path or --file_list must be provided")

    out_path = args.out_path
    # sample_name = args.sample_name
    study_name = args.study_name
    min_reads = args.min_reads
    min_uniq = args.min_uniq
    # celline_file = args.cell_line


    logger.info('Reading kraken sample denosing results', status='run')
    # Read the all krak report
    kraken_reports_all = read_kraken_reports(files, sample_names=None,study_name=study_name, min_reads=min_reads, min_uniq=min_uniq)


    logger.info('Finishing reading kraken sample denosing results', status='complete')
    logger.info('Checking sample number', status='run')


    kraken_reports_all_species = kraken_reports_all.copy()

    # candidate_species_all = kraken_reports_all.loc[kraken_reports_all["classification_rank"] == "S"]
    candidate_species_all = kraken_reports_all[kraken_reports_all['classification_rank'].str.startswith('S')]
    logger.info(f'Saving the result', status='run')
    # Save the filtered data to CSV
    # filter_kraken_reports_specific.to_csv(out_path, sep='\t', index=False)
    candidate_species_all.to_csv(args.out_path, sep='\t', index=False)

if __name__ == "__main__":
    main()