#!/usr/bin/env snakemake

import sys
from snakemake.utils import min_version
import os
import numpy as np
import pandas

## beta test
sys.path.append('/data/project/host-microbiome/microcat/microcat/')
import sample

wildcard_constraints:
    Patient = "[a-zA-Z0-9_]+", # Any alphanumeric characters and underscore
    # tissue = "S[0-9]+",  # S followed by any number
    lane = "L[0-9]{3}",  # L followed by exactly 3 numbers
    plate = "P[0-9]{3}",  # L followed by exactly 3 numbers
    library = "[0-9]{3}"  # Exactly 3 numbers


min_version("7.0")
shell.executable("bash")

class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)


PLATFORM = None

if config["params"]["host"]["starsolo"]["do"]:
    if "tenX" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "tenX"
    elif "Smartseq" in config["params"]["host"]["starsolo"]["assay"]:
        PLATFORM = "smartseq"
    else:
        raise ValueError("Platform must be either 'tenX' or 'smartseq'")
elif config["params"]["host"]["cellranger"]["do"]:
    PLATFORM = "tenX"
else:
    raise ValueError("Platform must be either 'tenX' or 'smartseq'")


# def parse_samples(samples_tsv):
#     # Load the sample.tsv file into a pandas DataFrame
#     df = pandas.read_csv(samples_tsv, sep='\t')
    
#     # Split the 'id' column into 'sample' and 'lane' columns at the last underscore '_'
#     # The resulting 'sample' and 'lane' columns are added to the DataFrame
#     df[['sample', 'SX', 'lane']] = df['id'].str.rsplit(pat='_', n=2, expand=True)

#     # Create a 'fastqs_dir' column that contains the directory of the fastq files
#     # This is done by removing the filename from the 'fq1' column
#     df['fastqs_dir'] = df['fq1'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
#     # If the 'lane' column contains any NaN values (due to the absence of a lane in 'id'), replace them with 'no_lane'
#     df['lane'] = df['lane'].fillna('no_lane')  
    
#     # Return a list of dictionaries where each dictionary contains the 'sample', 'lane', and 'fastqs_dir' for a sample
#     return df[['sample', 'lane', 'fastqs_dir']].to_dict(orient='records')

# def get_fastqs_dir(wildcards):

#     # Iterate over each sample in the global '_samples' list
#     for s in _samples:
#         # If the 'sample' field of the current sample dictionary matches the 'sample' wildcard
#         if s['sample'] == sample_name:
#             # Return the 'fastqs_dir' field of the matching sample
#             return s['fastqs_dir']

#     # If no matching sample is found in the '_samples' list, raise a ValueError with an appropriate message
#     raise ValueError(f"No fastqs_dir found for sample: {sample_name}")

    
# def get_starsolo_sample_id(samples, wildcards, read):
#     sample_reads = [s[read] for s in samples if s['sample'] == wildcards.sample]
#     if read == "fq1":
#         return ','.join(sorted(sample_reads))
#     elif read == "fq2":
#         return ' '.join(sorted(sample_reads))


try:
    SAMPLES = sample.parse_samples(config["params"]["samples"],platform = PLATFORM)
    SAMPLES_ID_LIST = SAMPLES.index.get_level_values("sample_id").unique()
except FileNotFoundError:
    warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)


# sample_names = [sample_dict["sample"] for sample_dict in _samples]

# include rules
# include: "../rules/common.smk",
include: "../rules/host.smk"
# include: "../rules/ERCC.smk"
include: "../rules/classfier.smk"



rule all:
    input:
        rules.host_all.input,
        # directory(os.path.join(config["output"]["host"], "unmapped_host/")),
        # os.path.join(
        #             config["output"]["host"],
        #             "starsolo_count/Aligned_out.bam")
        rules.classifier_all.input