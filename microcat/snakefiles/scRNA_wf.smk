#!/usr/bin/env snakemake

import sys
from snakemake.utils import min_version
import os
import numpy as np
import pandas

min_version("7.0")
shell.executable("bash")

class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)


def parse_samples(samples_tsv):
    # Load the sample.tsv file into a pandas DataFrame
    df = pandas.read_csv(samples_tsv, sep='\t')
    
    # Split the 'id' column into 'sample' and 'lane' columns at the last underscore '_'
    # The resulting 'sample' and 'lane' columns are added to the DataFrame
    df[['sample', 'SX', 'lane']] = df['id'].str.rsplit(pat='_', n=2, expand=True)

    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    # This is done by removing the filename from the 'fq1' column
    df['fastqs_dir'] = df['fq1'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # If the 'lane' column contains any NaN values (due to the absence of a lane in 'id'), replace them with 'no_lane'
    df['lane'] = df['lane'].fillna('no_lane')  
    
    # Return a list of dictionaries where each dictionary contains the 'sample', 'lane', and 'fastqs_dir' for a sample
    return df[['sample', 'lane', 'fastqs_dir']].to_dict(orient='records')

def get_fastqs_dir(wildcards):

    # Iterate over each sample in the global '_samples' list
    for s in _samples:
        # If the 'sample' field of the current sample dictionary matches the 'sample' wildcard
        if s['sample'] == sample_name:
            # Return the 'fastqs_dir' field of the matching sample
            return s['fastqs_dir']

    # If no matching sample is found in the '_samples' list, raise a ValueError with an appropriate message
    raise ValueError(f"No fastqs_dir found for sample: {sample_name}")

    
def get_starsolo_sample_id(samples, wildcards, read):
    sample_reads = [s[read] for s in samples if s['sample'] == wildcards.sample]
    if read == "fq1":
        return ','.join(sorted(sample_reads))
    elif read == "fq2":
        return ' '.join(sorted(sample_reads))


try:
    _samples = parse_samples(config["params"]["samples"])
except FileNotFoundError:
    warning(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)


sample_names = [sample_dict["sample"] for sample_dict in _samples]

# def get_all_inputs(wildcards):
#     return expand(
#         [
#             os.path.join(
#                 config["output"]["host"],
#                 "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"
#             ),
#             os.path.join(
#                 config["output"]["host"],
#                 "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq"
#             ),
#         ],
#         sample=sample_names,
#         barcode=[barcode for sample_name in sample_names for barcode in aggregate_CB_bam_output(sample_name)]
#     )
# include rules
# include: "../rules/common.smk",
include: "../rules/host.smk",
# include: "../rules/ERCC.smk"
# include: "../rules/classfier.smk"





# classifier_DO = any([
#     config["params"]["classifier"]["kraken2uniq"]["do"],
#     config["params"]["classifier"]["krakenuniq"]["do"]
# ])


# rule all:
#     input:
#         # expand(os.path.join(
#         #     config["output"]["host"],
#         #     "cellranger_count/{sample}/{sample}_unmappped2human_CB_sorted_bam.bam"),sample=[s['sample'] for s in _samples])
#         directory(expand(os.path.join(
#                 config["output"]["host"],
#                 "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/")),sample=[s['sample'] for s in _samples])
        # rules.classifier_all.input
# rule all:
#     input:
#         expand([
#             "{outdir}/krakenuniq_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_output.txt",
#             "{outdir}/krakenuniq_report/custom/{sample}/cell_level/{sample}_{barcode}_krakenuniq_report.txt",
#             "{outdir}/krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_classified.fq",
#             "{outdir}/krakenuniq_classified_output/{sample}/cell_level/{sample}_{barcode}_krakenuniq_unclassified.fq"
#         ], outdir=config["output"]["classifier"], sample=sample_names, barcode=get_barcodes(wildcards.sample))

rule all:
    input:
        rules.host_all.input

        # expand([
        #     os.path.join(
        #     config["output"]["host"],
        #     "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R1.fastq"),
        #     os.path.join(
        #     config["output"]["host"],
        #     "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}_R2.fastq")],
        #         sample=glob_wildcards(os.path.join(config["output"]["host"], "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}.bam")).sample,
        #         barcode=glob_wildcards(os.path.join(config["output"]["host"], "cellranger_count/{sample}/unmapped_bam_CB_demultiplex/CB_{barcode}.bam")).barcode
        #         )
