#!/usr/bin/env python

import pandas as pd
import re
import glob
import os
import sys
import json
from ruamel.yaml import YAML



SRA_HEADERS = {
    "PE": "sra_pe",
    "SE": "sra_se",
}


BAM_HEADERS = {
    "PE": "bam_pe",
    "SE": "bam_se",
}

FQ_HEADERS = {
    "PE_CB": "CB",
    "PE_cDNA": "cDNA",
    "SE": "single_reads",
}


HEADERS = {
    "SRA": SRA_HEADERS,
    "BAM": BAM_HEADERS,
    "FQ": FQ_HEADERS
}

def parse_samples(sample_tsv, platform):
    samples_df = pd.read_csv(sample_tsv, sep="\t")
    
    #Check if id, fq1, fq2 columns exist
    if not set(['id', 'fq1', 'fq2']).issubset(samples_df.columns):
        raise ValueError("Columns 'id', 'fq1', 'fq2' must exist in the sample.tsv")


    # Extract library, lane, and plate from id
    samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
    
    # Determine the platform and parse accordingly
    if platform == 'lane':
        samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
        samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
        # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
        # and is always the last part in patient_tissue
        samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        samples_df['tissue'] = samples_df['patient_tissue'].apply(lambda x: x.split('_')[-1])
        samples_df['patient'] = samples_df['patient_tissue'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
    elif platform == 'plate':
        samples_df['is_plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("P"))
        samples_df.loc[samples_df['is_plate'], 'plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
        samples_df['patient_tissue_cell'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
        # and is always the last part in patient_tissue
        samples_df['tissue'] = samples_df['patient_tissue_cell'].str.extract(r'(S\d+)_')
        # 提取patient和cell
        samples_df[['patient', 'cell']] = samples_df['patient_tissue_cell'].str.extract(r'(.+)_S\d+_(.+)')
        samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
        samples_df = samples_df.drop(columns=['patient_tissue_cell'])
        samples_df['patient_tissue'] = samples_df['patient'] + '_' + samples_df['tissue']
    else:
        raise ValueError("Platform must be either 'lane' or 'plate'")


    if samples_df[['patient_tissue', 'library']].isnull().values.any():
        raise ValueError(f"id column must follow the format '{{Patient}}_{{tissue}}_{{lane or plate}}_{{library}}' for platform {platform}")
    
    # Create sample identifier
    samples_df['sample_id'] = samples_df['patient_tissue']
    
    # Check if sample names contain "."
    if samples_df['sample_id'].str.contains("\\.").any():
        raise ValueError("Sample names must not contain '.', please remove '.'")
    
    # Determine if the sequencing is paired-end or single-end
    samples_df['seq_type'] = 'single-end'
    samples_df.loc[samples_df['fq2'].notnull(), 'seq_type'] = 'paired-end'
    
    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    samples_df['fastqs_dir'] = samples_df['fq1'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # Set index
    if platform == 'lane':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "lane", "library"])
    elif platform == 'plate':
        samples_df = samples_df.set_index(["sample_id","patient", "cell","tissue", "plate", "library"])

    # Check if fq1 and fq2 files exist
    for _, row in samples_df.iterrows():
        fq1_exists = os.path.isfile(row['fq1'])
        if not fq1_exists:
            raise FileNotFoundError(f"File not found: {row['fq1']}")
        
        if row['seq_type'] == 'paired-end':
            fq2_exists = os.path.isfile(row['fq2'])
            if not fq2_exists:
                raise FileNotFoundError(f"File not found: {row['fq2']}")

    return samples_df


def parse_bam_samples(sample_tsv, platform):
    samples_df = pd.read_csv(sample_tsv, sep="\t")
    
    is_bam = False

    # Check if id, fq1, fq2 columns exist
    if not set(['id', 'bam']).issubset(samples_df.columns):
        raise ValueError("Columns 'id', 'bam', must exist in the sample.tsv")

    
    # Extract library, lane, and plate from id
    samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
    
    # Determine the platform and parse accordingly
    if platform == 'lane':
        samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
        samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
        # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
        # and is always the last part in patient_tissue
        samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        samples_df['tissue'] = samples_df['patient_tissue'].apply(lambda x: x.split('_')[-1])
        samples_df['patient'] = samples_df['patient_tissue'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
    elif platform == 'plate':
        samples_df['is_plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("P"))
        samples_df.loc[samples_df['is_plate'], 'plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
        samples_df['patient_tissue_cell'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
        # and is always the last part in patient_tissue
        samples_df['tissue'] = samples_df['patient_tissue_cell'].str.extract(r'(S\d+)_')
        # 提取patient和cell
        samples_df[['patient', 'cell']] = samples_df['patient_tissue_cell'].str.extract(r'(.+)_S\d+_(.+)')
        samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
        samples_df = samples_df.drop(columns=['patient_tissue_cell'])
        samples_df['patient_tissue'] = samples_df['patient'] + '_' + samples_df['tissue']
    else:
        raise ValueError("Platform must be either 'lane' or 'plate'")

    if samples_df[['patient_tissue', 'library']].isnull().values.any():
        raise ValueError(f"id column must follow the format '{{Patient}}_{{tissue}}_{{lane or plate}}_{{library}}' for platform {platform}")
    
    # Create sample identifier
    samples_df['sample_id'] = samples_df['patient_tissue']
    
    # Check if sample names contain "."
    if samples_df['sample_id'].str.contains("\\.").any():
        raise ValueError("Sample names must not contain '.', please remove '.'")      

    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    samples_df['fastqs_dir'] = samples_df['fq1'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # Set index
    if platform == 'lane':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "lane", "library"])
    elif platform == 'plate':
        samples_df = samples_df.set_index(["sample_id","patient", "cell","tissue", "plate", "library"])
    
    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    samples_df['fastqs_dir'] = samples_df['bam'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # Set index
    if platform == 'tenX':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "lane", "library"])
    elif platform == 'smartseq':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "plate", "library"])

    for _, row in samples_df.iterrows():
        bam_exists = os.path.isfile(row['bam'])
        if not fq1_exists:
            raise FileNotFoundError(f"File not found: {row['bam']}")

    return samples_df



# def get_starsolo_sample_id(SAMPLES, wildcards, fq_column):
#     sample_id = wildcards.sample
#     try:
#         file_path = SAMPLES.loc[sample_id, fq_column]
#     except KeyError:
#         raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES DataFrame.")
#     return ','.join(sorted(file_path))

# def get_starsolo_sample_id(SAMPLES, wildcards, fq_column):
#     sample_id = wildcards.sample
#     try:
#         file_paths = SAMPLES.loc[sample_id, fq_column]
#     except KeyError:
#         raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES DataFrame.")
    
#     # If there is only one lane, the file_paths will be a string instead of a list
#     # In that case, convert it to a list
#     if isinstance(file_paths, str):
#         file_paths = [file_paths]
    
#     return ','.join(sorted(file_paths))
def get_starsolo_sample_id(SAMPLES, wildcards, fq_column):
    sample_id = wildcards.sample
    try:
        # file_paths = SAMPLES.loc[sample_id, fq_column]
        # sorted_paths = sorted(file_paths)
        # joined_paths = ','.join(sorted_paths)
        file_paths = SAMPLES.loc[sample_id, fq_column]
        sorted_paths = sorted(file_paths)
        joined_paths = ','.join(sorted_paths)
        return joined_paths
    except KeyError:
        raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES DataFrame.")
# def get_starsolo_sample_id(samples, wildcards, read):
#     sample_reads = [s[read] for s in samples if s['sample'] == wildcards.sample]
#     if read == "fq1":
#         return ','.join(sorted(sample_reads))
#     elif read == "fq2":
#         return ' '.join(sorted(sample_reads))


    

def get_sample_id(SAMPLES, wildcards, col):
    sample_id = wildcards.sample
    try:
        col_value = SAMPLES.loc[sample_id, col]
    except KeyError:
        raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES DataFrame.")
    return col_value

def get_fastqs_dir(SAMPLES, wildcards):
    """
    Get fastq dir belonging to a specific sample.
    """
    sample_id = wildcards.sample

    try:
        fastqs_dir = SAMPLES.loc[sample_id,"fastqs_dir"]
        fastqs_dir = fastqs_dir.unique()
    except KeyError:
        raise ValueError(f"Sample ID '{sample_id}' not found fastqs_dir in SAMPLES DataFrame.")
    return fastqs_dir

def get_samples_bax(wildcards,bam_dir, suffix="bam"):

    sample_id = wildcards.sample

    bam_file = os.path.join(
        bam_dir,
        "unmapped_host",
        sample_id,
        f"Aligned_sortedByCoord_unmapped_out.{suffix}")

    return bam_file


def get_samples_id_by_tissue(sample_df, tissue):
    """
    Get unique sample IDs belonging to a specific tissue.
    """
    return sample_df.loc[:, tissue, :, :].index.get_level_values("sample_id").unique()

def get_samples_id_by_patient(sample_df, patient):
    """
    Get unique sample IDs belonging to a specific patient.
    """
    return sample_df.loc[patient, :, :, :].index.get_level_values("sample_id").unique()

def get_samples_id_by_lane(sample_df, lane):
    """
    Get unique sample IDs belonging to a specific lane.
    """
    return sample_df.loc[:, :, lane, :].index.get_level_values("sample_id").unique()

def get_samples_id_by_library(sample_df, library):
    """
    Get unique sample IDs belonging to a specific library.
    """
    return sample_df.loc[:, :, :, library].index.get_level_values("sample_id").unique()


def get_tissue_by_patient(sample_df, patient):
    """
    Get unique tissues associated with a specific patient.
    """
    return sample_df.loc[patient, :, :, :].index.get_level_values("tissue").unique()

# def get_samples_id_by_plate(sample_df, plate):
#     """Get unique sample IDs belonging to a specific plate."""
#     return sample_df.loc[:, :, plate, :].index.get_level_values("sample").unique()


def get_SAMattrRGline_from_manifest(manifest_file):
    # read the manifest file
    df = pd.read_csv(manifest_file, sep="\t", header=None)
    # get the third column which contains the cell ids
    cell_ids = df[2]
    # generate the --outSAMattrRGline input format
    rgline = " , ".join([f"ID:{cell_id}" for cell_id in cell_ids])
    return rgline

def get_SAMattrRGline_by_sample(samples_df, wildcards):
    sample_id = wildcards.sample

    # Filter samples based on sample_id and plate
    samples_filtered = samples_df.loc[(samples_df.index.get_level_values("sample_id") == sample_id)]

    # samples_filtered = samples_df.loc[(samples_df.index.get_level_values("sample_id") == sample_id) &
    #                                 (samples_df.index.get_level_values("plate") == plate)]

    # Get the cell IDs for the filtered samples
    cell_ids = samples_filtered.index.get_level_values("cell").unique()

    # Generate the --outSAMattrRGline input format
    rgline = " , ".join([f"ID:{cell_id}" for cell_id in cell_ids])
    return rgline
