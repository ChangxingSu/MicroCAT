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
    
    # Check if id, fq1, fq2 columns exist
    if not set(['id', 'fq1', 'fq2']).issubset(samples_df.columns):
        raise ValueError("Columns 'id', 'fq1', 'fq2' must exist in the sample.tsv")
    
    # Extract library, lane, and plate from id
    samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
    
    # Determine the platform and parse accordingly
    if platform == 'tenX':
        samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
        samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
    elif platform == 'smartseq':
        samples_df['is_plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("P"))
        samples_df.loc[samples_df['is_plate'], 'plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
    else:
        raise ValueError("Platform must be either 'tenX' or 'smartseq'")
    
    # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
    # and is always the last part in patient_tissue
    samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    samples_df['tissue'] = samples_df['patient_tissue'].apply(lambda x: x.split('_')[-1])
    samples_df['patient'] = samples_df['patient_tissue'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])

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
    if platform == 'tenX':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "lane", "library"])
    elif platform == 'smartseq':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "plate", "library"])

    return samples_df


def parse_bam_samples(sample_tsv, platform):
    samples_df = pd.read_csv(sample_tsv, sep="\t")
    
    # Check if id, fq1, fq2 columns exist
    if not set(['id', 'bam_pe']).issubset(samples_df.columns):
        raise ValueError("Columns 'id', 'bam_pe', must exist in the sample.tsv")
    
    # Extract library, lane, and plate from id
    samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
    
    # Determine the platform and parse accordingly
    if platform == 'tenX':
        samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
        samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
    elif platform == 'smartseq':
        samples_df['is_plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("P"))
        samples_df.loc[samples_df['is_plate'], 'plate'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
    else:
        raise ValueError("Platform must be either 'tenX' or 'smartseq'")
    
    # Extract patient and tissue, using the fact that tissue is always "S" followed by a number
    # and is always the last part in patient_tissue
    samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    samples_df['tissue'] = samples_df['patient_tissue'].apply(lambda x: x.split('_')[-1])
    samples_df['patient'] = samples_df['patient_tissue'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])


    if samples_df[['patient_tissue', 'library']].isnull().values.any():
        raise ValueError(f"id column must follow the format '{{Patient}}_{{tissue}}_{{lane or plate}}_{{library}}' for platform {platform}")
    
    # Create sample identifier
    samples_df['sample_id'] = samples_df['patient_tissue']
    
    # Check if sample names contain "."
    if samples_df['sample_id'].str.contains("\\.").any():
        raise ValueError("Sample names must not contain '.', please remove '.'")
    
    # Determine if the sequencing is paired-end or single-end
    if BAM_HEADERS["SE"] in samples_df.columns:
        samples_df['seq_type'] = 'single-end'
    elif BAM_HEADERS["PE"] in samples_df.columns:
        samples_df['seq_type'] = 'paired-end'
    else:
        raise ValueError("Bam file must be either 'bam_pe' or 'bam_se'")
    
    # Create a 'fastqs_dir' column that contains the directory of the fastq files
    samples_df['fastqs_dir'] = samples_df['bam_pe'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
    # Set index
    if platform == 'tenX':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "lane", "library"])
    elif platform == 'smartseq':
        samples_df = samples_df.set_index(["sample_id","patient", "tissue", "plate", "library"])

    return samples_df



def get_starsolo_sample_id(SAMPLES, wildcards, fq_column):
    sample_id = wildcards.sample
    try:
        file_path = SAMPLES.loc[sample_id, fq_column]
    except KeyError:
        raise ValueError(f"Sample ID '{sample_id}' not found in SAMPLES DataFrame.")
    return file_path


def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]

def get_fastqs_dir(SAMPLES, wildcards):
    """
    Get fastq dir belonging to a specific sample.
    """
    sample_id = wildcards.sample

    try:
        fastqs_dir = SAMPLES.loc[sample_id,"fastqs_dir"]
    except KeyError:
        raise ValueError(f"Sample ID '{sample_id}' not found fastqs_dir in SAMPLES DataFrame.")
    return fastqs_dir


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

