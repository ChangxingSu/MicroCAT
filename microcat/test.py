
import pandas as pd
import re
import glob
import os
import sys
import json
from ruamel.yaml import YAML
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
    samples_df['fastqs_dir'] = samples_df['bam'].apply(lambda x: '/'.join(x.split('/')[:-1]))
    
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
        if not bam_exists:
            raise FileNotFoundError(f"File not found: {row['bam']}")

    return samples_df

def get_starsolo_sample_id(SAMPLES, wildcards, fq_column):
    sample_id = wildcards
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

SAMPLES = parse_bam_samples("/data/scRNA_analysis/benchmark/Haber2017/sample_infected.tsv",platform = "lane")

bam_reads = get_starsolo_sample_id(SAMPLES, "SH_Hpoly_Day10_Rep2_S1", "bam")

print(bam_reads)