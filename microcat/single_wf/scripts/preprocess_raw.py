import pandas as pd
import warnings
import os
import json
import subprocess
import shutil
import sys

try:
    samples_df = pd.read_csv(snakemake.input[0], sep="\t")
except FileNotFoundError:
    warnings.warn(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)


# Create an empty list to store processing information
processing_info = []

if not set(['id', 'fq1', 'fq2']).issubset(samples_df.columns):
    raise ValueError("Columns 'id', 'fq1', 'fq2' must exist in the sample.tsv")

# Extract library, lane, and plate from id
samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
samples_df['is_lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1].startswith("L"))
samples_df.loc[samples_df['is_lane'], 'lane'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: x.split('_')[-1])
samples_df['patient_tissue'] = samples_df['patient_tissue_lane_plate'].apply(lambda x: '_'.join(x.split('_')[:-1]))
samples_df = samples_df.drop(columns=['patient_tissue_lane_plate'])
samples_df = samples_df.loc[(samples_df["patient_tissue"] == snakemake.wildcards["sample"])]
samples_df = samples_df.reset_index()

# # Extract required columns from the parsed samples DataFrame
# manifest_df = samples_df[['fq1', 'fq2', 'cell']]

# # Filter out rows where fq2 is NaN
# manifest_df = manifest_df[manifest_df['fq2'].notna()]
output = str(snakemake.output[0])
outdir = os.path.dirname(output)

# Remove directory if exists
if os.path.exists(outdir):
    shutil.rmtree(outdir)

# Create output directory
os.makedirs(outdir)

# Iterate through samples and create symlinks with renamed files based on patient_tissue
for index, sample in samples_df.iterrows():
    fq1 = sample['fq1']
    fq2 = sample['fq2']
    patient_tissue = sample['patient_tissue']
    lane = sample['lane']
    library = sample['library']

    # Check file extension
    if (fq1.endswith(".fastq") and fq2.endswith(".fastq")) or (fq1.endswith(".fq") and fq2.endswith(".fq")):
        file_type = 'fastq'
        # Build new filenames
        new_fq1 = os.path.join(outdir, f'{patient_tissue}_S1_{lane}_R1_{library}.fastq')
        new_fq2 = os.path.join(outdir, f'{patient_tissue}_S1_{lane}_R2_{library}.fastq')

        # Create symlinks
        subprocess.call(['ln', '-s', fq1, new_fq1])
        subprocess.call(['ln', '-s', fq2, new_fq2])
        
        # Build processing info dictionary
        processing_info.append({
            'file_type': file_type,
            'Lane': lane,
            'Library': library,
            'Original_fq1': fq1,
            'New_fq1': new_fq1,
            'Original_fq2': fq2,
            'New_fq2': new_fq2
        })
    elif (fq1.endswith(".fastq.gz") and fq2.endswith(".fastq.gz")) or (fq1.endswith(".fq.gz") and fq2.endswith(".fq.gz")):
        file_type = 'fastq.gz'
        # Build new filenames
        new_fq1 = os.path.join(outdir, f'{patient_tissue}_S1_{lane}_R1_{library}.fastq.gz')
        new_fq2 = os.path.join(outdir, f'{patient_tissue}_S1_{lane}_R2_{library}.fastq.gz')

        # Create symlinks
        subprocess.call(['ln', '-s', fq1, new_fq1])
        subprocess.call(['ln', '-s', fq2, new_fq2])

        # Build processing info dictionary
        processing_info.append({
            'file_type': file_type,
            'Lane': lane,
            'Library': library,
            'Original_fq1': fq1,
            'New_fq1': new_fq1,
            'Original_fq2': fq2,
            'New_fq2': new_fq2
        })
    else:
        raise ValueError("The file extension must be either .fastq or .fastq.gz")

# Write processing information to JSON file
with open(snakemake.output[0], 'w') as json_file:
    json.dump(processing_info, json_file, indent=4)