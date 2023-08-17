import pandas as pd
import warnings

try:
    samples_df = pd.read_csv(snakemake.input[0], sep="\t")
except FileNotFoundError:
    warnings.warn(f"ERROR: the samples file does not exist. Please see the README file for details. Quitting now.")
    sys.exit(1)

if not set(['id', 'fq1', 'fq2']).issubset(samples_df.columns):
    raise ValueError("Columns 'id', 'fq1', 'fq2' must exist in the sample.tsv")

# Extract library, lane, and plate from id
samples_df[['patient_tissue_lane_plate', 'library']] = samples_df['id'].str.rsplit("_", n=1, expand=True)
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

# samples_df = samples_df.loc[(samples_df.plate == snakemake.wildcards["plate"]) & (samples_df["patient_tissue"] == snakemake.wildcards["sample"])]
samples_df = samples_df.loc[(samples_df["patient_tissue"] == snakemake.wildcards["sample"])]

samples_df = samples_df.reset_index()

# Extract required columns from the parsed samples DataFrame
manifest_df = samples_df[['fq1', 'fq2', 'cell']]
# Check if fq2 is non-empty for any rows
if any(samples_df['fq2'].notna()):
    warnings.warn("WARNING: Detected samples with both single-end and paired-end sequencing data.")

# Remove rows where fq2 is not empty (not NaN)
manifest_df = manifest_df[manifest_df['fq2'].isna()]

manifest_df['fq2'] = manifest_df['fq2'].fillna('-')
manifest_df.to_csv(snakemake.output[0], sep='\t', index=False, header=False)