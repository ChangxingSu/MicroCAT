import argparse
import pandas as pd
import sys
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--krak_output_file", action="store", help="path to kraken output file")
parser.add_argument("--krak_study_denosing_file", action="store", help="path to krak_study_denosing file")
parser.add_argument("--out_krak2_denosing", action="store", help="output path to save files")
args = parser.parse_args()

# Read krak2 output file and create a copy
krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type', 'query_name', 'taxid_info', 'len', 'kmer_position'])
krak2_output_copy = krak2_output.copy()

# Extract 'taxa' and 'taxid' from 'taxid_info' column
krak2_output_copy[['taxa', 'taxid']] = krak2_output_copy['taxid_info'].str.extract(r'(.*) \(taxid (\d+)\)')
krak2_output_copy['taxid'] = krak2_output_copy['taxid'].str.replace(r'\)', '').str.strip()

# Read taxa file (krak_study_denosing)
krak_study_denosing = pd.read_csv(args.krak_study_denosing_file, sep="\t")

krak2_output_copy['taxid'] =krak2_output_copy['taxid'].astype(str)
krak_study_denosing['ncbi_taxa'] = krak_study_denosing['ncbi_taxa'].astype(str)
# Filter krak2_output_copy to keep only rows with taxid appearing in krak_study_denosing ncbi_taxa
krak2_output_filtered = krak2_output_copy[krak2_output_copy['taxid'].isin(krak_study_denosing['ncbi_taxa'])]

# Filter the original krak2_output based on the filtered krak2_output_copy
krak2_output_to_save = krak2_output[krak2_output['query_name'].isin(krak2_output_filtered['query_name'])]

# Save the filtered krak2_output to the specified output file
krak2_output_to_save.to_csv(args.out_krak2_denosing, sep='\t', index=False, header=False)
