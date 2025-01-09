import argparse
import pandas as pd
import sys
import re
# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--krak_output_file", action="store", help="path to kraken output file")
parser.add_argument("--krak_study_denosing_file", action="store", help="path to krak_study_denosing file")
parser.add_argument("--out_krak2_denosing", action="store", help="output path to save files")
args = parser.parse_args()

# # Read krak2 output file and create a copy
# krak2_output = pd.read_csv(args.krak_output_file, sep="\t", names=['type', 'query_name', 'taxid_info', 'len', 'kmer_position'])
# krak2_output_copy = krak2_output.copy()

# Read taxa file (krak_study_denosing)
krak_study_denosing = pd.read_csv(args.krak_study_denosing_file, sep="\t")

# krak2_output_copy['taxid'] =krak2_output_copy['taxid'].astype(str)
krak_study_denosing['ncbi_taxa'] = krak_study_denosing['ncbi_taxa'].astype(str)
desired_taxid_list = set(krak_study_denosing['ncbi_taxa'].unique())

with open(args.krak_output_file, 'r') as kfile_in:
    with open(args.out_krak2_denosing, 'w') as kfile_out:
        for kraken_line in kfile_in:
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                # To handle situation like: `Blattabacterium sp. (Nauphoeta cinerea) (taxid 1316444)`
                # kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
                read_type, query_name, taxid_info, read_len, kmer_position = kraken_line.strip().split('\t')
                tax_id = str(re.search(r'\(taxid (\d+)\)', taxid_info).group(1))
            # except:
            #     # in this case, something is wrong!
            #     print("Here is an error. Queryname: {}".format(query_name))
            #     # sys.exit()
            except (ValueError, KeyError) as e:
                print("Error occur:", e)
                continue
            if tax_id in desired_taxid_list:
                kfile_out.write(kraken_line)

# # Extract 'taxa' and 'taxid' from 'taxid_info' column
# krak2_output_copy[['taxa', 'taxid']] = krak2_output_copy['taxid_info'].str.extract(r'(.*) \(taxid (\d+)\)')
# krak2_output_copy['taxid'] = krak2_output_copy['taxid'].str.replace(r'\)', '').str.strip()


# Filter krak2_output_copy to keep only rows with taxid appearing in krak_study_denosing ncbi_taxa
# krak2_output_filtered = krak2_output_copy[krak2_output_copy['taxid'].isin(krak_study_denosing['ncbi_taxa'])]

# # Filter the original krak2_output based on the filtered krak2_output_copy
# krak2_output_to_save = krak2_output[krak2_output['query_name'].isin(krak2_output_filtered['query_name'])]

# Save the filtered krak2_output to the specified output file
krak2_output_to_save.to_csv(args.out_krak2_denosing, sep='\t', index=False, header=False)
