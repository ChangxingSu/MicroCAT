import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_path", required=True, help="file containing the starsolo data")
parser.add_argument("--output_cluster", required=True, help="path to the output hdf5 file")
args = parser.parse_args()

adata = sc.read_10x_mtx(args.input_path,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
adata.obs[['leiden']].to_csv(args.output_cluster, sep="\t")