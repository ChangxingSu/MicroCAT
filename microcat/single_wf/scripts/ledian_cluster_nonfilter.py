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

if os.path.isdir(args.input_path):
    adata = sc.read_10x_mtx(args.input_path,  # the directory with the `.mtx` file
        var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
        cache=True)                              # write a cache file for faster subsequent reading
elif os.path.isfile(args.input_path):
    adata = sc.read_h5ad(args.input_path,  # the directory with the `h5` file
                        )
else:
    raise FileNotFoundError(f"Invalid path: {args.input_path}")

adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

if adata.shape[0] < 300:
    open(args.output_cluster, 'a').close()
else:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    adata.obs[['leiden']].to_csv(args.output_cluster, sep="\t")

# open(args.output_cluster, 'a').close()