# import function
from cellbender.remove_background.downstream import anndata_from_h5
import numpy as np
import pandas as pd
import scanpy as sc
import argparse
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_hdf5", required=True, help="file containing the cellbender filtered hdf5 file")
parser.add_argument("--output_cluster", required=True, help="path to the output hdf5 file")
args = parser.parse_args()

# load the data
adata = anndata_from_h5(args.input_hdf5)

adata = adata[adata.obs['cell_probability'] > 0.5]
# compute a UMAP and do clustering using the cellbender latent gene expression embedding
sc.pp.neighbors(adata, use_rep='gene_expression_encoding', metric='euclidean', method='umap')
sc.tl.umap(adata)
sc.tl.leiden(adata)

adata.obs[['leiden']].to_csv(args.output_cluster, sep="\t")