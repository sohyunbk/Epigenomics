import scanpy as sc
import squidpy as sq
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_path', required=True)
args = parser.parse_args()

# Create the output directory if it does not exist
adata = sq.read.visium(args.input_path)

pd.set_option('display.max_columns', None)

adata.var_names_make_unique()
adata.var.set_index('gene_ids', inplace=True)
adata.var.index.name = 'gene_id'

adata.var['Mt'] = adata.var_names.str.startswith('chrMt')
adata.var['PT'] = adata.var_names.str.startswith('chrPt')
sc.pp.calculate_qc_metrics(adata, qc_vars=['Mt', 'PT'], inplace=True)

print("Current working directory:", os.getcwd())
adata.write("./adata.h5ad")
