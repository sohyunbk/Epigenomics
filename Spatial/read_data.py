import scanpy as sc
import squidpy as sq
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_path', required=True)
parser.add_argument('--output_path', required=True)
parser.add_argument('--output_name', required=True)
parser.add_argument('--MarkerGene', required=True)
args = parser.parse_args()

# Create the output directory if it does not exist
print(args.output_path)
if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)

adata = sq.read.visium(args.input_path)
output_path = args.output_path
output_name = args.output_name
MarkerGene = args.MarkerGene

pd.set_option('display.max_columns', None)
sc.settings.figdir = output_path

adata.var_names_make_unique()
adata.var.set_index('gene_ids', inplace=True)
adata.var.index.name = 'gene_id'

adata.var['Mt'] = adata.var_names.str.startswith('chrMt')
adata.var['PT'] = adata.var_names.str.startswith('chrPt')
sc.pp.calculate_qc_metrics(adata, qc_vars=['Mt', 'PT'], inplace=True)

adata.write(f"{output_path}/adata.h5ad")
