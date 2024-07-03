import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--input_file', required=True, help='Path to the input AnnData file')
parser.add_argument('--output_file', required=True, help='Path to save the output AnnData file')
parser.add_argument('--output_name', required=True, help='Name for the output plots')
args = parser.parse_args()

# Read the input file
adata = sc.read(args.input_file)

# Normalize and preprocess the data
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# Manifold embedding and clustering based on transcriptional similarity
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)

# Save UMAP plot
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4, save="_" + args.output_name)

# Write the output file
adata.write(args.output_file)
