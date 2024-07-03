import scanpy as sc
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', required=True)
args = parser.parse_args()

adata = sc.read(args.input_file)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

## 4) Manifold embedding and clustering based on transcriptional similarity
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(
    adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2
)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts","clusters"], wspace=0.4, save="_"+output_name)

adata.write("adata_norm.h5ad")
