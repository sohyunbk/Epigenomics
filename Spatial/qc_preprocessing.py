import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--output_name', required=True)
args = parser.parse_args()
print("Current working directory:", os.getcwd())
# Create the output directory if it does not exist
adata = sc.read("adata.h5ad")

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

plt.savefig(args.output_name+"_QC_Histogram.pdf") ## Save Figure

print(f'Before filtering:\n cell - {adata.n_obs}; gene - {adata.n_vars}')       # check how many genes X cells


sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_cells(adata, min_genes=50)
adata = adata[adata.obs["total_counts_Mt"] < 20].copy()
adata = adata[adata.obs["total_counts_PT"] < 20].copy()

print(f"#cells after MT filter: {adata.n_obs}")

plt.savefig(f"{args.output_name}_QC_Histogram.pdf")
adata.write("adata_qc.h5ad")
