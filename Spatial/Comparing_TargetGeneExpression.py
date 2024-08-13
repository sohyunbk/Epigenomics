import scanpy as sc
import squidpy as sq
import pandas as pd
import os
import argparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import math
import argparse
import os
import math
import matplotlib as mpl
import scipy.sparse

#### Designed for Pannel D - 2 WT and 2 bif3 together
#### Cluster 6 - WT CZ / Cluster 5 - Bif3 CZ

Adata = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/A619_D/adata_processed.h5ad")
Adata
print(Adata.obs.head())
print(Adata.obs['total_counts'])

cell_gene_count_matrix = Adata.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_df = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata.obs.index, columns=Adata.var.index)
print(cell_gene_count_matrix)
print(cell_gene_count_df)
