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
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

### WT 
Adata_raw = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/A619_D/adata_raw.h5ad")
Adata = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/A619_D/adata_processed.h5ad")
Adata
Meta = Adata.obs

####### Getting WT CZ cells we are interested in 
print(Meta.head())
print(Adata.obs['total_counts'])
WT_CZ = Meta[Meta['clusters'] == '12']
print(WT_CZ)
### Only use the right WT sample : 5 cells
### Barcode ID : ACTCGTCAGTAATCCC-1, CGAGGGACTGCGGTCG-1, CTTGTTGCTGAGTCAA-1, GCGATTGTTAACGTTA-1, GGTGATAAGGAGCAGT-1
WT_CZ_Barcodeslist = ["ACTCGTCAGTAATCCC-1", "CGAGGGACTGCGGTCG-1", "CTTGTTGCTGAGTCAA-1", "GCGATTGTTAACGTTA-1", "GGTGATAAGGAGCAGT-1"]

#WT_CZ_Barcodeslist = WT_CZ.index.tolist()

####### Getting Bif3 CZ cells we are interested in 

Bif3_CZ = Meta[Meta['clusters'] == '5' ]
print(Bif3_CZ)
# Barcode ID 10 cells: AAATCTAGCCCTGCTA-1, AGAGGCTTCGGAAACC-1, AGGAGGCCTTCGCGCG-1,ATTGCGATCAGTAACT-1, 
#ATTGGATTACAGCGTA-1, CCATATGGAAACTATA-1, CGCAGGCGATCCAAAC-1, GCGAAACGATCGGGAG-1,TCCGCGGCCCAATGAA-1, TCCGTTAAGCTAATAT-1  
Bif3_CZ_Barcodeslist = [
    "AAATCTAGCCCTGCTA-1", 
    "AGAGGCTTCGGAAACC-1", 
    "AGGAGGCCTTCGCGCG-1", 
    "ATTGCGATCAGTAACT-1", 
    "ATTGGATTACAGCGTA-1", 
    "CCATATGGAAACTATA-1", 
    "CGCAGGCGATCCAAAC-1", 
    "GCGAAACGATCGGGAG-1", 
    "TCCGCGGCCCAATGAA-1", 
    "TCCGTTAAGCTAATAT-1"
]

#Bif3_CZ_Barcodeslist = Bif3_CZ.index.tolist()

###################################
cell_gene_count_matrix = Adata_raw.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata.obs.index, columns=Adata.var.index)
print(cell_gene_count_raw)

cell_gene_count_matrix = Adata.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_normalized = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata.obs.index, columns=Adata.var.index)
print(cell_gene_count_normalized)

#### Get the CZ cells from the table ###

barcode_set = set(Bif3_CZ_Barcodeslist)
Bif3CZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(barcode_set)]

barcode_set = set(WT_CZ_Barcodeslist)
WTCZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(barcode_set)]

#### TMM normalization: EdgeR in python - rpy2 is needed.
pandas2ri.activate()
edgeR = importr('edgeR')

All_count = pd.concat([WTCZ_Count,Bif3CZ_Count]).T

## it was hard to work with rpy2... should move to R.
All_count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene/WTFive_Bif3Ten_count.csv', index=True, header=True)


#### Load the target genes ###########
Markergenes = pd.read_csv('/scratch/sb14489/3.scATAC/0.Data/MarkerGene/SpatialMarkerFinal.txt', sep='\t')
#filtered_genes = Markergenes[Markergenes['name'].str.startswith('arf')]
ARF_list = [
    "arftf4", "arftf30", "arftf18", "arftf3", "arftf20",  # Activator
    "arftf25", "arftf10", "arftf36",  # Repressor
    "arftf23", "arftf26"
]
ARF_geneID = Markergenes[Markergenes['name'].isin(ARF_list)]['geneID']
#len(gene_ids)
#Zm00001eb001720 = knox1
AllPlotgenes= ARF_geneID.tolist()
AllPlotgenes.append("Zm00001eb001720")

AllPlotgenes_symbol = ARF_list
AllPlotgenes_symbol.append("knox1")

Bif3CZ_Count['Zm00001eb433460']
WTCZ_Count['Zm00001eb433460']

Bif3CZ_Count['Zm00001eb066640']
WTCZ_Count['Zm00001eb066640']

Bif3CZ_Count['Zm00001eb001720']
WTCZ_Count['Zm00001eb001720']
