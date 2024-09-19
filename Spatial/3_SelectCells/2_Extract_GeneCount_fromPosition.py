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



def WriteTheCountData(InputDirName,BarcodeList,OutputName):
  Adata_raw = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/"+InputDirName+"/adata_raw.h5ad")
  ###################################
  cell_gene_count_matrix = Adata_raw.X
  cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
  cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata_raw.obs.index, columns=Adata_raw.var.index)
    #print(cell_gene_count_raw)
    #### Get the CZ cells from the table ###
  CZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(BarcodeList)]
  CZ_Count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/'+OutputName+'.csv', index=True, header=True)


A619_B= [
    "GCATTGTAATTCATAT-1",
    "AGCTTGATCTTAACTT-1",
    "CGCCTGGCCTACGTAA-1",
    "AAACGAGACGGTTGAT-1",
    "AAGAGGCATGGATCGC-1",
    "ATAACGCCGGAGGGTC-1",
    "TTGAATATGGACTTTC-1"
]

A619_C = [
    "TCGTATTACCCATTGC-1",
    "GGCGCAGGACATCTTC-1",
    "ATCCTGAATCGCTGCG-1",
    "ATAAATATTAGCAGCT-1",
    "CATGATGGAAGTTAGC-1",
    "AGTGAACAAACTTCTC-1",
    "AAGTGCCTTGACTGTA-1",
    "CCCAGTAAACTTGGGA-1",
    "AGAAGGTACACTTCAC-1",
    "GTAACATCTAAGATAA-1",
    "CCGACAATAGGCCGCC-1"
]

A619_D = [
    "AGGTACGATATTGCCA-1",
    "TTAAGCCGACAACTTC-1",
    "CCATAAACAACCCGAC-1",
    "ACTCGTCAGTAATCCC-1",
    "GGTGATAAGGAGCAGT-1",
    "CTTGTTGCTGAGTCAA-1"
]

### WT - B : Cluster 5
WriteTheCountData("A619_B",A619_B,"A619_B_MannuallySelected")
## WT-C :Cluster 1 + 7
WriteTheCountData("A619_C",A619_C,"A619_C_MannuallySelected")
## WT-D: WT and Bif3 Mixedd
WriteTheCountData("A619_D_Res2",A619_D,"A619_D_WT_MannuallySelected")

