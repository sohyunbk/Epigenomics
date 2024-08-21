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


def WriteTheCountData(InputDirName,ClusterName,OutputName):
  Adata_raw = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/"+InputDirName+"/adata_raw.h5ad")
  Adata = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/"+InputDirName+"/adata_processed.h5ad")
  #Adata
  Meta = Adata.obs
  CZ = Meta[Meta['clusters'] == ClusterName]
  #print(CZ)
  ###################################
  cell_gene_count_matrix = Adata_raw.X
  cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
  cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata_raw.obs.index, columns=Adata_raw.var.index)
  #print(cell_gene_count_raw)
  #### Get the CZ cells from the table ###
  barcode_set = CZ.index.tolist()
  CZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(barcode_set)]
  CZ_Count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene/'+OutputName+'.csv', index=True, header=True)

### WT - B : Cluster 5
WriteTheCountData("A619_B","5","A619_B_Cluser5")
## WT-C :Cluster 1 + 7
WriteTheCountData("A619_C","1","A619_C_Cluser1")
WriteTheCountData("A619_C","7","A619_C_Cluser7")


############################################
## WT-D : WT : special! 
InputDirName = "A619_D_Res2"
Adata_raw = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/"+InputDirName+"/adata_raw.h5ad")
Adata = sc.read("/scratch/sb14489/9.spatialRNAseq/2.QC_MarkerGene/"+InputDirName+"/adata_processed.h5ad")
#Adata
Meta = Adata.obs
Meta
array_row_list = Meta[Meta['clusters'] == "15"]['array_row'].astype(int).tolist()
## The max : 38 is upto the WT
Meta['array_row'] = Meta['array_row'].astype(int)
Meta['Genotype'] = Meta['array_row'].apply(lambda x: 'WT' if x <= 39 else 'Bif3')

## WT CentralZone
CZ1 = Meta[(Meta['clusters'] == "15") & (Meta['Genotype'] == "WT")]
CZ2 = Meta[(Meta['clusters'] == "1") & (Meta['Genotype'] == "WT")]
CZ3 = Meta[(Meta['clusters'] == "7") & (Meta['Genotype'] == "WT")]
CZ = pd.concat([CZ1,CZ2,CZ3], axis=0)
##################################
cell_gene_count_matrix = Adata_raw.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata_raw.obs.index, columns=Adata_raw.var.index)
#print(cell_gene_count_raw)
#### Get the CZ cells from the table ###
barcode_set = CZ.index.tolist()
CZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(barcode_set)]
CZ_Count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene/A619_D_WTCZOnly.csv', index=True, header=True)

##### Bif3 Special
#Adata
Meta = Adata.obs
Meta
array_row_list = Meta[Meta['clusters'] == "15"]['array_row'].astype(int).tolist()
## The max : 38 is upto the WT
Meta['array_row'] = Meta['array_row'].astype(int)
Meta['Genotype'] = Meta['array_row'].apply(lambda x: 'WT' if x <= 39 else 'Bif3')

## WT CentralZone
CZ1 = Meta[(Meta['clusters'] == "6") & (Meta['Genotype'] == "Bif3")]
CZ2 = Meta[(Meta['clusters'] == "7") & (Meta['Genotype'] == "Bif3")]
CZ = pd.concat([CZ1,CZ2], axis=0)
##################################
cell_gene_count_matrix = Adata_raw.X
cell_gene_count_matrix_dense = cell_gene_count_matrix.toarray() if isinstance(cell_gene_count_matrix, scipy.sparse.spmatrix) else cell_gene_count_matrix
cell_gene_count_raw = pd.DataFrame(cell_gene_count_matrix_dense, index=Adata_raw.obs.index, columns=Adata_raw.var.index)
#print(cell_gene_count_raw)
#### Get the CZ cells from the table ###
barcode_set = CZ.index.tolist()
CZ_Count = cell_gene_count_raw.loc[cell_gene_count_raw.index.isin(barcode_set)]
CZ_Count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene/A619_D_Bif3CZOnly.csv', index=True, header=True)


