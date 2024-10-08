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
  CZ_Count.to_csv('/scratch/sb14489/9.spatialRNAseq/3.TargetGene_IMOnly/'+OutputName+'.csv', index=True, header=True)

### WT - B : Cluster 5
WriteTheCountData("A619_B","5","A619_B_Cluser5")
## WT-C :Cluster 1 + 7
WriteTheCountData("A619_C","7","A619_C_Cluser7")
## WT-D: WT and Bif3 Mixedd
WriteTheCountData("A619_D_Res2","15","A619_D_WT_Cluser15")
WriteTheCountData("A619_D_Res2","6","A619_D_Bif3_Cluser6")

### Bif3 - A : Cluster 1 and 10
WriteTheCountData("Bif3_A_HighDepth_Res1.5","1","Bif3_A_Cluser1")
### Bif3 - B : Cluster 6 and 5
WriteTheCountData("Bif3_B_HighDepth_Res1.5","6","Bif3_B_Cluser6")
### Bif3 - C :*** its not Res1.5
WriteTheCountData("Bif3_C_HighDepth","11","Bif3_C_Cluser11")
### Bif3 - D : Cluster *** its not Res1.5
WriteTheCountData("Bif3_D_HighDepth","7","Bif3_D_Cluser7")
