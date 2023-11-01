#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score
import numpy



data = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/NonRedundantACRs_18Cells.500bp_DanQ/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/NonRedundantACRs_18Cells.500bp_DanQ/test_targets.npz')
unique, counts = numpy.unique(target['data'], return_counts=True)
dict(zip(unique, counts))

Table = np.loadtxt("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/NonRedundantACRs_18Cells.500bp_DanQ/test_performance.txt", skiprows=1, usecols=(1, 2))
# Calculate the averages for roc_auc and average_precision columns
average_roc_auc = np.mean(Table[:, 0])
average_average_precision = np.mean(Table[:, 1])


for i in range(6,18):
    data = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR'+str(i)+'CT_DanQ/test_predictions.npz')
    target = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR'+str(i)+'CT_DanQ/test_targets.npz')
    unique, counts = numpy.unique(target['data'], return_counts=True)
    #dict(zip(unique, counts))
    Table = np.loadtxt("/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR"+str(i)+"CT_DanQ/test_performance.txt", skiprows=1, usecols=(1, 2))
    # Calculate the averages for roc_auc and average_precision columns
    average_roc_auc = np.mean(Table[:, 0])
    average_average_precision = np.mean(Table[:, 1])
    print(str(i)+'\t'+str(dict(zip(unique, counts))[0])+"\t"+str(dict(zip(unique, counts))[1])+"\t"+str(average_roc_auc)+"\t"+str(average_average_precision))
