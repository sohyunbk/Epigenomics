#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score
import numpy

data_control_Mutated = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/control_SNVs_curated_RandomSelectSNPperACR_Mutated_DanQ/test_predictions.npz')
target_control_Mutated = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/control_SNVs_curated_RandomSelectSNPperACR_Mutated_DanQ/test_targets.npz')

data_control_NonMutated = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/control_SNVs_curated_RandomSelectSNPperACR_NotMutated_DanQ/test_predictions.npz')
target_control_NonMutated = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/control_SNVs_curated_RandomSelectSNPperACR_NotMutated_DanQ/test_targets.npz')

print(list(data_control_NonMutated.keys()))
data_control_NonMutated['data']
len(data_control_NonMutated['data'])

unique, counts = numpy.unique(target_control_Mutated['data'], return_counts=True)
dict(zip(unique, counts))

unique, counts = numpy.unique(target_control_NonMutated['data'], return_counts=True)
dict(zip(unique, counts))

PredictionArray = data_control_Mutated['data'][:,2]
PredictionArray2 = data_control_NonMutated['data'][:,2]

#TargetArray = target['data'][:,14]
#PredictionArray = data['data'][:,14]

from sklearn.metrics import average_precision_score
from sklearn import metrics
average_precision_score(TargetArray, PredictionArray)
fpr, tpr, thresholds = metrics.roc_curve(TargetArray, PredictionArray, pos_label=1)

roc_auc_score(TargetArray, PredictionArray)
##################
##################

from numpy import load
import numpy
from sklearn.metrics import confusion_matrix

cutoff = 0.5
y_pred_classes = np.zeros(len(PredictionArray))
y_pred_classes[PredictionArray > cutoff] = 1

np.count_nonzero(TargetArray)
len(TargetArray)-np.count_nonzero(TargetArray)
confusion_matrix(TargetArray, y_pred_classes)
count_ones = np.count_nonzero(TargetArray == 1)

# Count the number of 0s
count_zeros = np.count_nonzero(TargetArray == 0)

print("Number of 1s:", count_ones)
print("Number of 0s:", count_zeros)


TP = 0
TN = 0
FP = 0
FN = 0
for i in range(0,len(data['data'])):
    if data['data'][i][7] <=0.5 and target['data'][i][7] ==0:
        TN+=1
    elif data['data'][i][7] > 0.5 and target['data'][i][7] ==1:
        TP +=1
    elif data['data'][i][7] <=0.5 and target['data'][i][7] ==1:
        FN +=1
    elif data['data'][i][7] > 0.5 and target['data'][i][7] ==0:
        FP +=1
print(TP)
print(TN)
print(FP)
print(FN)



###### Matfile
import scipy.io
data = scipy.io.loadmat('/scratch/sb14489/8.ML_ACR/DeepFormer_Ex/DeepFormer/data/test.mat')
