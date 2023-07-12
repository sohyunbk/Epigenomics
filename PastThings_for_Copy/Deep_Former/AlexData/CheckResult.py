#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score

data = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_RemoveRedundantACR_Try500bp/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_RemoveRedundantACR_Try500bp/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_withBigN/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_withBigN/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_SameNumberNegative/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_SameNumberNegative/test_targets.npz')

unique, counts = numpy.unique(target['data'], return_counts=True)
dict(zip(unique, counts))

#TargetArray=np.zeros(shape=len(data['data']))
#PredictionArray=np.zeros(shape=len(data['data']))

#for i in range(0,len(data['data'])):
#    TargetArray[i] = target['data'][i][7]
#    PredictionArray[i] = data['data'][i][7]

TargetArray = target['data'][:,0]
PredictionArray = data['data'][:,0]

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

cutoff = 0.8
y_pred_classes = np.zeros(len(PredictionArray))
y_pred_classes[PredictionArray > cutoff] = 1

np.count_nonzero(TargetArray)
len(TargetArray)-np.count_nonzero(TargetArray)
confusion_matrix(TargetArray, y_pred_classes)



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
