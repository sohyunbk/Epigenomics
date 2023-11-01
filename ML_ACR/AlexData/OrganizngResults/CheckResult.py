#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score
import numpy

data = load('/scratch/sb14489/8.ML_ACR/2.MaizeEar/2.Selene/DanQ/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.MaizeEar/2.Selene/DanQ/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_withBigN/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome_withBigN/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_WithNegative/test_targets.npz')

target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/500bp_MappableRegions_DanQ_withoutCuda_SameNumberNegative/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_18Classes/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ/DanQ_18Classes/test_targets.npz')

data = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/NonRedundantACRs_18Cells.500bp_DanQ/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/NonRedundantACRs_18Cells.500bp_DanQ/test_targets.npz')
#{0.0: 139701, 1.0: 399435}
#>>> confusion_matrix(TargetArray, y_pred_classes)
## cut off 0.5
#array([[ 2829,  5238],
#       [ 2139, 19746]])

#Cut off 0.7
#array([[ 5583,  2484],
#       [ 7971, 13914]])

data = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR16CT_DanQ/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR16CT_DanQ/test_targets.npz')
# Cut off 0.7
# {0.0: 235953, 1.0: 303183}
#array([[10715,  2915],
#       [10714,  5608]])

# Cut off 0.5
#>>> confusion_matrix(TargetArray, y_pred_classes)
#array([[ 6495,  7135],
#       [ 4781, 11541]])

unique, counts = numpy.unique(target['data'], return_counts=True)
dict(zip(unique, counts))

#TargetArray=np.zeros(shape=len(data['data']))
#PredictionArray=np.zeros(shape=len(data['data']))

#for i in range(0,len(data['data'])):
#    TargetArray[i] = target['data'][i][7]
#    PredictionArray[i] = data['data'][i][7]

TargetArray = target['data'][:,2]
PredictionArray = data['data'][:,2]

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
