from numpy import load
import numpy

data = load('test_predictions.npz')
target = load('test_targets.npz')
unique, counts = numpy.unique(target['data'], return_counts=True)
dict(zip(unique, counts))
#{0.0: 117394, 1.0: 2606}

## This is too slow
TP = 0
TN = 0
FP = 0
FN = 0
for i in range(0,len(data['data'])):
    if data['data'][i][0] <0.5 and target['data'][i][0] ==0:
        TN+=1
    elif data['data'][i][0] >= 0.5 and target['data'][i][0] ==1:
        TP +=1
    elif data['data'][i][0] <0.5 and target['data'][i][0] ==1:
        FN +=1
    elif data['data'][i][0] >= 0.5 and target['data'][i][0] ==0:
        FP +=1
print(TP)
print(TN)
print(FP)
print(FN)

#data['data'][0][0] = 0.0005001723
## Check which sample has "1" value
for i in range(0,len(target['data'])):
    for s in target['data'][i]:
        if s == 1:
            print(target['data'][i])

NegativeSample = 0
PositiveSample = 0

print(NegativeSample)
print(PositiveSample)

##################
##################
### using sklearn

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np

data = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome/test_predictions.npz')
target = load('/scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test/500bp_AllGenome/test_targets.npz')

data = load('test_predictions.npz')
target = load('test_targets.npz')

cutoff = 0.5
y_pred_classes = np.zeros_like(data['data'])
y_pred_classes[data['data'] > cutoff] = 1

## This works only for one target
confusion_matrix(target['data'], y_pred_classes)

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

TargetArray=np.zeros(shape=len(data['data']))
PredictionArray=np.zeros(shape=len(data['data']))

for i in range(0,len(data['data'])):
    TargetArray[i] = target['data'][i][7]
    PredictionArray[i] = data['data'][i][7]

from sklearn.metrics import average_precision_score
average_precision_score(TargetArray, PredictionArray)
fpr, tpr, thresholds = metrics.roc_curve(TargetArray, PredictionArray, pos_label=1)

roc_auc_score(TargetArray, PredictionArray)
