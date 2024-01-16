#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score
import numpy
import pandas as pd

WD = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"

####################################################################################################################################
## 1) Load Files
CT = [
    "L1_leaf_primordia_boundary",
    "abaxial_bundle_sheath",
    "adaxial_leaf_primordia",
    "bundle_sheath",
    "cortex",
    "dividing_leaf_primordia",
    "ground_meristem",
    "guard_mother_cell",
    "hypodermal_sclerenchyma",
    "leaf_primordia",
    "mesophyll",
    "mesophyll_precursors",
    "phloem_SE_procambial_precursors",
    "pith_parenchyma",
    "procambial_meristem",
    "protodermal_cell",
    "protophloem_SE",
    "xylem_parenchyma"
]

### 1) Control Group:

data_control_Mutated = load(WD+'/control_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/test_predictions.npz')['data']
target_control_Mutated = load(WD+'/control_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/test_targets.npz')['data']
Sample_control_Mutated = pd.read_csv(WD+'/control_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/SampledData/validate_data.bed', \
sep='\t', header=None, names=['chromosome', 'start', 'end', 'strand', 'info'])

data_control_NonMutated = load(WD+'/control_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/test_predictions.npz')['data']
target_control_NonMutated = load(WD+'/control_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/test_targets.npz')['data']
Sample_control_NonMutated = pd.read_csv(WD+'/control_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/SampledData/validate_data.bed', \
sep='\t', header=None, names=['chromosome', 'start', 'end', 'strand', 'info'])

control_change = data_control_Mutated  - data_control_NonMutated
Control_output = pd.DataFrame( index=range(18), columns=range(2))
Control_Dic = {}
for nACR in range(0,target_control_NonMutated.shape[0]): #target_control_NonMutated.shape[0]= sampleNumber
    for nCT in range(0,target_control_NonMutated.shape[1]):
        Control_Dic.setdefault(CT[nCT],{})
        Control_Dic[CT[nCT]].setdefault(target_control_NonMutated[nACR][nCT],[])
        Control_Dic[CT[nCT]][target_control_NonMutated[nACR][nCT]].append(control_change[nACR][nCT])
## Arrange Control_Dic
for i in range(0,len(CT)):
    sCT  = CT[i]
    Control_output.iloc[i, 0] = str(np.mean(Control_Dic[sCT][0.0]))+"+-"+str(np.std(Control_Dic[sCT][0.0]))
    Control_output.iloc[i, 1] = str(np.mean(Control_Dic[sCT][1.0]))+"+-"+str(np.std(Control_Dic[sCT][1.0]))
    Control_output.rename(index={i: sCT}, inplace=True)

Control_output.to_csv('/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/3.SNPDataSummary/Control_output.csv', index=True)
