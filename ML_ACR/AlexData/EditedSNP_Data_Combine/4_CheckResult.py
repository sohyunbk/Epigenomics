#conda activate ARF_ML_sklearnUp

from sklearn.metrics import confusion_matrix
from numpy import load
import numpy as np
from sklearn.metrics import roc_auc_score
import numpy
import pandas as pd

WD = "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData//"

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

PredictedValue = load(WD+'/Control/test_predictions.npz')['data']
PredictedValue.shape
## Input bed?

bed = pd.read_csv(WD+'/converted_control_SNVs.v2.curated1000bp.bed', \
sep='\t', header=None)
bed[1] = bed[1] + 250  # Add 250 to the second column
bed[2] = bed[2] - 250  # Subtract 250 from the third column
bed[5] = bed[0].astype(str) + "_" + bed[1].astype(str) + "_" + bed[2].astype(str)
bed_order = bed[5].rename('acrID')

Sample = pd.read_csv(WD+'/control_SNVs.v2.curated.txt', \
sep='\t', header=0)
Sample_ordered = pd.merge(Sample, bed_order, how='right', on='acrID')
CT_modified_NonSNPChange = [ct + ".Prediction.NonSNPChange" for ct in CT]
predicted_df = pd.DataFrame(PredictedValue, columns=CT_modified_NonSNPChange)
Sample_Combined = pd.concat([Sample_ordered.reset_index(drop=True), predicted_df], axis=1)

PredictedValue_SNPChange = load(WD+'/Control_SNPChange/test_predictions.npz')['data']
CT_modified_SNPChange = [ct + ".Prediction.SNPChange" for ct in CT]
predicted_df_SNPChange = pd.DataFrame(PredictedValue_SNPChange, columns=CT_modified_SNPChange)
Sample_Combined = pd.concat([Sample_Combined.reset_index(drop=True), predicted_df_SNPChange], axis=1)

Sample_Combined.to_csv(WD+"/control_predcition.csv", index=False)


### 2) Test Group:

PredictedValue = load(WD+'/Test/test_predictions.npz')['data']
PredictedValue.shape
## Input bed?

bed = pd.read_csv(WD+'/test_SNVs.v2.curated1000bp.bed', \
sep='\t', header=None)
bed[1] = bed[1] + 250  # Add 250 to the second column
bed[2] = bed[2] - 250  # Subtract 250 from the third column
bed[5] = bed[0].astype(str) + "_" + bed[1].astype(str) + "_" + bed[2].astype(str)
bed_order = bed[5].rename('acrID')

Sample = pd.read_csv(WD+'/test_SNVs.v2.curated.txt', \
sep='\t', header=0)
Sample_ordered = pd.merge(Sample, bed_order, how='right', on='acrID')
CT_modified_NonSNPChange = [ct + ".Prediction.NonSNPChange" for ct in CT]
predicted_df = pd.DataFrame(PredictedValue, columns=CT_modified_NonSNPChange)
Sample_Combined = pd.concat([Sample_ordered.reset_index(drop=True), predicted_df], axis=1)

PredictedValue_SNPChange = load(WD+'/Test_SNPChange/test_predictions.npz')['data']
CT_modified_SNPChange = [ct + ".Prediction.SNPChange" for ct in CT]
predicted_df_SNPChange = pd.DataFrame(PredictedValue_SNPChange, columns=CT_modified_SNPChange)
Sample_Combined = pd.concat([Sample_Combined.reset_index(drop=True), predicted_df_SNPChange], axis=1)

Sample_Combined.to_csv(WD+"/test_predcition.csv", index=False)
