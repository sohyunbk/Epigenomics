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

### 2) Test Group:


data_Test_Mutated = load(WD+'/test_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/test_predictions.npz')['data']
target_Test_Mutated = load(WD+'/test_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/test_targets.npz')['data']
Sample_Test_Mutated = pd.read_csv(WD+'/test_SNVs_curated_RandomSelectSNPperACR_Mutated_500bp_DanQ/SampledData/validate_data.bed', \
sep='\t', header=None, names=['chromosome', 'start', 'end', 'strand', 'info'])

data_Test_NonMutated = load(WD+'/test_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/test_predictions.npz')['data']
target_Test_NonMutated = load(WD+'/test_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/test_targets.npz')['data']
Sample_Test_NonMutated = pd.read_csv(WD+'/test_SNVs_curated_RandomSelectSNPperACR_NotMutated_500bp_DanQ/SampledData/validate_data.bed', \
sep='\t', header=None, names=['chromosome', 'start', 'end', 'strand', 'info'])

Test_change = data_Test_Mutated  - data_Test_NonMutated
Test_output = pd.DataFrame( index=range(18), columns=range(2))
Test_Dic = {}
for nACR in range(0,target_Test_NonMutated.shape[0]): #target_control_NonMutated.shape[0]= sampleNumber
    for nCT in range(0,target_Test_NonMutated.shape[1]):
        Test_Dic.setdefault(CT[nCT],{})
        Test_Dic[CT[nCT]].setdefault(target_Test_NonMutated[nACR][nCT],[])
        Test_Dic[CT[nCT]][target_Test_NonMutated[nACR][nCT]].append(Test_change[nACR][nCT])
## Arrange Control_Dic
for i in range(0,len(CT)):
    sCT  = CT[i]
    Test_output.iloc[i, 0] = str(np.mean(Test_Dic[sCT][0.0]))+"+-"+str(np.std(Test_Dic[sCT][0.0]))
    Test_output.iloc[i, 1] = str(np.mean(Test_Dic[sCT][1.0]))+"+-"+str(np.std(Test_Dic[sCT][1.0]))
    Test_output.rename(index={i: sCT}, inplace=True)


def parse_midpoint(acrID):
    # Extracts the numeric parts and computes the midpoint
    parts = acrID.split('_')
    start = int(parts[1])
    end = int(parts[2])
    return (start + end) / 2

SNPData_test =  "/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/test_SNVs_curated_RandomSelectSNPperACR.txt"
SNPData_test_df = pd.read_csv(SNPData_test, sep='\t', header=0)

for nACR in range(0,target_control_NonMutated.shape[0]): #target_control_NonMutated.shape[0]= sampleNumber
    for nCT in range(0,target_control_NonMutated.shape[1]):
        target_acrID = str(Sample_Test_Mutated.iloc[nACR]['chromosome'])+ \
        "_"+str(Sample_Test_Mutated.iloc[nACR]['start'])+ \
        "_"+str(Sample_Test_Mutated.iloc[nACR]['end'])
        target_chromosome = target_acrID.split('_')[0]

        # Compute the midpoint of the target acrID
        target_midpoint = parse_midpoint(target_acrID)

        # Filter the DataFrame for the same chromosome
        same_chromosome_df = SNPData_test_df[SNPData_test_df['acrID'].str.startswith(target_chromosome)]
        # Calculate the absolute difference in midpoints for the same chromosome
        same_chromosome_df['midpoint'] = same_chromosome_df['acrID'].apply(parse_midpoint)
        same_chromosome_df['difference'] = abs(same_chromosome_df['midpoint'] - target_midpoint)

        # Find the acrID with the smallest difference on the same chromosome
        closest_acrID = same_chromosome_df.loc[same_chromosome_df['difference'].idxmin(), 'acrID']
        row_with_specified_acrID = SNPData_test_df[SNPData_test_df['acrID'] == closest_acrID]
        row = row_with_specified_acrID.iloc[0]
        columns_with_1 = [col for col in row.index if row[col] == 1]
