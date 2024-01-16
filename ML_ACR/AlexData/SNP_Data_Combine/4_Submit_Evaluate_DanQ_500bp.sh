#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=32             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/withoutCUDA.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/withoutCUDA.%j.err    # Standard error log
#SBATCH --array=0-3

SampleName=(
control_SNVs_curated_RandomSelectSNPperACR
test_SNVs_curated_RandomSelectSNPperACR
control_SNVs_curated_RandomSelectSNPperACR
test_SNVs_curated_RandomSelectSNPperACR)

OutPutName=(
control_SNVs_curated_RandomSelectSNPperACR_NotMutated
test_SNVs_curated_RandomSelectSNPperACR_NotMutated
control_SNVs_curated_RandomSelectSNPperACR_Mutated
test_SNVs_curated_RandomSelectSNPperACR_Mutated)

FastaFile=(
/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa
/scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa
/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/ControlSNPChange_MaizeV5_RandomSNPSelection.fa
/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/TestSNPChange_MaizeV5_RandomSNPSelection.fa)

 /home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Run_Selene_MakeYMLFile_Evaluate.py \
  --wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Evaluate_DanQ_withoutCuda_SeqLength500.yml \
 --learningRate 0.0005 \
 --bedfile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz \
 --featurefile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${SampleName[SLURM_ARRAY_TASK_ID]}"_distinctfeatures.txt \
 --OutwmlfileName /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${OutPutName[SLURM_ARRAY_TASK_ID]}"_500bp.wml \
 --NewOutputDir /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${OutPutName[SLURM_ARRAY_TASK_ID]}"_500bp_DanQ \
 --fasta "${FastaFile[SLURM_ARRAY_TASK_ID]}" \
 --TrainModelFile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/Seedling_18Celltypes.500.RestrictACR11CT_DanQ_500bp/best_model.pth.tar
