#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/withoutCUDA.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/withoutCUDA.%j.err    # Standard error log
#SBATCH --array=0-17

SampleName=(NonRedundantACRs_18Cells.500bp
  Seedling_18Celltypes.500.RestrictACR2CT
  Seedling_18Celltypes.500.RestrictACR3CT
  Seedling_18Celltypes.500.RestrictACR4CT
  Seedling_18Celltypes.500.RestrictACR5CT
  Seedling_18Celltypes.500.RestrictACR6CT
  Seedling_18Celltypes.500.RestrictACR7CT
  Seedling_18Celltypes.500.RestrictACR8CT
  Seedling_18Celltypes.500.RestrictACR9CT
  Seedling_18Celltypes.500.RestrictACR10CT
  Seedling_18Celltypes.500.RestrictACR11CT
  Seedling_18Celltypes.500.RestrictACR12CT
  Seedling_18Celltypes.500.RestrictACR13CT
  Seedling_18Celltypes.500.RestrictACR14CT
  Seedling_18Celltypes.500.RestrictACR15CT
  Seedling_18Celltypes.500.RestrictACR16CT
  Seedling_18Celltypes.500.RestrictACR17CT)

 /home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Run_Selene_MakeYMLFile.py \
  --wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Standard_DanQ_WithoutCuda_SeqLength500bp.yml \
 --learningRate 0.0005 \
 --bedfile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/"${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz \
 --featurefile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/"${SampleName[SLURM_ARRAY_TASK_ID]}"_distinctfeatures.txt \
 --OutwmlfileName /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${SampleName[SLURM_ARRAY_TASK_ID]}".wml \
 --NewOutputDir /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${SampleName[SLURM_ARRAY_TASK_ID]}"_DanQ
