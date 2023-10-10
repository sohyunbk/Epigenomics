#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/SortBed_MakeFeatureFile.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/SortBed_MakeFeatureFile.%j.err    # Standard error log
#SBATCH --array=0-16                   # Array range

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

WD="/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/"

module load tabix/0.2.6-GCCcore-8.3.0

sort -k1V -k2n -k3n "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}".bed > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed

bgzip -c "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed  > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz
tabix -p bed "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz
cut -f 4 "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed | sort -u > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_distinctfeatures.txt
