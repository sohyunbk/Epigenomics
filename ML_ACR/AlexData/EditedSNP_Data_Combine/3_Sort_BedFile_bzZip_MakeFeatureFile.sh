#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/SortBed_MakeFeatureFile.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/SortBed_MakeFeatureFile.%j.err    # Standard error log
#SBATCH --array=0-1                   # Array range

SampleName=(control_SNVs.v2.curated1000bp
test_SNVs.v2.curated1000bp)

WD="/scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"

ml tabix/0.2.6-GCCcore-11.3.0

sort -k1V -k2n -k3n "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}".bed > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed

bgzip -c "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed  > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz
tabix -p bed "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz
cut -f 4 "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed | sort -u > "$WD""${SampleName[SLURM_ARRAY_TASK_ID]}"_distinctfeatures.txt
