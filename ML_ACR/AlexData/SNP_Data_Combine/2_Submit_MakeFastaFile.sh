#!/bin/bash
#SBATCH --job-name=MakeMutantFasta        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MakeMutantFasta.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MakeMutantFasta.%j.err    # Standard error log
#SBATCH --array=0-1


Fasta=(TestSNPChange_MaizeV5_RandomSNPSelection.fa ControlSNPChange_MaizeV5_RandomSNPSelection.fa TestSNPChange_MaizeV5.fa ControlSNPChange_MaizeV5.fa)
SNPFiles=(test_SNVs_curated_RandomSelectSNPperACR.txt control_SNVs_curated_RandomSelectSNPperACR.txt test_SNVs_curated.txt control_SNVs_curated.txt)

~/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/AlexData/SNP_Data_Combine/2_MakeMutatedFasta.py \
-FastaFile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${Fasta[SLURM_ARRAY_TASK_ID]}" \
-SNPFile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${SNPFiles[SLURM_ARRAY_TASK_ID]}"
