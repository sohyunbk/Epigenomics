#!/bin/bash
#SBATCH --job-name=MakeMutantFasta        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MakeMutantFasta.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MakeMutantFasta.%j.err    # Standard error log
#SBATCH --array=0-1


Fasta=(test_SNVs.v2.curated.fa control_SNVs.v2.curated.fa)
SNPFiles=(test_SNVs.v2.curated.txt control_SNVs.v2.curated.txt)

~/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/AlexData/SNP_Data_Combine/2_MakeMutatedFasta.py \
-FastaFile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${Fasta[SLURM_ARRAY_TASK_ID]}" \
-SNPFile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/0.SNPData/"${SNPFiles[SLURM_ARRAY_TASK_ID]}"
