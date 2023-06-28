#!/bin/bash
#SBATCH --job-name=11_memeChip        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=60gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/11_memeChip.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/11_memeChip.%j.err    # Standard error log
#SBATCH --array=0-1

module load MEME/5.4.1-foss-2019b-Python-3.7.4
Sample=(IM-OC.Bif3Higher IM-OC.A619Higher)

meme-chip //scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_PromoterIntergenic/"${Sample[SLURM_ARRAY_TASK_ID]}".fa \
 -minw 4 -maxw 15 -o //scratch/sb14489/3.scATAC/2.Maize_ear/11.dACRs/A619_vs_Bif3_PromoterIntergenic/"${Sample[SLURM_ARRAY_TASK_ID]}"
