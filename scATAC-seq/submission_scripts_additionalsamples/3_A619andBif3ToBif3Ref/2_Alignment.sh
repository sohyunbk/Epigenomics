#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=167:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours -for mapping CellRanger
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3                   # Array range

#sbatch /home/sb14489/1.scATAC-seq/1_scATAC-seq/workflow_scripts/2_Alignment/2-1_Alignment_CellRangerv2.sh
cd /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/1.Mapped
Sample=(1_A619 3_bif3 1_A619_2 3_bif3_2)
module load CellRanger-ATAC/2.1.0

cellranger-atac count \
   --id="${Sample[SLURM_ARRAY_TASK_ID]}"  \
   --reference=/scratch/sb14489/0.Reference/Maize_Ki3/Zm-Ki3_OnlyChr_scATACCellRangerv2_Bif3  \
   --fastqs=/scratch/sb14489/3.scATAC/2.Maize_ear/1.Rawdata/"${Sample[SLURM_ARRAY_TASK_ID]}"  \
   --localcores=24 --localmem=200 --localvmem=100
