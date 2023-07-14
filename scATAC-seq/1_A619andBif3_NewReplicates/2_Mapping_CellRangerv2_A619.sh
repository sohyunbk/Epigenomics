#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30             # Number of CPU cores per task
#SBATCH --mem=600gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=30:00:04               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-1

Sample=(Sohyun_A619-1 Sohyun_A619-2)

module load CellRanger-ATAC/2.0.0

/home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Mapping.sh -SampleName "${Sample[SLURM_ARRAY_TASK_ID]}" \
-Ref /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scATACCellRangerv2 \
-InputPath /scratch/sb14489/3.scATAC/2.Maize_ear/1.Rawdata/ -OutfilePath /scratch/sb14489/3.scATAC/2.Maize_ear/2.Mapped_CellRangerv2/
