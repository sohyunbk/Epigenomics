#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=120:00:04               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/2_Mapping.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/2_Mapping.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-5

InputPath=(1.RawData_Re1 1.RawData_Re1 1.RawData_Re1 1.RawData_Re4 1.RawData_Re4 1.RawData_Re4)
Sample=(yinxin_Scanlon_early yinxin_Scanlon_mid yinxin_Scanlon_late Scanlon_early Scanlon_mid Scanlon_later-mid)

cd /scratch/sb14489/3.scATAC/6.LeafSheath/2.Mapped
module load CellRanger-ATAC/2.1.0

cellranger-atac count \
   --id="${Sample[SLURM_ARRAY_TASK_ID]}"  \
   --reference=/scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scATACCellRangerv2  \
   --fastqs=/scratch/sb14489/3.scATAC/6.LeafSheath/"${InputPath[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}"  \
   --localcores=24 --localmem=300 --localvmem=300
