#!/bin/bash
#SBATCH --job-name=Alignment        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=24             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=30:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/RemoveMultiMap_Deduplication.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

OGSampleNameList=(Sohyun_A619-1 Sohyun_A619-2 Sohyun_BIF3-1 Sohyun_BIF3-2)
NewSampleNameList=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)

/home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/RemoveMultiMap_Deduplication.sh --path /scratch/sb14489/3.scATAC/2.Maize_ear \
--MappedDir 2.Mapped_CellRangerv2  --OGSampleName "${OGSampleNameList[SLURM_ARRAY_TASK_ID]}" \
 --NewSampleName_forBam "${NewSampleNameList[SLURM_ARRAY_TASK_ID]}"
