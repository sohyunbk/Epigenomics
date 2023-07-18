#!/bin/bash
#SBATCH --job-name=FixingBarcode_GetTn5        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcode_GetTn5.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcode_GetTn5.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

SampleNameList=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)


module load Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

module load  SAMtools/1.10-iccifort-2019.5.281
samtools index -@ 20 /scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam

python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/FindTn5Insertion_fromBam.py  \
-BAM //scratch/sb14489/3.scATAC/2.Maize_ear/3.SortedBam/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"_Rmpcr.bam \
-exp_name "${SampleNameList[SLURM_ARRAY_TASK_ID]}" -threads 20 \
-output_file /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${SampleNameList[SLURM_ARRAY_TASK_ID]}"
