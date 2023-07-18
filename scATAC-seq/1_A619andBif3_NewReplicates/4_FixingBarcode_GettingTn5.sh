#!/bin/bash
#SBATCH --job-name=FixingBarcodes        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=400gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : 80 hours
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcodes.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcodes.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --array=0-3

SampleNameList=(A619_Re3 A619_Re4 bif3_Re3 bif3_Re4)
