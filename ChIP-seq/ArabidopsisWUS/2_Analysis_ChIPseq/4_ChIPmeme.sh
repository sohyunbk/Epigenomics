#!/bin/bash
#SBATCH --job-name=memeChip        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/MemeChip.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MemeChip.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

module load MEME/5.4.1-foss-2019b-Python-3.7.4
module load BEDTools/2.30.0-GCC-10.2.0

cd /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling

bedtools getfasta -fi /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
 -bed WUS_GS_NonModel_peaks.narrowPeak  > WUS_GS_NonModel_peaks.narrowPeak.fa

meme-chip WUS_GS_NonModel_peaks.narrowPeak.fa \
 -minw 4 -maxw 15 -fimo-skip  -o /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/4.MemeChip/Arabiopsis_WUS_GS_NonModel_peaks.narrowPeak
