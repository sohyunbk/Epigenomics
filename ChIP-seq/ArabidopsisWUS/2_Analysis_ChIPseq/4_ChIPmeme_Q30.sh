#!/bin/bash
#SBATCH --job-name=memeChip        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request #For normal fastq : 600gb
#SBATCH --time=10:00:00               # Time limit hrs:min:sec #For normal fastq : more than at least five days with 14 cpu 80 hours
#SBATCH --output=/scratch/sb14489/0.log/MemeChip.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MemeChip.%j.err    # Standard error log
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-12.2.0


bedtools getfasta -fi /scratch/sb14489/0.Reference/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
 -bed /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling_HighCutoff/WUS_GS_Q0.001_peaks.narrowPeak_Q30  \
 > /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling_HighCutoff/WUS_GS_Q0.001_peaks.narrowPeak_Chr_Q30.fa

meme-chip /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/3.PeakCalling_HighCutoff/WUS_GS_Q0.001_peaks.narrowPeak_Chr_Q30.fa \
 -minw 4 -maxw 15 -fimo-skip  -o /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/4.MemeChip/WUS_GS_Q0.001_peaks.narrowPeak_Chr_Q30
