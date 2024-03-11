#!/bin/bash
#SBATCH --job-name=Fimo        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Fimo.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Fimo.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

module load MEME/5.5.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.3.0

MemeMotifDB=/scratch/sb14489/3.scATAC/0.Data/Plant_Motif_PWM/TGAATGAA_TAAT.txt
Infile_bed=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak/A619Bif3_IM-OC_MergedPeak.bed
OutfilePathName=/scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/3.fimo/IM_OC
Infile_fasta="${Infile_bed%.bed}.fasta"

#bedtools getfasta -fi /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa -bed "$Infile_bed" -fo "$Infile_fasta"
fimo --o "$OutfilePathName" "$MemeMotifDB" "$Infile_fasta"
