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

#meme-chip /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/HB67_WUS1_B73v5_Q30_qval5_finalBl/HB67_WUS1_B73v5_Q30_qval5_finalBl.GEM_events.narrowPeak.fa \
# -minw 4 -maxw 15 -fimo-skip  -o /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.meme-chip/Maize_ZmWUS1_DAPseq

meme-chip /scratch/sb14489/7.DAPorChIP/CHIPseq_Ara_WUS/1.SRADownload/GSM3474971_WUS_Treatment_vs_Control.fa \
 -minw 4 -maxw 15 -fimo-skip  -o /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/2.meme-chip/Arabiopsis_WUS_Treatment_vs_Control
