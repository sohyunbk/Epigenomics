#!/bin/bash
#SBATCH --job-name=FindMotifs        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/FindMotifs.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/FindMotifs.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

module load MEME/5.4.1-foss-2019b-Python-3.7.4
module load BEDTools/2.30.0-GCC-10.2.0

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic.bed \
-b /scratch/sb14489/7.DAPorChIP/DAPseq_WUS/HB67_WUS1_B73v5_Q30_qval5_finalBl/HB67_WUS1_B73v5_Q30_qval5_finalBl.GEM_events.narrowPeak \
-wa > /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic_OverlappedWithWUS1DAP.bed


sh FindMotif_FromPeak.sh -BedFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic_OverlappedWithWUS1DAP.bed \
 -Ref /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa  \
 -Outfile /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/
