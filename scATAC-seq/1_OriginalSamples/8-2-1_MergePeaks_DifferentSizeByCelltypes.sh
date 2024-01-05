#!/bin/bash
#SBATCH --job-name=MergePeak        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MergePeak.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MergePeak.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

/home/sb14489/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8-2_Peaks_Cut/MergeCutPeaks_toFixHumpedPeaks.py -method Method1 \
    -inputpath1 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619 \
     -inputpath2 /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/Bif3 \
      -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai \
      -OutputFileName A619Bif3 -OutputPath /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_Bif3_MergedDifferentSizePeak \
      -Ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
      -Sparse1 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined_Sorted_k12.bed_OnlyChr \
      -Sparse2 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined_Sorted_k12.bed_OnlyChr
