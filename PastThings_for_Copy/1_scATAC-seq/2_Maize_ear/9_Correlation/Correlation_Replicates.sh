#!/bin/bash
#SBATCH --job-name=Correlation        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Correlation.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Correlation.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/9_Correlation/9_Correlation_ByReplicates.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/9.Correlation \
 --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks.sparse \
 --SampleName bif3 \
 --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt \
 --PeakI /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks_Intergenic.bed \
 --PeakG /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/bif3/bif3.500bp_peaks_Genic.bed

#Sparsefile_A619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks.sparse"
#MetaFileA619 <- "/scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt"
#SampleName <- "A619"
#Peak_Inter <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Intergenic.bed",header=F)
#Peak_Genic <-read.table("/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/A619/A619.500bp_peaks_Genic.bed",header=F)
