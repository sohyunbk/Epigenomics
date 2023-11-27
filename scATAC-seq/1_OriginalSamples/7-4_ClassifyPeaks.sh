#!/bin/bash
#SBATCH --job-name=ClassifyingPeaks        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ClassifyingPeaks.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ClassifyingPeaks.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate r_env

## It makes bw file too!

~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/0_Functions/ClassifyPeaks_Intergenic_Genic.py \
 -PeakFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed  \
 -BlackListFile /scratch/sb14489/0.Reference/Maize_B73/Zm.final_blaclist.Mito_Chloro_Chip.txt \
 -FaiFile /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
 -Ann_bed  /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
 -OutfilePath /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ \
 -OutfileName ComA619rel2.unique500bpPeaks
