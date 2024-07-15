#!/bin/bash
#SBATCH --job-name=ClassifyingPeaks        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ClassifyingPeaks.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ClassifyingPeaks.%j.err    # Standard error log
#SBATCH --array=0-2

source activate r_env

OutPutFilesDir=(A619_Bif3_500bpCommonPeak A619_rel2_500bpCommonPeak A619_relk1_500bpCommonPeak)
OutPutFiles=(ComA619Bif3.unique500bpPeaks ComA619rel2.unique500bpPeaks ComA619relk1.unique500bpPeaks)

## It makes bw file too!

~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/0_Functions/ClassifyPeaks_Intergenic_Genic.py \
 -PeakFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${OutPutFilesDir[SLURM_ARRAY_TASK_ID]}"/"${OutPutFiles[SLURM_ARRAY_TASK_ID]}".bed  \
 -BlackListFile /scratch/sb14489/0.Reference/Maize_B73/Zm.final_blaclist.Mito_Chloro_Chip.txt \
 -FaiFile /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai \
 -Ann_bed  /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_OnlyGene_Chr.bed \
 -OutfilePath /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${OutPutFilesDir[SLURM_ARRAY_TASK_ID]}"/ \
 -OutfileName "${OutPutFiles[SLURM_ARRAY_TASK_ID]}"
