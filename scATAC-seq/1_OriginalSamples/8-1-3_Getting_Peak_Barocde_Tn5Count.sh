#!/bin/bash
#SBATCH --job-name=Peak_Barcode_Tn5        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/7_Peak_Barcode_Tn5.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/7_Peak_Barcode_Tn5.%j.err    # Standard error log
#SBATCH --array=0-1

## If the file is small it takes 2 hour
## If it's the big file... it takes 6-7 hours
Tn5BedNames=(3_bif3 1_A619 2_rel2 1_A619)
PeakFiles=(A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed
A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed
A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed
A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed)
OutFileNames=(A619_Bif3_500bpCommonPeak/ComA619Bif3_Bif3Barcode_Tn5Count
A619_Bif3_500bpCommonPeak/ComA619Bif3_A619Barcode_Tn5Count
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_Rel2Barcode_Tn5Count
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_A619Barcode_Tn5Count)

ml Anaconda3/2022.10
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/PeakCalling_byCellTypes/Peak_Barcode_Tn5NotBianry_Sparse.R \
 --Peak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${PeakFiles[SLURM_ARRAY_TASK_ID]}" \
 --Re1_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Tn5BedNames[SLURM_ARRAY_TASK_ID]}"_Unique.bed \
 --Re2_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Tn5BedNames[SLURM_ARRAY_TASK_ID]}"_2_Unique.bed \
 --OutFileName /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${OutFileNames[SLURM_ARRAY_TASK_ID]}".sparse