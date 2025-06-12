#!/bin/bash
#SBATCH --job-name=Peak_Barcode_Tn5        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/7_Peak_Barcode_Tn5.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/7_Peak_Barcode_Tn5.%j.err    # Standard error log
#SBATCH --array=0-10

## If the file is small it takes 2 hour
## If it's the big file... it takes 6-7 hours
Tn5BedNames=(1_A619 3_bif3 3_bif3 1_A619 2_rel2 1_A619 4_relk1 1_A619)
PeakFiles=(
A619/A619.500bp_peaks_Intergenic.bed
Bif3/Bif3.500bp_peaks_Intergenic.bed
A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed
A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed
A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed
A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed
A619_relk1_500bpCommonPeak/ComA619relk1.unique500bpPeaks.bed
A619_relk1_500bpCommonPeak/ComA619relk1.unique500bpPeaks.bed)
OutFileNames=(
A619/A619Peak500bp_A619Barcode_Tn5Count
Bif3/Bif3Peak500bp_Bif3Barcode_Tn5Count
A619_Bif3_500bpCommonPeak/ComA619Bif3_Bif3Barcode_Tn5Count
A619_Bif3_500bpCommonPeak/ComA619Bif3_A619Barcode_Tn5Count
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_Rel2Barcode_Tn5Count
/A619_rel2_500bpCommonPeak/ComPeakA619rel2_A619Barcode_Tn5Count
/A619_relk1_500bpCommonPeak/ComPeakA619relk1_Relk1Barcode_Tn5Count
/A619_relk1_500bpCommonPeak/ComPeakA619relk1_A619Barcode_Tn5Count)

ml Anaconda3/2022.10
source activate r_env

Rscript ../workflow_scripts/PeakCalling_byCellTypes/Peak_Barcode_Tn5NotBianry_Sparse.R \
 --Peak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${PeakFiles[SLURM_ARRAY_TASK_ID]}" \
 --Re1_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Tn5BedNames[SLURM_ARRAY_TASK_ID]}"_Unique.bed \
 --Re2_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Tn5BedNames[SLURM_ARRAY_TASK_ID]}"_2_Unique.bed \
 --OutFileName /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${OutFileNames[SLURM_ARRAY_TASK_ID]}".sparse
