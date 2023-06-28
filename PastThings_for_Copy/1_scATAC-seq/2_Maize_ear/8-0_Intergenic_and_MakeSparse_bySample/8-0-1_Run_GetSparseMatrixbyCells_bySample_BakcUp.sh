#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/8_GetSparse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/8_GetSparse.%j.err    # Standard error log
#SBATCH --array=0-1

Sample=(A619 bif3)
BedFileName=(1_A619_Combined_Sorted_k12.bed 3_bif3_Combined_Sorted_k12.bed)

source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

PeakDir=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak/

#awk '$1 ~ /^chr/' "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.bed \
# > "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.FilteredOrgs.bed

#sort -k1,1V -k2,2n "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.FilteredOrgs.bed \
# > "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.FilteredOrgs_Sorted.bed
#https://bioconda.github.io/recipes/perl-sort-naturally/README.html?highlight=perl
sort -k1,1V -k2,2n "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.bed \
 > "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks_Sorted.bed

## Common peak with A619 Bif3
bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${BedFileName[SLURM_ARRAY_TASK_ID]}" \
 -b "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.FilteredOrgs_Sorted.bed \
 -wa -wb -g  /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai_sorted \
 -sorted | /home/sb14489/.conda/envs/Jbrowse/bin/perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8-1_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl - > "$PeakDir"/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${Sample[SLURM_ARRAY_TASK_ID]}".500bp_peaks.sparse
