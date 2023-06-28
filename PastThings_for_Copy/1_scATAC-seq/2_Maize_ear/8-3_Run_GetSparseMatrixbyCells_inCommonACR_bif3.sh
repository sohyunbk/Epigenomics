#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/8_GetSparse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/8_GetSparse.%j.err    # Standard error log
source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

cd /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/bif3_AnnV3

sort -k1,1V -k2,2n bif3.unique500bpPeaks.bed > bif3.unique500bpPeaks_sorted.bed
#https://bioconda.github.io/recipes/perl-sort-naturally/README.html?highlight=perl

## Common peak with A619 Bif3
bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined_Sorted_k12.bed \
 -b bif3.unique500bpPeaks_sorted.bed  -wa -wb -g  /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai_sorted \
 -sorted | /home/sb14489/.conda/envs/Jbrowse/bin/perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl - > ../A619_bif3_AnnV3/bif3_CommonpeakwithA619_sorted.sparse
