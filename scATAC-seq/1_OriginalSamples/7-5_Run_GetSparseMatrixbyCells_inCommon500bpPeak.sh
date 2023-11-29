#!/bin/bash
#SBATCH --job-name=GettingSparseMatrix        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/8_GetSparse.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/8_GetSparse.%j.err    # Standard error log

source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0

cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak
#~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-2_adjustPeaks.sh ComA619Bif3
#sort -k1,1V -k2,2n ComA619Bif3.unique500bpPeaks.bed > ComA619Bif3.unique500bpPeaks_sorted.bed
#sort -k1,1 -k2,2n ComA619Bif3.unique500bpPeaks.bed > ComA619Bif3.unique500bpPeaks_sorted.bed

#https://bioconda.github.io/recipes/perl-sort-naturally/README.html?highlight=perl

##Upload bedfile to jbrowse
#~/jbrowse/bin/flatfile-to-json.pl --bed ComA619Bif3.unique500bpPeaks_sorted.bed --trackLabel ComA619Bif3.unique500bpPeaks_sorted.bed  --out ./

### Check Where the erro come from!!! it's from -g
#bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_Combined_AnnV3.bed  \
#  -b ComA619Bif3.unique500bpPeaks_sorted.bed  -wa -wb  \
#  -sorted  > BedToolIntersect_A619Tn5_toComPeak.txt

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined_AnnV3.bed  \
  -b ComA619rel2.unique500bpPeaks.bed  -wa -wb \
  -sorted  | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl - > A619_toComPeak.sparse

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_Combined_AnnV3.bed \
    -b ComA619rel2.unique500bpPeaks.bed  -wa -wb \
    -sorted | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl - > Bif3_toComPeak.sparse
