source activate /home/sb14489/.conda/envs/Jbrowse

module load BEDTools/2.30.0-GCC-10.2.0
cd /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619

sort -k1,1V -k2,2n A619.500bp_peaks.bed > A619.500bp_peaks_sorted.bed
#https://bioconda.github.io/recipes/perl-sort-naturally/README.html?highlight=perl

bedtools intersect -a /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/1_A619_Combined_Sorted.bed \
  -b A619.500bp_peaks_sorted.bed  -wa -wb -g  /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai_sorted \
  -sorted | perl ~/1.scATAC-seq/1_scATAC-seq/0_CoreScript/8_FindCommonACRPos/8-3_fastSparse.nonbinary.peak.pl - > A619.500bp_peaks.sparse
