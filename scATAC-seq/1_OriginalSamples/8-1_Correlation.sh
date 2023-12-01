#!/bin/bash
#SBATCH --job-name=Correlation        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Correlation.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Correlation.%j.err    # Standard error log


ml Anaconda3/2022.10
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Correlation/Correlation_CommonPeak_Intergenic_MostVariable_FixVarErrorChange.R \
  --S1_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComPeakA619rel2_A619Barcode_Tn5Count.sparse \
  --S2_Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComPeakA619rel2_Rel2Barcode_Tn5Count.sparse \
  --S1Name A619 \
  --S2Name rel2 \
  --S1_Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt \
  --S2_Meta //scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/rel2/rel2_AnnV4.txt \
  --ClusterColumnName Ann_v4 \
  --S1and2_500bpPeak //scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed \
  --S1and2_500bpInterPeak //scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks_Intergenic.bed \
  --OutPath /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/1.Correlation \
  --OutFileName A619andrel2
