#!/bin/bash
#SBATCH --job-name=ChromVar        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/ChromVar.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/ChromVar.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate JASPAR_act

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MotifDeviation/MotifAnalysis_ChromVar_AllMotifs_SmoothData.R  \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AnnV4 \
 --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619/A619Peak500bp_A619Barcode_Tn5Count.sparse \
 --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt \
 --Markov /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_IncludingZmCLE7/A619_IncludingZmCLE7.MarkovMatrix.rds \
 --SampleName A619 --IGPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619/A619.500bp_peaks_Intergenic.bed

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MotifDeviation/MotifAnalysis_ChromVar_AllMotifs_SmoothData.R \
  --WD /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AnnV4 \
  --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/Bif3/Bif3Peak500bp_Bif3Barcode_Tn5Count.sparse \
  --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Bif3/Bif3_AnnV4_metadata.txt \
  --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
  --Markov /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_IncludingZmCLE7/Bif3_IncludingZmCLE7.MarkovMatrix.rds \
  --SampleName Bif3 --IGPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/Bif3/Bif3.500bp_peaks_Intergenic.bed
