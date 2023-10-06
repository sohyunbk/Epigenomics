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

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/10_MotifAnalysis/10-1_MotifAnalysis_ChromVar_AllMotifs_SmoothData.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif \
 --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/A619_toComPeak.sparse \
 --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_RemoveBLonlyMitoChloroChIP.REF_CELLs.reduced_dimensions.txt \
 --Markov /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/A619_Markov/A619.MarkovMatrix.rds \
 --SampleName A619 --IGPeak /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_Intergenic.bed

 Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MotifAnalysis_ChromVar_AllMotifs_SmoothData.R \
  --WD /scratch/sb14489/3.scATAC/2.Maize_ear/10.MotifAnalysis/1.ChromVar/AddOCdACRMotif \
  --Sparse /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/Bif3_toComPeak.sparse \
  --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt \
  --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/bif3_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
  --Markov /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/Bif3_Markov/Bif3.MarkovMatrix.rds \
  --SampleName Bif3 --IGPeak /scratch/sb14489/3.scATAC/2.Maize_ear/8.CommonACRs/A619_Bif3/ComA619Bif3.unique500bpPeaks_Intergenic.bed
