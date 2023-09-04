#!/bin/bash
#SBATCH --job-name=MarkerGene_bif3        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/bif3_Re3Re4_TSS35_FRiP55_MinT0.007_MaxT0.005 \
 --Name bif3_Re3Re4_TSS35_FRiP55_MinT0.007_MaxT0.005 \
 --meta //scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/bif3/bif3_Tn5Cut1000_Binsize500_MinT0.007_MaxT0.005_PC100_FeaturesN2e+05_k50_res0.9.AfterHarmony.metadata.txt \
 --geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_bif3_Re3Re4.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS35_FRiP55/bif3/bif3_Tn5Cut1000_Binsize500_MinT0.007_MaxT0.005_PC100.AfterHarmony.PCA.txt \
 --markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt
