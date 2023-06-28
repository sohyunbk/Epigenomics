#!/bin/bash
#SBATCH --job-name=MarkerGene        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/6-1_Annotation/6-1_MarkerGenes_Tfidf_SubmitScript.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/relk1_MarkerGene221130 \
 --Name relk1 \
 --meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP_k50_res0.9.AfterHarmony.metadata.txt \
 --geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_relk1.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/rel2/rel2_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100_RemoveBLonlyMitoChloroChIP.AfterHarmony.PCA.txt \
 --markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/221130_EarMarker.txt
