#!/bin/bash
#SBATCH --job-name=MarkerGene_bif3        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=40             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:38:00               # Time limit hrs:min:sec ## Tales ;ole 20-30min
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/MarkerGenes_Tfidf_SubmitScript.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/1.MarkerGene/bif3_Re1Re2Re3Re4 \
 --Name bif3_Re1Re2Re3Re4 \
 --meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/4Replicates/bif3/bif3_Re1234_FeaturesN163511_k50_res0.9.AfterHarmony.metadata.txt \
 --geneact /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_bif3_Re1Re2Re3Re4_withUpdatedGTF_sorted.txt \
 --pcs /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/4Replicates/bif3/bif3_Re1234.AfterHarmony.PCA.txt \
 --markers /scratch/sb14489/3.scATAC/0.Data/MarkerGene/230426_EarMarker.txt
