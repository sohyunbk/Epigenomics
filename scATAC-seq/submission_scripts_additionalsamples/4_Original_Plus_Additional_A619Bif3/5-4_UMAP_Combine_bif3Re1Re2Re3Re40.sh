#!/bin/bash
#SBATCH --job-name=4_UMAP_bif3        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-4_UMAP_Comb4Re_bif3.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-4_UMAP_Comb4Re_bif3.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript ../workflow_scripts/UMAP_4Replicates.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/4Replicates/ \
 --OldRDS /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/CombineAll/Combined_Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100.rds \
 --Re1 bif3_Re1  --Re2 bif3_Re2 --Re3 bif3_Re3 --Re4 bif3_Re4 \
 --SampleS bif3 --PreOptions_forRe3Re4 Tn5Cut1000_Binsize500_MinT0.005_MaxT0.01_PC100 \
 --WD_forRe3Re4 /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/AdditionalSample_TSS3_FRiP4/
