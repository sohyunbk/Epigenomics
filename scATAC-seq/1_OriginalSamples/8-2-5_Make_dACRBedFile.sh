#!/bin/bash
#SBATCH --job-name=MakeBedFileFordACR        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/MakeBedFileFordACR.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/MakeBedFileFordACR.%j.err    # Standard error log

ml Anaconda3/2023.09-0
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/dACR/MakeBedfile_fromdACRResult.R \
--DEGDir /scratch/sb14489/3.scATAC/2.Maize_ear/8.Comparative_Analysis/2.dACR/A619_vs_Bif3_AnnV4/ \
--DEGFileName  .EdgeRResult_PseudoReplicate_withPromoterRegion.txt 
