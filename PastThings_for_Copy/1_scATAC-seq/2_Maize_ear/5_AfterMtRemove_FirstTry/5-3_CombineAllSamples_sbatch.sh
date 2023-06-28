#!/bin/bash
#SBATCH --job-name=QC        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-3_QC.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-3_QC.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/5_CellClustering/5-3_CombineAllSamples_CleanData.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/AfterMtMapping/ \
 --NewDir CombineAll --PreFix Tn5Cut1000_Binsize500_Mt0.1 --MinPeakNumber 100 --MinT 0.001 --MaxT 0.005
