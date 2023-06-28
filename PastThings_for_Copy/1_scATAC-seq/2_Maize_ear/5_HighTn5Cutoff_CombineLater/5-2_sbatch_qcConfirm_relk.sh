#!/bin/bash
#SBATCH --job-name=QC        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-2_QC.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-2_QC.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0                   # Array range

List=(1_A619)

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/5_CellClustering_CombineLater/5-2_QCConfirm.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/ \
 --BinSize 500 --Name "${List[SLURM_ARRAY_TASK_ID]}" \
 --MinTn5 1000 --MtRatio 0.2
