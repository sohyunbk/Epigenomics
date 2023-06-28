#!/bin/bash
#SBATCH --job-name=3_QC        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/5-3_QC.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/5-3_QC.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-7                   # Array range

List=(1_A619 1_A619_2 2_rel2 2_rel2_2 3_bif3 3_bif3_2  4_relk1 4_relk1_2)

ml Anaconda3/2020.02
source activate r_env

Rscript /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/5_CellClustering_CombineLater/5-3_CleanFeature_DetectDoublets.R \
 --WD /scratch/sb14489/3.scATAC/2.Maize_ear/6.CellClustering/Organelle5Per_CombineLater/ \
 --Name "${List[SLURM_ARRAY_TASK_ID]}" --PreFix Tn5Cut1000_Binsize500_Mt0.05 --MinT 0.001 --MaxT 0.005 --nPC 100
