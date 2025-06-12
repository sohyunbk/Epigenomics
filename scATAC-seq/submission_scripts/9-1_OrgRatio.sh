#!/bin/bash
#SBATCH --job-name=OrgRatio        # Job name
#SBATCH --partition=batch         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/OrgRatio.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/OrgRatio.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-5

Samples=(1_A619, 1_A619_2, 2_rel2, 2_rel2_2, 3_bif3, 3_bif3_2)

ml Anaconda3/2022.10
source activate r_env

Rscript  /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Viualization/Organelle_QC_Plot.R \
--SampleName "${Samples[SLURM_ARRAY_TASK_ID]}" \
--Path /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/
