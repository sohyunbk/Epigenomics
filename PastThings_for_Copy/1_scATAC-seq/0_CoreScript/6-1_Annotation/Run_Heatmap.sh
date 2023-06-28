#!/bin/bash
#SBATCH --job-name=Heatmap        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/HeatMap.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/HeatMap.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate r_env

Rscript 6-2_All_MarkerGene_Expression_dotplot_HeatMap.R
