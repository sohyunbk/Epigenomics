#!/bin/bash
#SBATCH --job-name=Python        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=50:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Python.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Python.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate r_env

/home/sb14489/.conda/envs/r_env/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/AlexData/ReduceNegative_fromMappableRegions.py
