#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=32             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/withoutCUDA2.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/withoutCUDA2.%j.err    # Standard error log

cd /scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ

/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/AlexData/Run.py \
 -wmlFile /home/sb14489/Epigenomics/ML_ACR/Deep_Former/AlexData/500bp_MappableRegions_DanQ_withoutCuda_FiveClasses.yml \
 -learningRate 0.0005
#Segmentation fault (core dumped)
