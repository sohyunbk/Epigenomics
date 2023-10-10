#!/bin/bash
#SBATCH --job-name=Pytorch_Seedling        # Job name
#SBATCH --partition=gpu_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=50:00:00               # Time limit hrs:min:sec
#SBATCH --gres=gpu:A100:1
#SBATCH --output=/scratch/sb14489/0.log/Pytorch_Seedling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Pytorch_Seedling.%j.err    # Standard error log


ml Anaconda3/2020.02
## 2) conda activate
source activate pytorch

## 3) module load
module load CUDA/11.1.1-GCC-10.2.0

cd /scratch/sb14489/8.ML_ACR/2.DeepFormer_DanQ

/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/DanQ_train.py
