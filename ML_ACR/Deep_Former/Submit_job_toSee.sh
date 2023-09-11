#!/bin/bash
#SBATCH --job-name=Pytorch        # Job name
#SBATCH --partition=gpu_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=43:00:00               # Time limit hrs:min:sec
#SBATCH --gres=gpu:P100:1
#SBATCH --output=/scratch/sb14489/0.log/Pytorch.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Pytorch.%j.err    # Standard error log


ml Anaconda3/2020.02
## 2) conda activate
source activate pytorch

## 3) module load
module load CUDA/11.1.1-GCC-10.2.0

cd /scratch/sb14489/8.ML_ACR/DeepFormer_Ex/DeepFormer/maize_code

/home/sb14489/miniconda3/envs/pytorch/bin/python DanQ_train.py
