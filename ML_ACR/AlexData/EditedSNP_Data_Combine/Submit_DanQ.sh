#!/bin/bash
#SBATCH --job-name=schmitz_hm_p        # Job name
#SBATCH --partition=gpu_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --gres=gpu:P100:1
#SBATCH --output=/scratch/sb14489/0.log/DanQ_predict.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/DanQ_predict.%j.err    # Standard error log

module load CUDA/11.1.1-GCC-10.2.0

/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/Selene_Ex_RunFast/Run.py \
 -wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_Ex_RunFast/Standard_DanQ_WithoutCuda_SeqLength500bp_bedFile.yml \
 -learningRate 0.0005

# -wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_Ex_RunFast/Evaludate_test_bed.yml
