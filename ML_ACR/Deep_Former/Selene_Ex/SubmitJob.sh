#!/bin/bash
#SBATCH --job-name=Selene_ex        # Job name
#SBATCH --partition=gpu_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=8:00:00               # Time limit hrs:min:sec
#SBATCH --gres=gpu:A100:1
#SBATCH --output=/scratch/sb14489/0.log/SeleneEx_predict.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/SeleneEx_predict.%j.err    # Standard error log


module load CUDA/11.1.1-GCC-10.2.0

cd /scratch/sb14489/8.ML_ACR/SeleneEx/

#/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/Selene_Ex/Run.py -wmlFile simple_train.yml
/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/Selene_Ex/Run.py -wmlFile 4_simple_train_MultipleClass_withSampleNegativeFalse.yml
#/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/Selene_Ex/Run.py -wmlFile 2_simple_train_DifferentLength.yml
#/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/Selene_Ex/Run.py -wmlFile 3_simple_train_ChangeSampleN.yml
