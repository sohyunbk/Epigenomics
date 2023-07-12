#!/bin/bash
#SBATCH --job-name=DanQ_predict        # Job name
#SBATCH --partition=gpu_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=3:00:00               # Time limit hrs:min:sec
#SBATCH --gres=gpu:P100:1
#SBATCH --output=/scratch/sb14489/0.log/DanQ_predict.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/DanQ_predict.%j.err    # Standard error log
#SBATCH --array=0-1

module load CUDA/11.1.1-GCC-10.2.0

cd /scratch/sb14489/8.ML_ACR/2.DeeperDeepSea_Test

#ymllist=(200bp_AllGenome.yml 204bp_AllGenome.yml 204bp_MappableRegions.yml 500bp_AllGenome.yml)
ymllist=(500bp_AllGenome_DanQ.yml 500bp_MappableRegions_DanQ.yml)
learningRatelist=(0.0005 0.0005)

/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/AlexData/Run.py \
 -wmlFile "${ymllist[SLURM_ARRAY_TASK_ID]}" \
 -learningRate "${learningRatelist[SLURM_ARRAY_TASK_ID]}"
