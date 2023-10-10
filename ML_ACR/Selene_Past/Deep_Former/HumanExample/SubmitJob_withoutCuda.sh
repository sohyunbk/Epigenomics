#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=32             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/withoutCUDA.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/withoutCUDA.%j.err    # Standard error log
#SBATCH --array=0

#ymllist=(500bp_MappableRegions_DanQ_withoutCuda.yml)
#learningRatelist=(0.0005)

ymllist=(/home/sb14489/ACR_ML_caQTLs/Deep_Former/HumanExample/DanQ_JASPAR_YAML.yml)
learningRatelist=(0.005)

/home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/ACR_ML_caQTLs/Deep_Former/AlexData/Run.py \
 -wmlFile "${ymllist[SLURM_ARRAY_TASK_ID]}" \
 -learningRate "${learningRatelist[SLURM_ARRAY_TASK_ID]}"

#srun --pty  -p gpu_p --gres=gpu:A100:1   --mem=1G --nodes=1 --ntasks-per-node=1 --time=1:00:00 --job-name=qlogin /bin/bash -l
