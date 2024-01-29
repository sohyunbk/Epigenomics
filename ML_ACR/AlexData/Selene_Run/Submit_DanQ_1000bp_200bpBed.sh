#!/bin/bash
#SBATCH --job-name=withoutCUDA        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=100             # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=30:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/withoutCUDA.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/withoutCUDA.%j.err    # Standard error log
#SBATCH --array=0-16

SampleName=(
  NonRedundantACRs_18Cells.200bp
  )

 /home/sb14489/miniconda3/envs/pytorch/bin/python /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Run_Selene_MakeYMLFile.py \
  --wmlFile /home/sb14489/Epigenomics/ML_ACR/Selene_AllModels_BashRun_ymlFiles/Standard_DanQ_WithoutCuda_SeqLength1000bp_200bpBed.yml \
 --learningRate 0.0005 \
 --bedfile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/"${SampleName[SLURM_ARRAY_TASK_ID]}"_Sorted.bed.gz \
 --featurefile /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/1.InputBed/"${SampleName[SLURM_ARRAY_TASK_ID]}"_distinctfeatures.txt \
 --OutwmlfileName /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${SampleName[SLURM_ARRAY_TASK_ID]}"_1000bp_200bpBed.wml \
 --NewOutputDir /scratch/sb14489/8.ML_ACR/1.MaizeGenotypes_Alex/2.Selene/"${SampleName[SLURM_ARRAY_TASK_ID]}"_DanQ_1000bp_200bpBed
