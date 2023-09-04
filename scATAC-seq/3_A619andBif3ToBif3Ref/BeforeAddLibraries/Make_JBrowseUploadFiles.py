#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=70gb                   # Job memory request
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Bigwig.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate /home/sb14489/.conda/envs/ucsc

mkdir BwFiles

SampleName=3_bif3_Re1
bedSort  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM.bdg  ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg
bedGraphToBigWig ./"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}"_treat_pileup_CPM_Sorted.bdg /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai ./BwFiles/"$SampleName"_"${ClusterN[SLURM_ARRAY_TASK_ID]}".bw
