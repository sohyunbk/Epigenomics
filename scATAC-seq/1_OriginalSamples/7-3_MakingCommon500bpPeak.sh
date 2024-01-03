#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/7_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/7_PeakCalling.%j.err    # Standard error log
#SBATCH --array=0-1

ml Anaconda3/2022.10
source activate r_env

Path1=(A619 A619)
Path2=(Bif3 rel2)
OutPutFilesNames=(A619_Bif3_500bpCommonPeak/ComA619Bif3.unique500bpPeaks.bed A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed)

~/.conda/envs/r_env/bin/python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/PeakCalling_byCellTypes/MergingACRs_500bpFixed_withHighestTn5.py \
 -Path1_forCellPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${Path1[SLURM_ARRAY_TASK_ID]}" \
 -Path2_forCellPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${Path2[SLURM_ARRAY_TASK_ID]}" \
 -OutPutFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${OutPutFilesNames[SLURM_ARRAY_TASK_ID]}"
