#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/7_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/7_PeakCalling.%j.err    # Standard error log

ml Anaconda3/2022.10
source activate r_env

cd $Path
~/.conda/envs/r_env/bin/python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/PeakCalling_byCellTypes/MergingACRs_500bpFixed_withHighestTn5.py \
 -Path1_forCellPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619 \
 -Path2_forCellPeak /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/rel2 \
 -OutPutFile /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/A619_rel2_500bpCommonPeak/ComA619rel2.unique500bpPeaks.bed
