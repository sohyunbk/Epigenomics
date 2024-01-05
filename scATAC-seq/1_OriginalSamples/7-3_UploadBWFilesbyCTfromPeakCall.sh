#!/bin/bash
#SBATCH --job-name=Bigwig        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=10gb                   # Job memory request
#SBATCH --time=0:30:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9-3_Bigwig.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9-3_Bigwig.%j.err    # Standard error log
#SBATCH --array=0-2

ml Anaconda3/2023.09-0
source activate /home/sb14489/.conda/envs/ucsc

WorkingDirs=(
A619
Bif3
rel2
)

#find /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/"${WorkingDirs[SLURM_ARRAY_TASK_ID]}" \
# -name "*.macs_treat_pileup.normalized.bdg" -exec cp {} /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/BWFiles \;

cd /scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4/BWFiles

for file in "${WorkingDirs[SLURM_ARRAY_TASK_ID]}"*.normalized.bdg; do
    sorted_bdg=${file%.bdg}.sorted.bdg
    sorted_bw=${file%.bdg}.sorted.bw
    # Check if the output file does not exist
    if [ ! -f "$sorted_bw" ]; then
        sort -k1,1 -k2,2n "$file" > "$sorted_bdg"
        bedGraphToBigWig "$sorted_bdg" /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai "$sorted_bw"
    fi
done
