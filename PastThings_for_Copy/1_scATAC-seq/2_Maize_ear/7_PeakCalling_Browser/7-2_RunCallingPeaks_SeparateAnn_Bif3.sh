#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9_PeakCalling.%j.err    # Standard error log
#SBATCH --array=0

Path=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_SeparateAnn_Mutant
Combined=(3_bif3)
Sample=(bif3)

ml Anaconda3/2020.02
source activate r_env

cd $Path

~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/7_PeakCalling_Browser/7-2_call_scACRs.py \
 -bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted.bed \
 -meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt \
 -col Ann_v3 -base "${Sample[SLURM_ARRAY_TASK_ID]}" -outdir "${Sample[SLURM_ARRAY_TASK_ID]}" \
 -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai -bw TRUE
