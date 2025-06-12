#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/7_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/7_PeakCalling.%j.err    # Standard error log
#SBATCH --array=0-2

Path=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V4
Combined=(1_A619 3_bif3 2_rel2)
Sample=(A619 Bif3 rel2)
MetaData=(Ref_AnnV4_metadata.txt Bif3_AnnV4_metadata.txt rel2_AnnV4.txt)
ml Anaconda3/2022.10
source activate r_env

cd $Path
~/.conda/envs/r_env/bin/python ../workflow_scripts/PeakCalling_byCellTypes/Call_scACRs_byCellTypes_WithoutFakePeak.py \
 -bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted_k12.bed_OnlyChr \
 -meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${Sample[SLURM_ARRAY_TASK_ID]}"/"${MetaData[SLURM_ARRAY_TASK_ID]}" \
 -col Ann_v4 -base "${Sample[SLURM_ARRAY_TASK_ID]}" -outdir "${Sample[SLURM_ARRAY_TASK_ID]}" \
 -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_OnlyChr.fa.fai -bw TRUE
