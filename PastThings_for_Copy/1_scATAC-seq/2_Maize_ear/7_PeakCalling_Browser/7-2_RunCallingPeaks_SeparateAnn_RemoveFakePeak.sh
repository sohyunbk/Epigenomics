#!/bin/bash
#SBATCH --job-name=PeakCalling        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=200gb                   # Job memory request
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/9_PeakCalling.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/9_PeakCalling.%j.err    # Standard error log
#SBATCH --array=0-1

Path=/scratch/sb14489/3.scATAC/2.Maize_ear/7.PeakCalling/Ann_V3_RemoveFakePeak
Combined=(1_A619 3_bif3)
Sample=(A619 bif3)
MetaData=(Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt)
ml Anaconda3/2020.02
source activate r_env
#source activate /home/sb14489/.conda/envs/ucsc
#module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

cd $Path

~/.conda/envs/r_env/bin/python /home/sb14489/1.scATAC-seq/1_scATAC-seq/0_CoreScript/7_PeakCalling_Browser/7-2_call_scACRs_WithoutFakePeak.py \
 -bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Combined[SLURM_ARRAY_TASK_ID]}"_Combined_Sorted.bed \
 -meta /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/"${MetaData[SLURM_ARRAY_TASK_ID]}" \
 -col Ann_v3 -base "${Sample[SLURM_ARRAY_TASK_ID]}" -outdir "${Sample[SLURM_ARRAY_TASK_ID]}" \
 -fai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai -bw TRUE
