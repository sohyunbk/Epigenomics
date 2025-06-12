#!/bin/bash
#SBATCH --job-name=12_GettingBed        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=10             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/12_GettingBed.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/12_GettingBed.%j.err    # Standard error log
#SBATCH --array=0-5                  # Array range

 BedList=(1_A619_2_Unique.bed 1_A619_Unique.bed 3_bif3_2_Unique.bed 3_bif3_Unique.bed A619_Re1andRe2.bed Bif3_Re1andRe2.bed)
 MetaData=(A619/Ref_AnnV4_metadata.txt A619/Ref_AnnV4_metadata.txt Bif3/Bif3_AnnV3_metadata.txt Bif3/Bif3_AnnV3_metadata.txt A619/Ref_AnnV4_metadata.txt Bif3/Bif3_AnnV3_metadata.txt)
 OutFileList=(1_A619_Re2 1_A619_Re1 3_bif3_Re2 3_bif3_Re1 A619_Re1andRe2 Bif3_Re1andRe2)

 ml Anaconda3/2022.10
 source activate r_env

## I do not understand why the last sample has problem when I run this code.
## If I run it twice, things look fine,,.

 ~/.conda/envs/r_env/bin/python ../workflow_scripts/PeakCalling_byCellTypes/CellTypeBW_fromBed_KnownCT.py \
 -BedFile /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/4.Bam_FixingBarcode/"${BedList[SLURM_ARRAY_TASK_ID]}" \
 -Outfile /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/5.Jbrowse_MACS2/"${OutFileList[SLURM_ARRAY_TASK_ID]}" \
 -MetaFile /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${MetaData[SLURM_ARRAY_TASK_ID]}" \
 --Thread 10 --Fai /scratch/sb14489/0.Reference/Maize_Ki3/Zm-Ki3-REFERENCE-NAM-1.0_OnlyChr_Bif3.fa.fai
