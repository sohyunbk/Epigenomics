#!/bin/bash
#SBATCH --job-name=12_GettingBed        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=300gb                   # Job memory request
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/12_GettingBed.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/12_GettingBed.%j.err    # Standard error log
#SBATCH --array=0-3                  # Array range

BedList=(1_A619_2_Unique.bed 1_A619_Unique.bed 3_bif3_2_Unique.bed 3_bif3_Unique.bed)
MetaData=(Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt Ref_AfterMt0.5Cutoff/Tn5Cut1000_Binsize500_Mt0.05_MinT0.01_MaxT0.05_PC100/Ref_AnnV3_metadata.txt  Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt Organelle5Per_CombineLater/bif3/Bif3_AnnV3_metadata.txt)
OutFileList=(1_A619_Re2 1_A619_Re1 3_bif3_Re2 3_bif3_Re1)

python /home/sb14489/SingleCell_Pipelines/scATAC-seq/0_CoreScript/UploadJbrowser/Tn5InsertionByCelltype_RunningMacs2.py \
-BedFile /scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode/"${BedList[SLURM_ARRAY_TASK_ID]}" \
-Outfile /scratch/sb14489/3.scATAC/4.Bif3Ref/5.Jbrowse_MACS2/"${OutFileList[SLURM_ARRAY_TASK_ID]}" \
-MetaFile /scratch/sb14489/3.scATAC/2.Maize_ear/5.CellClustering/"${MetaData[SLURM_ARRAY_TASK_ID]}"
