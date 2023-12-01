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
 MetaData=(A619/Ref_AnnV4_metadata.txt A619/Ref_AnnV4_metadata.txt Bif3/Bif3_AnnV3_metadata.txt Bif3/Bif3_AnnV3_metadata.txt)
 OutFileList=(1_A619_Re2 1_A619_Re1 3_bif3_Re2 3_bif3_Re1)

 ~/.conda/envs/r_env/bin/python /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/PeakCalling_byCellTypes/WithouPeakCalling_Getting_Bed_byCelltype_Macs2.py \
 -BedFile /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/4.Bam_FixingBarcode/"${BedList[SLURM_ARRAY_TASK_ID]}" \
 -Outfile /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3/5.Jbrowse_MACS2/"${OutFileList[SLURM_ARRAY_TASK_ID]}" \
 -MetaFile /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${MetaData[SLURM_ARRAY_TASK_ID]}" \
