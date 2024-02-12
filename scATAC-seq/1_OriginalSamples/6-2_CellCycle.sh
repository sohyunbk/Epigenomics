#!/bin/bash
#SBATCH --job-name=Denovo        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Denovo.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Denovo.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-1

### Wanted to add other sample too!
ml Anaconda3/2020.02
source activate r_env
SampleName=(A619 Bif3)
MetaFiles=(A619/Ref_AnnV4_metadata.txt Bif3/Bif3_AnnV4_metadata.txt)
Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/CellCycle.R \
 --Meta /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/"${MetaFiles[SLURM_ARRAY_TASK_ID]}" \
 --Gene /scratch/sb14489/3.scATAC/0.Data/CellCycle/CellCycle.bed \
 --GAFile /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_"${SampleName[SLURM_ARRAY_TASK_ID]}"_includingZmCLE7.txt \
 --CO /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt \
 --OutputANDPreFix /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.CellCycle/"${SampleName[SLURM_ARRAY_TASK_ID]}"

#meta <-"/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/A619/Ref_AnnV4_metadata.txt"
#gene <- "/scratch/sb14489/3.scATAC/0.Data/CellCycle/CellCycle.bed"
#GA <- "/scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_includingZmCLE7.txt"
#CellOrderFile <- "/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/0.AnnotatedMeta/Ann_v4_CellType_order_forA619Bif3.txt"
#OutputPathFileName <-"/scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/2.CellCycle/A619"
