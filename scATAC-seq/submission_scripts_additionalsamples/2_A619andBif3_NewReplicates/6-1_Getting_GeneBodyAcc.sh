#!/bin/bash
#SBATCH --job-name=MarkerGene_A619        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail
#SBATCH --array=0-1                  # Array range

ml Anaconda3/2020.02
source activate r_env
#Zm00001eb999999
Re1=(A619_Re3_Unique.bed bif3_Re3_Unique.bed)
Re2=(A619_Re4_Unique.bed bif3_Re4_Unique.bed)
Sample=(A619_Re3Re4 bif3_Re3Re4)

Rscript ../workflow_scripts/GeneBodyA_Submit.R \
 --Ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf_AddZmCLE7.gtf \
 --ChrFai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
 --Re1_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Re1[SLURM_ARRAY_TASK_ID]}" \
 --Re2_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/"${Re2[SLURM_ARRAY_TASK_ID]}" \
 --OutFileName /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_"${Sample[SLURM_ARRAY_TASK_ID]}".txt
