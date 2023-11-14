#!/bin/bash
#SBATCH --job-name=MarkerGene        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=2:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

ml Anaconda3/2020.02
source activate r_env
#Zm00001eb999999

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/GeneBodyA_Submit.R \
 --Ann /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_MtPtAdd_Rsf_AddZmCLE7.gtf \
 --ChrFai /scratch/sb14489/0.Reference/Maize_B73/Zm-B73-REFERENCE-NAM-5.0_MtPtAdd_Rsf.fa.fai \
 --Re1_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/2_rel2_Unique.bed \
 --Re2_bed /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/3_bif3_2_Unique.bed \
 --OutFileName /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_rel2_includingZmCLE7.txt
