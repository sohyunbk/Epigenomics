#!/bin/bash
#SBATCH --job-name=GACombining       # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=40             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_GACombining.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_GeneGACombining.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

##Just combinding
#cat GA_A619_Re.txt GA_A619_Re3Re4.txt > GA_A619_Re1Re2Re3Re4_withUpdatedGTF.txt
#sort -k1,1 GA_A619_Re1Re2Re3Re4_withUpdatedGTF.txt > GA_A619_Re1Re2Re3Re4_withUpdatedGTF_sorted.txt
sort -k1,1 --parallel=40 /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re1Re2Re3Re4_withUpdatedGTF.txt \
 -o /scratch/sb14489/3.scATAC/2.Maize_ear/4.Bam_FixingBarcode/GA_A619_Re1Re2Re3Re4_withUpdatedGTF_sorted.txt
