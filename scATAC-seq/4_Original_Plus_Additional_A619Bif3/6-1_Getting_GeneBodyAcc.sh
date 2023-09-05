#!/bin/bash
#SBATCH --job-name=GACombining       # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=90gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=5:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/6_MarkerGene.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/6_MarkerGene.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

##Just combinding
cat GA_A619_Re.txt GA_A619_Re3Re4.txt > GA_A619_Re1Re2Re3Re4_withUpdatedGTF.txt
sort -k1,1 GA_A619_Re1Re2Re3Re4_withUpdatedGTF.txt > GA_A619_Re1Re2Re3Re4_withUpdatedGTF_sorted.txt
