#!/bin/bash
#SBATCH --job-name=Meme_motif        # Job name
#SBATCH --partition=highmem_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request ## Should have more than 300 here
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Meme_motif.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Meme_motif.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail


ml Anaconda3/2022.10
source activate r_env

Rscript /home/sb14489/Epigenomics/scATAC-seq/0_CoreScript/Annotation_Cluster/GOTerm_forDenovoMarkers.R \
 --OutPutDir /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/6.GOTerm/AnnV4 \
 --FileNameFix _deseq_2_results.tsv \
 --WDir_forDenovo /scratch/sb14489/3.scATAC/2.Maize_ear/6.Annotation/3.Denovo/AnnV4/ \
 --GOData /scratch/sb14489/0.Reference/Maize_B73/GOTerm/B73_GO.RemoveBlank.out
