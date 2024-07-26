#!/bin/bash
#SBATCH --job-name=Alignment_lowTime        # Job name
#SBATCH --partition=schmitz_hm_p         # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=100             # Number of CPU cores per task
#SBATCH --mem=500gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/starsolo.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/starsolo.%j.err    # Standard error log

ml Anaconda3/2020.02
source activate r_env
#mkdir /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNAstarsolo
STAR  --runMode genomeGenerate --runThreadN 10 --genomeDir /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNAstarsolo \
 --genomeFastaFiles /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNACellRanger/fasta/genome.fa \
  --sjdbGTFfile /scratch/sb14489/0.Reference/Maize_B73/Maize_B73_V5_withMtPt_scRNACellRanger/genes/genes.gtf
