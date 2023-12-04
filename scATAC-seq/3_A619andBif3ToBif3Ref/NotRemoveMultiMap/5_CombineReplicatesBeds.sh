#!/bin/bash
#SBATCH --job-name=FixingBarcode        # Job name
#SBATCH --partition=batch             # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=12             # Number of CPU cores per task
#SBATCH --mem=50gb                   # Job memory request
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/4_FixingBarcodes.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/4_FixingBarcodes.%j.err    # Standard error log
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Sohyun.Bang@uga.edu  # Where to send mail

cd /scratch/sb14489/3.scATAC/4.Bif3Ref_Ki3_NotRemoveMultiMap//4.Bam_FixingBarcode/
cat 1_A619_Unique.bed 1_A619_2_Unique.bed > A619_Re1andRe2.bed
cat 3_bif3_Unique.bed 3_bif3_2_Unique.bed > Bif3_Re1andRe2.bed
