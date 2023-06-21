#!/bin/bash
#SBATCH --job-name=FindSamLines        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=14:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/FindSamLines.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/FindSamLines.%j.err    # Standard error log

#!/bin/bash
Path="/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/"
# Read A.txt line by line
while IFS=$'\t' read -r _ _ _ _ _ _ line_a
do
    # Find lines in B.sam that start with the third column of A.txt and save them to Output
    grep "^$line_a" /scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode/3_bif3_2_BarcodeFixed.sam >> "$Path"Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.sam
done < "$Path"Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.intersect


while IFS=$'\t' read -r _ _ _ _ _ _ line_a
do
    # Find lines in B.sam that start with the third column of A.txt and save them to Output
    grep "^$line_a" /scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode/3_bif3_2_BarcodeFixed.sam >> "$Path"Bif3_Re2_ToBif3Ref_Added500bp.sam
done < "$Path"Bif3_Re2_ToBif3Ref_Added500bp.intersect
