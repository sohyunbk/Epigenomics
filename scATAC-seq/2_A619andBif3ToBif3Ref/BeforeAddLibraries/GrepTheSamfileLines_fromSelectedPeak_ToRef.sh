#!/bin/bash
#SBATCH --job-name=FindSamLines        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/FindSamLines.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/FindSamLines.%j.err    # Standard error log

#!/bin/bash
Path="/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/"
# Read A.txt line by line
#!/bin/bash

cd $Path
# Input files
input_file_a="Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.intersect"
input_file_b="/scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode/3_bif3_2_BarcodeFixed.sam"

# Output file
output_file="Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.sam"

# Extract unique values from the 6th column of A.txt and create a pattern file
awk '{print $6}' "$input_file_a" | sort | uniq > patterns2.txt

# Define the number of parallel processes
num_processes=30

# Search for matching lines in B.sam using multiple parallel processes
xargs -P "$num_processes" -a patterns.txt -I {} grep -w {} "$input_file_b" >> "$output_file"
xargs -P 3 -a patterns2.txt -I {} grep -w {} Bif3_Re2_ToA619Ref_ZmWUS1PromoterPeak.sam

# Cleanup: Remove the temporary pattern file
rm "$Path"patterns.txt
#    grep "^$line_a" /scratch/sb14489/3.scATAC/4.Bif3Ref/4.Bam_FixingBarcode/3_bif3_2_BarcodeFixed.sam >> "$Path"Bif3_Re2_ToBif3Ref_Added500bp.sam
#done < "$Path"Bif3_Re2_ToBif3Ref_Added500bp.intersect
