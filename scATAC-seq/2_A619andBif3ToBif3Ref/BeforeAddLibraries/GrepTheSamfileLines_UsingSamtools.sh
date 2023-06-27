#!/bin/bash
#SBATCH --job-name=Try        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Try.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Try.%j.err    # Standard error log

module load SAMtools/1.16.1-GCC-11.3.0
#!/bin/bash
Path="/scratch/sb14489/3.scATAC/4.Bif3Ref/6.Compare_Reads_inTwoRegions/"
# Read A.txt line by line
#!/bin/bash

cd $Path
# Input files

samtools view -@ 30  --output-fmt SAM -h /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.bam  \
 -L ZmWUS1PromoterRegions.bed \
 -o Bif3Re2_ToBif3Ref.sam
