#!/bin/bash
#SBATCH --job-name=Try2        # Job name
#SBATCH --partition=batch        # Partition (queue) name
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=80gb                   # Job memory request
#SBATCH --time=60:00:00               # Time limit hrs:min:sec
#SBATCH --output=/scratch/sb14489/0.log/Try2.%j.out   # Standard output log
#SBATCH --error=/scratch/sb14489/0.log/Try2.%j.err    # Standard error log

module load SAMtools/1.16.1-GCC-11.3.0
#!/bin/bash
# Read A.txt line by line
#!/bin/bash


# Input files
#samtools view -@ 30 /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.bam "chr2:3679079-3679413"
samtools index -@  30 /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.bam
samtools view -@ 30 /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.bam \
 "chr2:3679079-3679413" > /scratch/sb14489/3.scATAC/4.Bif3Ref/3.SortedBam/3_bif3_2_Markingpcr.chr2_3679079_3679413
